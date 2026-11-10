function [table,velocity,report] = refine_cps_modal_reference(table,velocity,native,frequency,wave,matrix,cfg,folder)
% Independent CPS search of intervals containing unmatched HVf modes.
% A global sign-change grid can skip a narrow pair. CPS's documented
% -cmin/-cmax bounds let its own secular equation resolve that interval.
% Native values define brackets only; no native root or energy is inserted.
previous=pwd; restore=onCleanup(@() cd(previous)); %#ok<NASGU>
native=native(:); velocity=velocity(:); report=struct([]);
assert(all(diff(native)>0),'thbi:preflight','Native roots must be ordered.');
if isempty(velocity)
    missing=(1:numel(native))';
else
    distance=abs(native-velocity')./native;
    missing=find(min(distance,[],2)>1e-6);
end
assert(~isempty(missing),'thbi:preflight','Unresolved CPS count discrepancy without isolated missing roots.');
first=[1;find(diff(missing)>1)+1]; last=[first(2:end)-1;numel(missing)];
for g=1:numel(first)
    indices=missing(first(g):last(g)); lo=indices(1); hi=indices(end);
    left=native(lo)*(1-1e-3); right=native(hi)*(1+1e-3);
    if lo>1, left=mean(native(lo-1:lo)); end
    if hi<numel(native), right=mean(native(hi:hi+1)); end
    target=sum(native>left & native<right);
    suffix='ray'; ascii='TREGN.ASC'; eigen='tregn96';
    if wave=='L', suffix='lov'; ascii='TLEGN.ASC'; eigen='tlegn96'; end
    matched=false;
    for factor=[1 4 16]
        work=fullfile(folder,sprintf('interval_%03d_factor_%02d',g,factor));
        if ~exist(work,'dir'), mkdir(work); end
        cd(work); write_elastic_cps_model('start.mod',matrix);
        fid=fopen('periods.txt','w'); assert(fid>=0,'fiveobs:io','Cannot write CPS interval period.');
        fprintf(fid,'%.15g\n',1/frequency); fclose(fid);
        commands={sprintf('%s -M start.mod -PARR periods.txt -%s -NMOD %d -FACR %d -FACL %d', ...
            quote_shell_path(fullfile(cfg.cps_bin,'tprep96')),wave,max(20,target+4),factor,factor), ...
            sprintf('%s -cmin %.12g -cmax %.12g',quote_shell_path(fullfile(cfg.cps_bin,'tdisp96')),left,right), ...
            quote_shell_path(fullfile(cfg.cps_bin,eigen)), ...
            sprintf('%s -%s -C -PER -S -ASC',quote_shell_path(fullfile(cfg.cps_bin,'tdpegn96')),wave)};
        for j=1:numel(commands)
            [status,message]=system([commands{j} ' > stage_' num2str(j) '.log 2>&1']);
            assert(status==0,'thbi:preflight','CPS interval stage%d failed: %s',j,message);
        end
        binary=read_cps_tdisp(['tdisp96.' suffix]);
        assert(numel(binary.period)==1 && abs(binary.period-1/frequency)<1e-5);
        c=binary.velocity{1};
        if isempty(c), continue; end
        rows=dlmread(ascii,'',1,0); [~,order]=sort(rows(:,1)); rows=rows(order,:);
        assert(size(rows,1)==numel(c),'thbi:preflight','CPS interval modal schema mismatch.');
        keep=c>left & c<right; c=c(keep); rows=rows(keep,:);
        if numel(c)==target
            [c,order]=sort(c); rows=rows(order,:);
            expected=native(native>left & native<right);
            if all(abs(c-expected)./expected<3e-4), matched=true; break; end
        end
    end
    assert(matched,'thbi:preflight', ...
        'Independent CPS could not resolve %s interval %.9g to %.9g at %.8g Hz.',wave,left,right,frequency);
    keep=velocity<=left | velocity>=right;
    table=[table(keep,:);rows]; velocity=[velocity(keep);c]; %#ok<AGROW>
    [velocity,order]=sort(velocity); table=table(order,:);
    report(g)=struct('frequency_Hz',frequency,'wave',wave,'velocity_bounds',[left right], ...
        'factor',factor,'native_count',target,'cps_velocity',c,'commands',{commands},'folder',work); %#ok<AGROW>
    cd(previous);
end
assert(numel(velocity)==numel(native),'thbi:preflight','CPS refined mode count still differs from HVf.');
end
