function report = test_thbi_hvsr_convergence(par,model,folder,frequencies,reuse_saved)
% Complete modes, root/quadrature convergence and independent CPS surface
% Green amplitudes. Body-wave physics has one implementation; doubling its
% quadrature is a numerical check, not independent physics validation.
if nargin<4, frequencies=[.01;.1;.5;1;5;10;15;18;19;20;22;25]; end
if nargin<5, reuse_saved=false; end
frequencies=frequencies(:); assert(all(diff(frequencies)>0));
if ~exist(folder,'dir'), mkdir(folder); end
previous=pwd; restore=onCleanup(@() cd(previous)); %#ok<NASGU>
signature=check_native_vti_build(par.forc);
matrix=refine_layer_model(model.matrix,par.forc.max_layer_thickness_km);
base=par.forc.hvsr; assert(base.require_complete,'thbi:preflight','Require complete modes.');
fine=base; fine.nmr=2*base.nmr; fine.nml=2*base.nml;
fine.nks=2*base.nks; fine.prec=base.prec/10; fine.dphase=base.dphase/2;
settings={base,fine}; values=zeros(numel(frequencies),2); modes=cell(1,2);
report=struct('passed',false,'frequency_Hz',frequencies,'settings',{settings}, ...
    'native_signature',signature,'note','Independent CPS surface amplitudes; body quadrature convergence only.');
if reuse_saved
    saved=load(fullfile(folder,'convergence.mat'),'report','model');
    assert(isequaln(saved.model,model) && isequaln(saved.report.settings,settings) && ...
        isequaln(saved.report.native_signature,signature) && ...
        isequaln(saved.report.frequency_Hz,frequencies),'thbi:preflight','Saved native fixture changed.');
    values=saved.report.hvsr; modes=saved.report.mode_counts;
end

%% =================== CONFIGURED AND TIGHTER NATIVE SEARCHES ===================
for j=1:2*double(~reuse_saved)
    cfg=settings{j}; work=fullfile(folder,sprintf('resolution_%d',j)); mkdir(work); cd(work);
    fid=fopen('model.txt','w'); assert(fid>=0);
    fprintf(fid,'%d\n',size(matrix,1)); exchange=matrix; exchange(:,1:4)=1000*exchange(:,1:4);
    fprintf(fid,[repmat('%.12g ',1,7) '\n'],exchange'); fclose(fid);
    fid=fopen('frequencies.txt','w'); assert(fid>=0);
    fprintf(fid,'%.15g\n',frequencies); fclose(fid);
    command=sprintf('%s -vti -complete -f model.txt -ff frequencies.txt -nmr %d -nml %d -nks %d -prec %.15g -dphase %.15g -apsv %.15g -ash %.15g -ph -gr -rep -hv', ...
        quote_shell_path(fullfile(par.forc.hvf_vti_bin,'HVf_VTI')), ...
        cfg.nmr,cfg.nml,cfg.nks,cfg.prec,cfg.dphase,cfg.apsv,cfg.ash);
    [status,message]=system([command ' > HV.dat 2> HVf.log']);
    assert(status==0,'thbi:preflight','HVSR convergence calculation failed: %s. See %s',message,work);
    a=sscanf(fileread('HV.dat'),'%f'); assert(numel(a)==2*numel(frequencies));
    a=reshape(a,2,[])'; assert(all(isfinite(a(:,2)) & a(:,2)>0));
    assert(max(abs(a(:,1)-frequencies))<1e-5);
    values(:,j)=a(:,2); modes{j}=dlmread('modal_counts.dat','',1,0);
    cd(previous);
end
report.hvsr=values; report.mode_counts=modes;
report.relative_change=abs(values(:,1)-values(:,2))./abs(values(:,2));
report.maximum_relative_change=max(report.relative_change);
save(fullfile(folder,'convergence.mat'),'report','model','-v7');
assert(isequal(modes{1}(:,2:3),modes{2}(:,2:3)),'thbi:preflight','Mode counts change under refinement.');
assert(report.maximum_relative_change<.02,'thbi:preflight','HVSR root/quadrature refinement exceeds 2%%.');

%% =================== INDEPENDENT CPS ROOTS AND SURFACE-WAVE AMPLITUDES ===================
work=fullfile(folder,'independent_cps'); mkdir(work); cd(work);
write_elastic_cps_model('start.mod',matrix);
fid=fopen('periods.txt','w'); assert(fid>=0); fprintf(fid,'%.15g\n',1./frequencies); fclose(fid);
cap=2000; % official CPS330 tio.f/getegn arrays have MAXMOD=2000
commands={sprintf('%s -M start.mod -PARR periods.txt -R -L -NMOD %d -FACR 40 -FACL 40', ...
    quote_shell_path(fullfile(par.forc.cps_bin,'tprep96')),cap), ...
    quote_shell_path(fullfile(par.forc.cps_bin,'tdisp96')), ...
    quote_shell_path(fullfile(par.forc.cps_bin,'tregn96')), ...
    [quote_shell_path(fullfile(par.forc.cps_bin,'tdpegn96')) ' -R -C -PER -S -ASC'], ...
    quote_shell_path(fullfile(par.forc.cps_bin,'tlegn96')), ...
    [quote_shell_path(fullfile(par.forc.cps_bin,'tdpegn96')) ' -L -C -PER -S -ASC']};
for j=1:numel(commands)*double(~reuse_saved)
    [status,message]=system([commands{j} ' > stage_' num2str(j) '.log 2>&1']);
    assert(status==0,'thbi:preflight','CPS multimode stage%d failed: %s',j,message);
end
R=dlmread('TREGN.ASC','',1,0); L=dlmread('TLEGN.ASC','',1,0);
assert(size(R,2)==9 && size(L,2)==8,'thbi:preflight','Unexpected CPS modal table.');
nr=numel(frequencies); cps_green=zeros(nr,3); cps_counts=zeros(nr,2);
phase_error=zeros(nr,2); group_error=zeros(nr,2);
excluded_counts=zeros(nr,2);
boundary_failures=zeros(nr,2);
for w=1:2
    wave='RL'; table=R; if w==2, table=L; end
    suffix={'ray','lov'}; binary=read_cps_tdisp(['tdisp96.' suffix{w}]);
    [slow,valid]=read_modes(fullfile(folder,'resolution_1',[wave(w) 'ph.dat']),nr);
    [group,~]=read_modes(fullfile(folder,'resolution_1',[wave(w) 'gr.dat']),nr);
    for k=1:nr
        a=table(abs(table(:,4)-frequencies(k))<6e-5*max(1,frequencies(k)),:);
        [~,order]=sort(a(:,1)); a=a(order,:);
        ix=find(abs(1./binary.period-frequencies(k))<1e-6*max(1,frequencies(k)));
        assert(isscalar(ix),'thbi:preflight','CPS binary period mismatch.');
        c=binary.velocity{ix};
        assert(size(a,1)==numel(c) && numel(c)<cap && isequal(a(:,1),(0:numel(c)-1)'), ...
            'thbi:preflight','CPS ASCII/binary modes differ or reached capacity.');
        % This discrete surface sum includes trapped modes only. CPS can emit
        % boundary roots and radiating rows above its half-space shear speed.
        % Use its binary tensor/velocities, never rounded ASCII or HVf counts.
        trapped=c<binary.halfspace_shear(w)*(1-1e-9);
        excluded_counts(k,w)=sum(~trapped); a=a(trapped,:); c=c(trapped);
        % At long periods, independently screen CPS TI repeated-root artefacts.
        % High-frequency roots are checked individually against HVf below;
        % applying expm/QR to every deep high-frequency mode is unnecessary.
        if w==1 && frequencies(k)<=.1
            bc=arrayfun(@(v) cps_rayleigh_boundary_check(model.matrix,1/frequencies(k),v),c);
            valid_bc=bc<1e-3; boundary_failures(k,w)=sum(~valid_bc);
            a=a(valid_bc,:); c=c(valid_bc);
        end
        cps_counts(k,w)=numel(c);
        used=find(valid(k,:));
        assert(numel(used)==size(a,1) && numel(used)<cap,'thbi:preflight', ...
            'CPS/HVf %s count mismatch at %.8g Hz (%d/%d).',wave(w),frequencies(k),size(a,1),numel(used));
        if isempty(used), continue; end % physically absent Love modes
        phase=1./slow(k,used)'/1000; velocity=1./group(k,used)'/1000;
        phase_error(k,w)=max(abs(phase-c)./c);
        group_error(k,w)=sum(abs(velocity-a(:,6)).*abs(a(:,7)))/sum(abs(a(:,6).*a(:,7)));
        % CPS ARE/ALE=1/(2*c*U*I0). KGS->SI contributes1e-12.
        % HVf prints azimuth-averaged per-horizontal Green components.
        if w==1
            cps_green(k,1)=-.25*sum(a(:,7).*a(:,9).^2)*1e-12;
            cps_green(k,3)=-.5*sum(a(:,7))*1e-12;
        else
            cps_green(k,2)=-.25*sum(a(:,7))*1e-12;
        end
    end
end
raw=fileread(fullfile(folder,'resolution_1','G123.dat'));
raw=regexprep(raw,'[^\n]*FREQ\.[^\n]*\n','');
a=sscanf(raw,'%f'); assert(mod(numel(a),4*nr)==0);
a=reshape(a,4,nr,[]); native_green=squeeze(sum(a(2:4,:,:),3))';
green_error=abs(native_green-cps_green)./max(abs(cps_green),1e-30);
report.cps_counts=cps_counts; report.phase_relative_error=phase_error;
report.cps_boundary_or_radiating_counts=excluded_counts;
report.cps_failed_surface_condition_counts=boundary_failures;
report.energy_weighted_group_error=group_error;
report.native_surface_green=native_green; report.cps_surface_green=cps_green;
report.surface_green_relative_error=green_error;
cd(previous); save(fullfile(folder,'convergence.mat'),'report','model','-v7');
assert(all(phase_error(:)<3e-4) && all(group_error(:)<.015) && all(green_error(:)<.01), ...
    'thbi:preflight','CPS/HVf complete surface-mode comparison failed; inspect convergence.mat.');

%% =================== NEGATIVE CONTROL: SILENT MODE TRUNCATION IS FORBIDDEN ===================
work=fullfile(folder,'insufficient_cap'); mkdir(work); cd(work);
copyfile(fullfile(folder,'resolution_1','model.txt'),'model.txt');
copyfile(fullfile(folder,'resolution_1','frequencies.txt'),'frequencies.txt');
command=sprintf('%s -vti -complete -f model.txt -ff frequencies.txt -nmr 1 -nml 1 -nks 0 -prec 1e-8 -rep -hv', ...
    quote_shell_path(fullfile(par.forc.hvf_vti_bin,'HVf_VTI')));
[status,~]=system([command ' > HV.dat 2> HVf.log']);
assert(status~=0 && isempty(strtrim(fileread('HV.dat'))) && ...
    contains(fileread('HVf.log'),'mode cap reached'),'thbi:preflight','Truncated H/V escaped the native guard.');
cd(previous); report.passed=true;
save(fullfile(folder,'convergence.mat'),'report','model','-v7');
fprintf('PASS HVSR: refinement %.3g, CPS surface Green %.3g max relative error; complete-mode counts agree.\n', ...
    report.maximum_relative_change,max(green_error(:)));
end

function [slow,valid]=read_modes(filename,n)
token=regexp(strtrim(fileread(filename)),'\s+','split');
nx=str2double(token{1}); nm=str2double(token{2});
assert(nx==n && numel(token)==2+2*nx*nm,'thbi:preflight','Malformed modal file: %s',filename);
slow=reshape(str2double(token(3:2+nx*nm)),nx,nm);
valid=reshape(strcmpi(token(3+nx*nm:end),'T'),nx,nm);
assert(all(isfinite(slow(valid)) & slow(valid)>0),'thbi:preflight','Nonphysical modal slowness.');
end
