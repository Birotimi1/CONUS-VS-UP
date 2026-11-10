function report = test_thbi_hvsr_convergence(par,model,folder)
% Numerical HVSR convergence at low/high frequencies and across the thin cap peak.
% Same damped elastic model; vary body quadrature and surface-mode truncation.
% This is not an independent body-wave Green-function implementation.
if ~exist(folder,'dir'), mkdir(folder); end
previous=pwd; restore=onCleanup(@() cd(previous)); %#ok<NASGU>
matrix=refine_layer_model(model.matrix,par.forc.max_layer_thickness_km);
frequencies=[.01;.1;.5;1;5;10;15;18;19;20;22;25];
settings=[10 1000;10 2000;20 2000;40 4000]; % [Rayleigh/Love modes, body nodes]
values=zeros(numel(frequencies),size(settings,1));
for j=1:size(settings,1)
    work=fullfile(folder,sprintf('modes_%d_nodes_%d',settings(j,:))); mkdir(work); cd(work);
    fid=fopen('model.txt','w'); assert(fid>=0);
    fprintf(fid,'%d\n',size(matrix,1)); exchange=matrix; exchange(:,1:4)=1000*exchange(:,1:4);
    fprintf(fid,[repmat('%.12g ',1,7) '\n'],exchange'); fclose(fid);
    fid=fopen('frequencies.txt','w'); assert(fid>=0);
    fprintf(fid,'%.15g\n',frequencies); fclose(fid);
    command=sprintf('%s -vti -f model.txt -ff frequencies.txt -nmr %d -nml %d -nks %d -prec 1e-4 -apsv 1e-3 -ash 1e-3 -rep -hv', ...
        quote_shell_path(fullfile(par.forc.hvf_vti_bin,'HVf_VTI')),settings(j,1),settings(j,1),settings(j,2));
    [status,message]=system([command ' > HV.dat 2> HVf.log']);
    assert(status==0,'thbi:preflight','HVSR convergence calculation failed: %s',message);
    a=sscanf(fileread('HV.dat'),'%f'); assert(numel(a)==2*numel(frequencies));
    a=reshape(a,2,[])'; assert(all(isfinite(a(:,2)) & a(:,2)>0));
    assert(max(abs(a(:,1)-frequencies))<1e-5);
    values(:,j)=a(:,2); cd(previous);
end
report=struct('frequency_Hz',frequencies,'settings',settings,'hvsr',values, ...
    'relative_change',abs(values(:,1:3)-values(:,4))./abs(values(:,4)), ...
    'note','Convergence for this fixture only; no independent VTI body-wave amplitude comparison.');
report.maximum_relative_change=max(report.relative_change,[],1);
report.passed=all(report.maximum_relative_change<.02);
save(fullfile(folder,'convergence.mat'),'report','model','-v7');
fprintf('HVSR maximum relative differences vs 40 modes/4000 nodes: %s\n',mat2str(report.maximum_relative_change,5));
assert(report.passed,'thbi:preflight','HVSR truncation/quadrature changes exceed 2%%. Inspect convergence.mat.');
end
