function report = test_thbi_body_quadrature(par,body_report,folder)
% Resolve integration accuracy separately from same-grid physics agreement.
% Only body integrals run here; no surface search, synthetic data or MCMC.
assert(body_report.passed && body_report.native_executed);
signature=check_native_vti_build(par.forc);
assert(isequaln(signature,body_report.native_signature),'thbi:preflight','Native identity changed after body reference.');
if ~exist(folder,'dir'), mkdir(folder); end
previous=pwd; restore=onCleanup(@() cd(previous)); %#ok<NASGU>
report=struct('passed',false,'native_signature',signature,'nks',[par.forc.hvsr.nks 2*par.forc.hvsr.nks], ...
    'apsv',par.forc.hvsr.apsv,'ash',par.forc.hvsr.ash,'cases',{{}});
for k=4:numel(body_report.cases) % truth, xi bounds, LVZ, and independent prior draws
    c=body_report.cases{k}; values=zeros(numel(c.frequency_Hz),3,2);
    fprintf('Body quadrature refinement: %s\n',c.name);
    for j=1:2
        work=fullfile(folder,c.name,sprintf('nks_%d',report.nks(j))); mkdir(work); cd(work);
        exchange=c.native_matrix; exchange(:,1:4)=1000*exchange(:,1:4);
        fid=fopen('model.txt','w'); assert(fid>=0); fprintf(fid,'%d\n',size(exchange,1));
        fprintf(fid,[repmat('%.15g ',1,7) '\n'],exchange'); fclose(fid);
        fid=fopen('frequencies.txt','w'); assert(fid>=0); fprintf(fid,'%.15g\n',c.frequency_Hz); fclose(fid);
        command=sprintf('%s -vti -f model.txt -ff frequencies.txt -nmr 0 -nml 0 -nks %d -apsv %.15g -ash %.15g -rep -hv', ...
            quote_shell_path(fullfile(par.forc.hvf_vti_bin,'HVf_VTI')),report.nks(j),report.apsv,report.ash);
        [status,message]=run_native_command([command ' > HV.dat 2> HVf.log'],par.forc,['body quadrature/' c.name]);
        assert(status==0,'thbi:preflight','Body quadrature native failure: %s',message);
        raw=dlmread('GBW.dat','',1,0); values(:,:,j)=raw(:,2:4);
    end
    relative=abs(values(:,:,1)-values(:,:,2))./abs(values(:,:,2));
    result=struct('name',c.name,'frequency_Hz',c.frequency_Hz,'matrix',c.native_matrix, ...
        'body_integrals',values,'relative_change',relative,'passed',all(relative(:)<.02));
    report.cases{end+1}=result;
    save(fullfile(folder,'body_quadrature_report.mat'),'report','-v7');
    fprintf('  maximum component quadrature change %.6g\n',max(relative(:)));
    assert(result.passed,'thbi:preflight','Body quadrature exceeds 2%% in %s.',c.name);
end
report.passed=true; save(fullfile(folder,'body_quadrature_report.mat'),'report','-v7');
end
