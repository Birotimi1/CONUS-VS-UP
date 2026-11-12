function report = test_thbi_native_runtime(folder,run_native)
% Test monitoring/error behavior and, optionally, native numerical damping.
% Tiny two-layer fixture only: no production synthetic inversion is launched.
if nargin<2, run_native=false; end
if ~exist(folder,'dir'), mkdir(folder); end
previous=pwd; cleanup=onCleanup(@() cd(previous)); %#ok<NASGU>
cd(folder);
cfg=struct('native_heartbeat_seconds',1,'native_timeout_seconds',10,'verbose',true);
[status,message]=run_native_command('sleep 1.2; printf ''finished\n''',cfg,'heartbeat fixture');
assert(status==0 && contains(message,'still running') && contains(message,'finished'));
[status,~]=run_native_command('exit 17',cfg,'exit-code fixture'); assert(status==17);
cfg.native_timeout_seconds=1;
marker=fullfile(folder,'must_not_exist.txt');
if exist(marker,'file'), delete(marker); end
caught=false;
try
    run_native_command(['sleep 3; touch ' quote_shell_path(marker)],cfg,'timeout fixture');
catch problem
    caught=strcmp(problem.identifier,'thbi:native_timeout');
end
assert(caught); pause(3.2); assert(~exist(marker,'file'),'Native descendant survived timeout.');
report=struct('monitor_passed',true,'native_executed',run_native);
par=bayes_inv_parms('reference100_vti');
assert(par.forc.hvsr.apsv==1e-3 && par.forc.hvsr.ash==1e-3);
assert(isequal(par.conditions.monotonic,[true false false]));
if run_native
    report.native_signature=check_native_vti_build(par.forc);
    model=[.1 1.2 .5 1.8;0 6.2 3.5 2.7];
    frequencies=[.2;1;3]; damping=[1e-3 1e-2];
    values=zeros(3,2,2);
    cfg.native_timeout_seconds=120; cfg.native_heartbeat_seconds=10;
    for solver=1:2
        matrix=model; flag=''; bin=fullfile(par.forc.hvf_bin,'HVf');
        if solver==2
            matrix=[matrix ones(2,3)]; flag=' -vti'; bin=fullfile(par.forc.hvf_vti_bin,'HVf_VTI');
        end
        for j=1:2
            work=fullfile(folder,sprintf('solver_%d_damping_%d',solver,j)); mkdir(work); cd(work);
            exchange=matrix; exchange(:,1:4)=1000*exchange(:,1:4);
            fid=fopen('model.txt','w'); fprintf(fid,'%d\n',size(matrix,1));
            fprintf(fid,[repmat('%.12g ',1,size(matrix,2)) '\n'],exchange'); fclose(fid);
            fid=fopen('frequencies.txt','w'); fprintf(fid,'%.15g\n',frequencies); fclose(fid);
            command=sprintf('%s%s -f model.txt -ff frequencies.txt -nmr 1 -nml 1 -nks 4000 -prec 1e-4 -apsv %.12g -ash %.12g -rep -hv', ...
                quote_shell_path(bin),flag,damping(j),damping(j));
            [status,~]=run_native_command([command ' > HV.dat 2> HVf.log'],cfg,'native damping fixture');
            assert(status==0); a=reshape(sscanf(fileread('HV.dat'),'%f'),2,[])';
            assert(isequal(size(a),[3 2]) && all(isfinite(a(:,2)) & a(:,2)>0));
            values(:,j,solver)=a(:,2);
        end
    end
    report.frequencies=frequencies; report.damping=damping; report.values=values;
    report.isotropic_limit_error=max(abs(values(:,:,1)-values(:,:,2))./abs(values(:,:,1)),[],1);
    report.damping_response=max(abs(values(:,1,2)-values(:,2,2)));
    assert(all(report.isotropic_limit_error<.01) && report.damping_response>1e-8, ...
        'thbi:test','Damping parity/response failed.');
end
cd(previous); save(fullfile(folder,'runtime_damping_report.mat'),'report','-v7');
fprintf('PASS native progress/timeout/configuration; native damping executed=%d.\n',run_native);
end
