% VALIDATE_THBI_BEFORE_SYNTHETIC -- MATLAB implementation/native execution audit.
% No production inversion. Two 12-step native chains only exercise wiring,
% checkpoint continuation and postprocessing; they cannot assess recovery.

%% =================== SETTINGS, PATHS AND PRESERVED AUDIT OUTPUT ===================
repo=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(repo,'THBI'));
a0_STARTUP_BAYES(); par=bayes_inv_parms('reference100_vti');
par.res.show_figures=false;
par.res.dir=fullfile(par.res.root,['WORKFLOW_VALIDATION_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(par.res.dir,'dir'),'thbi:audit','Audit output already exists.');
mkdir(par.res.dir); audit_dir=par.res.dir;
diary(fullfile(audit_dir,'validation.log'));
fprintf('AUDIT_OUTPUT=%s\nMATLAB %s | %s\n',audit_dir,version,computer);
copyfile([mfilename('fullpath') '.m'],fullfile(audit_dir,'validation_launcher.m'));
validate_configuration(par);
report=struct('passed',false,'runtime',version,'started',datestr(now,30));

%% =================== MATLAB PRIOR, RJ, FULL DR, CACHE AND RESTART CHECKS ===================
test_thbi_reference100(); test_thbi_vti();
test_thbi_reference100('reference100_vti');
test_thbi_reference_chain('reference100_vti');
test_thbi_rj(); test_thbi_layered_dr();
report.discrete_dr=test_thbi_dr_stationarity('layered_path');
report.continuous_prior=test_thbi_prior_invariance(6000,'layered_path');
test_thbi_results(); test_thbi_failures('extended100');
report.algorithm_checks=true;
save(fullfile(audit_dir,'audit_progress.mat'),'par','report','-v7');

%% =================== NATIVE ISOTROPIC LIMIT AND INDEPENDENT VTI PARITY ===================
report.vti=test_thbi_vti_native(par,fullfile(audit_dir,'vti_comparison'));
assert(report.vti.passed);
save(fullfile(audit_dir,'audit_progress.mat'),'par','report','-v7');

%% =================== ACTUAL FULL-GRID SYNTHETIC DATA AND FIGURES ===================
started=tic;
[data,TRUEmodel,forward]=a2_LOAD_DATA(par); par.forc=forward.config;
report.full_grid_seconds=toc(started);
save(fullfile(audit_dir,'inputs.mat'),'par','data','TRUEmodel','forward','-v7');
plot_reference_synthetics(TRUEmodel,data,par); plot_RF_BANDS(forward,par,'true');
fprintf('Full configured synthetic grids generated in %.1f seconds.\n',report.full_grid_seconds);

%% =================== INDEPENDENT NATIVE Ps/Sp WAVEFORM CONTROLS ===================
control=bayes_inv_parms('extended100');
control.forc.cps_bin=par.forc.cps_bin;
report.rf_controls=test_thbi_receiver_functions(control,fullfile(audit_dir,'rf_controls'));
save(fullfile(audit_dir,'audit_progress.mat'),'par','report','-v7');

%% =================== SMALL NATIVE CHAIN FIXTURE: ALL SIX LIKELIHOOD BLOCKS ===================
% Retain the full RF products and representative coordinates from each curve.
% This reduces execution-audit cost without changing any production settings.
smoke=par; smoke_data=data;
for d=1:4
    name=par.inv.datatypes{d}; v=data.(name);
    ix=unique(round(linspace(1,v.N,7)));
    v.x=v.x(ix); v.obs=v.obs(ix); v.clean=v.clean(ix); v.sig=v.sig(ix);
    v.R=v.R(ix,ix); smoke_data.(name)=v;
end
smoke_data=a2_PREPARE_DATA(smoke_data,smoke);
smoke.inv.niter=12; smoke.inv.burnin=4; smoke.inv.nchains=2;
smoke.inv.parallel=false; smoke.inv.print_every=1; smoke.inv.checkpoint_every=4;
smoke.res.prediction_stride=2; smoke.res.nprior=20;
smoke.res.dir=fullfile(audit_dir,'native_chain_smoke'); mkdir(smoke.res.dir);
% Sensitivity/reversibility on the same representative curve grid as the
% execution check, retaining all RF products and coverage through100s.
report.sensitivity=test_thbi_reference_native(smoke,smoke_data,TRUEmodel, ...
    fullfile(audit_dir,'native_sensitivity'));
report.sensitivity_grid='Seven representative coordinates per curve; full RF products.';
prior=a3_BUILD_EMPIRICAL_PRIOR(smoke,smoke.res.nprior);
save(fullfile(smoke.res.dir,'prior.mat'),'prior','-v7');
run_one_chain(smoke,smoke_data,1,5); % stop and resume a real native checkpoint
chains={run_one_chain(smoke,smoke_data,1),run_one_chain(smoke,smoke_data,2)};
for c=1:2
    q=chains{c}; assert(q.last_iter==smoke.inv.niter && q.failures==0);
    assert(all(isfinite(q.logL)) && all(isfinite(q.logP)) && all(q.prediction_saved));
    assert(isequal(q.accepted,q.accepted_stage>0));
    expected=zeros(smoke.inv.niter,3);
    for z=1:3
        used=q.accepted_stage>0;
        expected(:,z)=double(used & q.move==4 & q.zone==z)-double(used & q.move==5 & q.zone==z) ...
            +double(q.accepted_stage==2 & q.dr_move==4 & q.dr_zone==z) ...
            -double(q.accepted_stage==2 & q.dr_move==5 & q.dr_zone==z);
    end
    assert(isequal(diff([q.initial.model.K;q.K],1,1),expected));
    for i=1:smoke.inv.niter, assert(a1_TEST_CONDITIONS(model_from_chain(q,i,smoke),smoke)); end
end
assert(~isequal(chains{1}.initial.model.matrix,chains{2}.initial.model.matrix));

%% =================== c1-c5, SAVED RESPONSES, FRESH LIKELIHOOD AND PLOTS ===================
summary=c1_PROCESS_RESULTS(chains,smoke_data,smoke);
[posterior,retained]=c2_BUILD_POSTERIOR(chains,smoke);
suite=c3_BUILD_MODEL_SUITE(retained,smoke);
final=c4_FINAL_MODEL(posterior,suite,retained,smoke);
assert(posterior.Nstored==16 && summary.bestfit_chain==final.representative_chain && ...
    summary.bestfit_iteration==final.representative_iteration);
fresh=c5_FINAL_FORWARD_MODEL(final,smoke,smoke_data);
c=final.representative_chain; i=final.representative_iteration;
report.final_E2_error=max(abs(fresh.misfit.E2-chains{c}.E2(i,:)')./(1+abs(chains{c}.E2(i,:)')));
report.final_logL_error=abs(fresh.logL-chains{c}.logL(i))/(1+abs(chains{c}.logL(i)));
assert(report.final_E2_error<1e-6 && report.final_logL_error<1e-6);
save(fullfile(smoke.res.dir,'smoke_results.mat'),'smoke','smoke_data','TRUEmodel', ...
    'chains','summary','posterior','suite','final','fresh','-v7');
c2_PLOT_RESULTS(smoke,smoke_data,TRUEmodel,chains,summary,posterior,suite,final);
plot_ANISOTROPY(chains,TRUEmodel,smoke); plot_RF_BANDS(fresh.native,smoke,'final');

%% =================== AUDIT RESULT: EXECUTION, NOT SCIENTIFIC RECOVERY ===================
report.native_chain_steps=[12 12]; report.completed=datestr(now,30);
report.note='Native execution and regression audit only. Short chains do not establish convergence or depth resolution.';
% A successful execution smoke test does not override numerical HVSR failure.
% Preserve this independent result even when the final numerical gate fails.
report.execution_passed=true;
try
    report.hvsr=test_thbi_hvsr_convergence(par,TRUEmodel.zoned_parameters, ...
        fullfile(audit_dir,'hvsr_convergence'));
    report.cutoff=test_thbi_hvsr_cutoff(par,fullfile(audit_dir,'love_cutoff'));
    report.passed=report.hvsr.passed && report.cutoff.passed;
catch problem
    report.hvsr=struct('passed',false,'identifier',problem.identifier,'message',problem.message);
    report.passed=false;
end
save(fullfile(audit_dir,'validation_report.mat'),'report','par','-v7');
fid=fopen(fullfile(audit_dir,'validation_report.json'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report)); fclose(fid);
fprintf('Native execution passed; complete validation passed: %d\n%s\n%s\n',report.passed,report.note,audit_dir);
diary off;
assert(report.passed,'thbi:preflight','Numerical HVSR validation failed; do not launch production.');
