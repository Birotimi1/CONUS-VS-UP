% RUN_ECM1_THIN_THBI_100S_PS_SP -- launch this script yourself in MATLAB.
% Published ECM1 10-m sediment column + AK135 reference mantle to250km.
% HVSR, ellipticity, phase/group velocity, Ps H-kappa and Sp CCP likelihoods.
% New data products/state schema require a NEW run, not an old checkpoint.

%% =================== a0: PATHS AND EXPLICIT TEST SETTINGS ===================
launcher=mfilename('fullpath'); repo=fileparts(fileparts(launcher));
addpath(fullfile(repo,'THBI')); [repo,paths]=a0_STARTUP_BAYES();
par=bayes_inv_parms('reference100_vti');
par.inv.nchains=2; par.inv.niter=15000; par.inv.burnin=8000;
par.inv.parallel=false; par.inv.seed=42000;
% Full-band native forwards can take minutes. Report/preserve every completed
% iteration for this first recovery run; do not wait for 100 costly steps.
par.inv.checkpoint_every=1; par.inv.print_every=1;
par.forc.verbose=true;          % announce each observable before its native call
par.res.nprior=1000; par.res.prediction_stride=10; par.res.show_figures=true;
fprintf('HVSR numerical damping: apsv=%g, ash=%g; native heartbeat=%ds.\n', ...
    par.forc.hvsr.apsv,par.forc.hvsr.ash,par.forc.native_heartbeat_seconds);
run_inversion=true;             % false: generate/plot data and run native preflight only
run_native_rf_controls=true;    % existing independent CPS Ps/Sp waveform controls
validate_configuration(par);
assert(strcmp(par.proposal.dr.mode,'layered_path') && par.proposal.dr.same_zone_weight==0);
assert(strcmp(par.mod.start.method,'prior'),'Use independent dispersed starts.');

%% =================== SOLVER-FREE IMPLEMENTATION CHECKS ===================
fprintf('Implementation checks: the short analytic chains below are not the recovery inversion.\n');
test_thbi_reference100();
test_thbi_vti();
test_thbi_reference100('reference100_vti');
test_thbi_reference_chain('reference100_vti');
test_thbi_rj();
test_thbi_layered_dr();
test_thbi_dr_stationarity('layered_path');

%% =================== FRESH OUTPUT DIRECTORY AND ACTUAL SETTINGS ===================
if ~exist(par.res.root,'dir'), mkdir(par.res.root); end
result_dir=fullfile(par.res.root,['ECM1_19178_AK135_VTI_HK_CCP_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(result_dir,'dir'),'thbi:run','Output directory already exists.');
mkdir(result_dir); par.res.dir=result_dir;
runtime=struct('version',version,'computer',computer,'started',datestr(now,30));
copyfile([launcher '.m'],fullfile(result_dir,'launcher.m'));
save(fullfile(result_dir,'test_settings.mat'),'par','runtime','-v7');
fprintf('\nOutput: %s\n%d chains x %d steps, burn-in %d per chain.\n', ...
    result_dir,par.inv.nchains,par.inv.niter,par.inv.burnin);

%% =================== MANDATORY NATIVE VTI VALIDATION BEFORE DATA OR INVERSION ===================
fprintf('Native validation gates: recovery chains have not started yet.\n');
% Run `sh native/HVf_VTI/build.sh` once after source changes. This checks
% the isotropic limit and independent CPS/HVf nonzero-anisotropy agreement.
vti_report=test_thbi_vti_native(par,fullfile(result_dir,'vti_preflight'));
assert(vti_report.passed,'thbi:preflight','VTI validation did not pass.');
body_report=test_thbi_body_wave_reference(par,fullfile(result_dir,'body_reference'));
assert(body_report.passed,'thbi:preflight','Independent body-wave validation did not pass.');
body_quadrature=test_thbi_body_quadrature(par,body_report,fullfile(result_dir,'body_quadrature'));
assert(body_quadrature.passed,'thbi:preflight','Body-wave quadrature validation did not pass.');
% Fundamental-mode parity alone cannot validate a multimode HVSR sum.
% This gate must pass before synthetic data generation or sampling.
truth_check=a2_BUILD_REFERENCE_TRUTH(par);
hvsr_report=test_thbi_hvsr_convergence(par,truth_check.zoned_parameters, ...
    fullfile(result_dir,'hvsr_convergence'));
cutoff_report=test_thbi_hvsr_cutoff(par,fullfile(result_dir,'love_cutoff'));

%% =================== a2: KNOWN EARTH AND ALL OBSERVED SYNTHETIC PRODUCTS ===================
fprintf('Generating the full synthetic dataset before starting the recovery chains.\n');
[data,TRUEmodel,forward]=a2_LOAD_DATA(par); par.forc=forward.config;
save(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','forward','runtime','-v7');
plot_reference_synthetics(TRUEmodel,data,par);
plot_RF_BANDS(forward,par,'true');
fprintf('%s\n',TRUEmodel.description);
fprintf('Vs below the Moho: %s km/s; half-space at250km: %.4f km/s\n', ...
    mat2str(TRUEmodel.zoned_parameters.mantmparm.VS',5),TRUEmodel.zoned_parameters.halfspaceVS);

%% =================== NATIVE FORWARD, SPLIT/MERGE AND DEPTH-SENSITIVITY CHECKS ===================
fprintf('Native sensitivity preflight: six full forward models; this stage can take hours.\n');
native_report=test_thbi_reference_native(par,data,TRUEmodel,fullfile(result_dir,'native_preflight'));
if run_native_rf_controls
    control_par=bayes_inv_parms('extended100');
    control_par.forc.cps_bin=par.forc.cps_bin; control_par.forc.hvf_bin=par.forc.hvf_bin;
    rf_controls=test_thbi_receiver_functions(control_par,fullfile(result_dir,'native_rf_controls'));
else
    rf_controls=struct('status','not_requested');
end
save(fullfile(result_dir,'preflight_report.mat'),'vti_report','body_report','body_quadrature','hvsr_report','cutoff_report','native_report','rf_controls','-v7');
if ~run_inversion
    fprintf('Synthetic inputs and native preflight saved. No inversion was requested.\n');
    return;
end

%% =================== a3: CONDITIONED PRIOR ===================
fprintf('Native preflight complete. Drawing the conditioned prior before chain initialization.\n');
if par.res.nprior>0
    prior=a3_BUILD_EMPIRICAL_PRIOR(par,par.res.nprior);
    save(fullfile(result_dir,'prior.mat'),'prior','-v7');
end

%% =================== b1-b9: TWO INDEPENDENT RJ CHAINS ===================
% b1 draws independent conditioned-prior starts; b2 samples eleven move families.
% Both DR stages use the ordinary random mixture, including physical-layer
% split/merge moves, with the complete reverse-path acceptance correction.
fprintf('Starting recovery chains: initial forward first, then progress/checkpoints every %d iteration(s).\n',par.inv.print_every);
chains=run_all_chains(par,data);
keep=(par.inv.burnin+1):par.inv.niter;
for c=1:numel(chains)
    chain=chains{c}; assert(chain.last_iter==par.inv.niter);
    assert(all(isfinite(chain.logL)) && all(isfinite(chain.logP)) && all(isfinite(chain.vpvs)));
    assert(isequal(chain.accepted,chain.accepted_stage>0));
    expected=zeros(par.inv.niter,3);
    for z=1:3
        used=chain.accepted_stage>0;
        expected(:,z)=double(used & chain.move==4 & chain.zone==z)-double(used & chain.move==5 & chain.zone==z);
        expected(:,z)=expected(:,z)+double(chain.accepted_stage==2 & chain.dr_move==4 & chain.dr_zone==z) ...
            -double(chain.accepted_stage==2 & chain.dr_move==5 & chain.dr_zone==z);
    end
    assert(isequal(diff([chain.initial.model.K;chain.K],1,1),expected), ...
        'thbi:run','Stored layer counts disagree with accepted RJ paths.');
    for i=1:par.inv.niter
        model=model_from_chain(chain,i,par); assert(a1_TEST_CONDITIONS(model,par));
    end
    assert(all(chain.prediction_saved));
    for d=1:numel(par.inv.datatypes)
        response=chain.predictions.(par.inv.datatypes{d});
        assert(all(isfinite(response(:))),'thbi:run','Missing saved model responses.');
    end
end

%% =================== c1-c4: ALL RETAINED STATES, POSTERIOR AND ACTUAL EARTH ===================
summary=c1_PROCESS_RESULTS(chains,data,par);
[posterior,retained]=c2_BUILD_POSTERIOR(chains,par);
suite_of_models=c3_BUILD_MODEL_SUITE(retained,par);
final_model=c4_FINAL_MODEL(posterior,suite_of_models,retained,par);
assert(summary.bestfit_chain==final_model.representative_chain && ...
    summary.bestfit_iteration==final_model.representative_iteration);
assert(posterior.Nstored==numel(keep)*numel(chains));
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');
save(fullfile(result_dir,'postprocessing.mat'),'posterior','suite_of_models','final_model','-v7.3');

%% =================== c5: FRESH FORWARD RECOMPUTATION OF RETAINED MODEL ===================
final_forward=c5_FINAL_FORWARD_MODEL(final_model,par,data);
c=final_model.representative_chain; i=final_model.representative_iteration;
error_E2=max(abs(final_forward.misfit.E2-chains{c}.E2(i,:)')./(1+abs(chains{c}.E2(i,:)')));
error_logL=abs(final_forward.logL-chains{c}.logL(i))/(1+abs(chains{c}.logL(i)));
assert(error_E2<1e-6 && error_logL<1e-6,'thbi:run','Fresh forward differs from saved likelihood.');
plot_RF_BANDS(final_forward.native,par,'final');
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary,posterior,suite_of_models,final_model);
plot_ANISOTROPY(chains,TRUEmodel,par);

%% =================== SAVE DIAGNOSTICS; FINISHING DOES NOT ESTABLISH RECOVERY ===================
test_report=struct('status','Execution checks passed; inspect convergence, RF products and depth resolution.', ...
    'vti_preflight',vti_report,'native_preflight',native_report,'rf_controls',rf_controls,'final_E2_relative_error',error_E2, ...
    'final_logL_relative_error',error_logL,'rank_diagnostics',summary.rank_diagnostics, ...
    'forward_consistency_question','? OPEN -- dispersion-only extension/shared-Earth comparison deferred by user');
save(fullfile(result_dir,'test_report.mat'),'test_report','-v7');
fid=fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',result_dir); fclose(fid);
fprintf('\n%s\nResults and figures: %s\n',test_report.status,result_dir);
% After interruption with identical code/settings: resume_thbi(result_dir).
% Never use this updated schema to resume an older waveform-likelihood run.
