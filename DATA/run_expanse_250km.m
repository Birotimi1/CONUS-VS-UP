% RUN_EXPANSE_250KM -- separate prepared data, independent chain jobs and plots.
% Native calculations run on compute nodes. Settings mirror the local 250-km
% VTI launcher: two chains, 15k steps, 8k burn-in, damped HVSR 0.1-25Hz.
% CONUS_STAGE, CONUS_RUN_DIR and CONUS_RUNTIME_SCRATCH are set by the job script.

%% =================== PROJECT, JOB AND DURABLE OUTPUT PATHS ===================
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();
stage=getenv('CONUS_STAGE'); result_dir=getenv('CONUS_RUN_DIR');
scratch=getenv('CONUS_RUNTIME_SCRATCH');
assert(any(strcmp(stage,{'prepare','benchmark','chains','postprocess'})),'thbi:batch','Unknown batch stage.');
assert(~isempty(result_dir) && result_dir(1)=='/' && exist(scratch,'dir')==7, ...
    'thbi:batch','Absolute persistent run and existing node-local scratch paths are required.');
if ~exist(result_dir,'dir'), mkdir(result_dir); end
log_dir=fullfile(result_dir,'logs'); if ~exist(log_dir,'dir'), mkdir(log_dir); end
job_id=getenv('SLURM_JOB_ID'); if isempty(job_id), job_id='manual'; end
chain_id=str2double(getenv('CONUS_CHAIN_ID'));
diary(fullfile(log_dir,sprintf('matlab_%s_%s_%g.log',stage,job_id,chain_id)));
fprintf('MATLAB %s | %s | stage=%s | %s\n',version,computer,stage,datestr(now,31));
prepared_file=fullfile(result_dir,'prepared.mat');

if strcmp(stage,'prepare')
    %% =================== CANONICAL CONFIGURATION AND IMMUTABLE IDENTITY ===================
    par=bayes_inv_parms('reference100_vti');
    par.inv.nchains=2; par.inv.niter=15000; par.inv.burnin=8000;
    par.inv.parallel=false; par.inv.seed=42000;
    par.inv.checkpoint_every=1; par.inv.print_every=1;
    par.forc.verbose=true; par.res.nprior=1000;
    par.res.prediction_stride=10; par.res.show_figures=false; par.res.dir=result_dir;
    par.res.scratch_root='/tmp'; % portable fallback; each job supplies its actual scratch
    validate_configuration(par);
    signature=check_native_vti_build(par.forc);
    identity=struct('par',par,'source',{source_snapshot()},'native_signature',signature);
    identity_file=fullfile(result_dir,'preparation_identity.mat');
    if exist(identity_file,'file')
        previous=load(identity_file,'identity');
        assert(isequaln(identity,previous.identity),'thbi:batch','Preparation code/configuration/binaries changed. Choose a new RUN_ID.');
    else
        save(identity_file,'identity','-v7');
    end
    if exist(prepared_file,'file')
        fprintf('Prepared bundle already complete and identity verified; no preflights repeated.\n'); diary off; return;
    end
    fprintf('HVSR damping apsv=%g ash=%g; native heartbeat=%ds, timeout=%ds.\n', ...
        par.forc.hvsr.apsv,par.forc.hvsr.ash,par.forc.native_heartbeat_seconds,par.forc.native_timeout_seconds);
    workpar=par; workpar.res.scratch_root=scratch;

    %% =================== ALGORITHM CHECKS (NO SCIENTIFIC INVERSION) ===================
    test_thbi_reference100(); test_thbi_vti(); test_thbi_reference100('reference100_vti');
    test_thbi_reference_chain('reference100_vti'); test_thbi_rj();
    test_thbi_layered_dr(); test_thbi_dr_stationarity('layered_path');

    %% =================== NATIVE VALIDATION, CACHED ONLY WITH THIS IDENTITY ===================
    gates_file=fullfile(result_dir,'native_gates.mat');
    if exist(gates_file,'file')
        load(gates_file,'vti_report','hvsr_report','cutoff_report');
    else
        vti_report=test_thbi_vti_native(workpar,fullfile(result_dir,'vti_preflight'));
        truth=a2_BUILD_REFERENCE_TRUTH(workpar);
        hvsr_report=test_thbi_hvsr_convergence(workpar,truth.zoned_parameters,fullfile(result_dir,'hvsr_convergence'));
        cutoff_report=test_thbi_hvsr_cutoff(workpar,fullfile(result_dir,'love_cutoff'));
        save(gates_file,'vti_report','hvsr_report','cutoff_report','-v7');
    end
    assert(vti_report.passed && hvsr_report.passed && cutoff_report.passed,'thbi:batch','Native gate failed.');

    %% =================== TRUE EARTH AND SHARED OBSERVATIONS ===================
    inputs_file=fullfile(result_dir,'inputs.mat');
    if exist(inputs_file,'file')
        loaded=load(inputs_file,'data','TRUEmodel','forward');
        data=loaded.data; TRUEmodel=loaded.TRUEmodel; forward=loaded.forward;
    else
        [data,TRUEmodel,forward]=a2_LOAD_DATA(workpar);
        runtime=struct('version',version,'computer',computer,'started',datestr(now,30));
        par.forc=forward.config;
        save(inputs_file,'par','data','TRUEmodel','forward','runtime','-v7');
    end
    par.forc=forward.config; workpar.forc=par.forc;
    plot_reference_synthetics(TRUEmodel,data,workpar); plot_RF_BANDS(forward,workpar,'true');

    %% =================== SENSITIVITY AND RF CONTROLS, ONCE PER PREPARED RUN ===================
    checks_file=fullfile(result_dir,'preflight_report.mat');
    if exist(checks_file,'file')
        load(checks_file,'native_report','rf_controls');
    else
        native_report=test_thbi_reference_native(workpar,data,TRUEmodel,fullfile(result_dir,'native_preflight'));
        control=bayes_inv_parms('extended100'); control.forc.cps_bin=par.forc.cps_bin;
        control.forc.hvf_bin=par.forc.hvf_isotropic_bin;
        rf_controls=test_thbi_receiver_functions(control,fullfile(result_dir,'native_rf_controls'));
        save(checks_file,'vti_report','hvsr_report','cutoff_report','native_report','rf_controls','-v7');
    end
    if par.res.nprior>0 && ~exist(fullfile(result_dir,'prior.mat'),'file')
        prior=a3_BUILD_EMPIRICAL_PRIOR(par,par.res.nprior);
        save(fullfile(result_dir,'prior.mat'),'prior','-v7');
    end
    prepared=struct('source',{source_snapshot()},'native_signature',check_native_vti_build(par.forc), ...
        'created',datestr(now,30),'par',par);
    save_checkpoint(prepared_file,prepared); % atomic .mat writer, variable name checkpoint
    fprintf('PREPARED: %s. Next submit benchmark to measure both dispersed initial models.\n',result_dir);
else
    %% =================== VERIFY THE SHARED BUNDLE BEFORE EACH JOB ===================
    assert(exist(prepared_file,'file')==2,'thbi:batch','Run prepare successfully first.');
    loaded=load(prepared_file,'checkpoint'); prepared=loaded.checkpoint; par=prepared.par;
    assert(isequaln(prepared.source,source_snapshot()),'thbi:batch','Source changed after preparation; use a new RUN_ID.');
    assert(isequaln(prepared.native_signature,check_native_vti_build(par.forc)), ...
        'thbi:batch','Native binaries changed after preparation; use a new RUN_ID.');
    load(fullfile(result_dir,'inputs.mat'),'data','TRUEmodel','runtime');
    if any(strcmp(stage,{'benchmark','chains'}))
        %% =================== ONE INDEPENDENT CHAIN PER ARRAY TASK ===================
        assert(isscalar(chain_id) && isfinite(chain_id) && chain_id==fix(chain_id) && ...
            chain_id>=1 && chain_id<=par.inv.nchains,'thbi:batch','Invalid array chain ID.');
        seconds=str2double(getenv('CONUS_CHAIN_SECONDS'));
        assert(isfinite(seconds) && seconds>0,'thbi:batch','Set CONUS_CHAIN_SECONDS.');
        operation=struct('max_wall_seconds',seconds,'scratch_root',scratch);
        stop_iter=par.inv.niter; if strcmp(stage,'benchmark'), stop_iter=0; end
        started=tic;
        chain=run_one_chain(par,data,chain_id,stop_iter,[],operation);
        fprintf('Chain %d returned at %d/%d after %.1fs; native calls=%d.\n', ...
            chain_id,chain.last_iter,par.inv.niter,toc(started),chain.nforward);
        if chain.last_iter<par.inv.niter
            fprintf('This chain is partial. Submit chains again with the same RUN_ID to continue.\n');
        end
    else
        %% =================== POSTPROCESS ONLY COMPLETE INDEPENDENT CHAINS ===================
        chains=cell(par.inv.nchains,1);
        for c=1:par.inv.nchains
            filename=fullfile(result_dir,sprintf('chain_%02d',c),'chain.mat');
            assert(exist(filename,'file')==2,'thbi:batch','Chain %d has not finished; resume it first.',c);
            one=load(filename,'chain'); chains{c}=one.chain;
            assert(chains{c}.last_iter==par.inv.niter && all(chains{c}.prediction_saved), ...
                'thbi:batch','Incomplete chain %d.',c);
        end
        par.res.scratch_root=scratch; % operational only, no checkpoint is rewritten
        summary=c1_PROCESS_RESULTS(chains,data,par);
        [posterior,retained]=c2_BUILD_POSTERIOR(chains,par);
        suite_of_models=c3_BUILD_MODEL_SUITE(retained,par);
        final_model=c4_FINAL_MODEL(posterior,suite_of_models,retained,par);
        final_forward=c5_FINAL_FORWARD_MODEL(final_model,par,data);
        c=final_model.representative_chain; i=final_model.representative_iteration;
        error_E2=max(abs(final_forward.misfit.E2-chains{c}.E2(i,:)')./(1+abs(chains{c}.E2(i,:)')));
        error_logL=abs(final_forward.logL-chains{c}.logL(i))/(1+abs(chains{c}.logL(i)));
        assert(error_E2<1e-6 && error_logL<1e-6,'thbi:batch','Fresh final likelihood disagrees.');
        save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');
        save(fullfile(result_dir,'postprocessing.mat'),'posterior','suite_of_models','final_model','-v7.3');
        c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary,posterior,suite_of_models,final_model);
        plot_ANISOTROPY(chains,TRUEmodel,par); plot_RF_BANDS(final_forward.native,par,'final');
        fprintf('Postprocessing complete; inspect convergence and recovery in %s\n',result_dir);
    end
end
diary off;
