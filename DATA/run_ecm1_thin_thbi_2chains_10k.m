% RUN_ECM1_THIN_THBI_2CHAINS_10K -- run this script yourself in MATLAB.
% Known thin Earth -> five noisy synthetics -> zoned transdimensional HBI.
% Two independent chains, 10,000 iterations each, discard iterations 1:4000.
% All 6,000 remaining states per chain are retained, including repeats.
% This exercises the implementation; convergence/recovery still need review.

%% =================== a0: PROJECT AND SOLVER PATHS ===================
launcher=mfilename('fullpath');
repo=fileparts(fileparts(launcher));
addpath(fullfile(repo,'THBI'));
[repo,paths]=a0_STARTUP_BAYES();

%% =================== EXPLICIT TEST SETTINGS ===================
par=bayes_inv_parms('pilot');
par.inv.profile='2chains_10k_burn4k';
par.inv.nchains=2;
par.inv.niter=10000;
par.inv.burnin=4000;
par.inv.parallel=false;          % two serial chains; no Parallel Toolbox needed
par.inv.seed=42000;              % sampler uses seed + 1009*chain_number
par.inv.checkpoint_every=100;
par.inv.print_every=100;
par.inv.prior_only=false;        % real CPS/HVf likelihoods for all five datasets
par.res.prediction_stride=10;    % 600 cached response draws per chain; no state thinning
par.res.nprior=1000;
par.res.show_figures=true;
run_implementation_checks=true; % short algebraic/toy checks BEFORE native work
validate_configuration(par);
assert(all(par.mod.kmax>par.mod.kmin) && all(par.proposal.probability(4:5)>0), ...
    'thbi:test','This test requires variable layer counts and birth/death moves.');

%% =================== NEW OUTPUT DIRECTORY AND LAUNCHER RECORD ===================
if ~exist(par.res.root,'dir'), mkdir(par.res.root); end
result_dir=fullfile(par.res.root,['ECM1_thin_THBI_' par.inv.profile '_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(result_dir,'dir'),'thbi:test','Result directory already exists; wait a second and rerun.');
mkdir(result_dir); par.res.dir=result_dir;
copyfile([launcher '.m'],fullfile(result_dir,'launcher.m'));
runtime=struct('version',version,'computer',computer,'started',datestr(now,30));
save(fullfile(result_dir,'test_settings.mat'),'par','runtime','-v7');
fprintf('\nOutput: %s\n',result_dir);
fprintf('2 chains x 10000 iterations; burn-in 4000; 6000 retained/chain (12000 total).\n');
fprintf('Chains run serially, with independent seeds and checkpoints every 100 steps.\n');

%% =================== IMPLEMENTATION CHECKS BEFORE EXPENSIVE SAMPLING ===================
% The failure check deliberately raises AND catches errors. Use dbstop if
% error for debugging; dbstop if caught error also stops at these expected errors.
if run_implementation_checks
    test_thbi_rj();
    test_thbi_chain();
    test_thbi_failures();
    test_thbi_results();
end

%% =================== a2: KNOWN THIN EARTH AND ALL FIVE SYNTHETICS ===================
[data,TRUEmodel,forward]=a2_LOAD_DATA(par);
par.forc=forward.config;
save(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','forward','runtime','-v7');
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; v=data.(name);
    assert(numel(v.clean)==v.N && all(isfinite(v.clean)) && all(isfinite(v.obs)), ...
        'thbi:test','Invalid true/noisy synthetic data for %s.',name);
    fprintf('Synthetic %-5s: %d finite samples.\n',name,v.N);
end

%% =================== NATIVE SPLIT/MERGE CHECK ON THE ACTUAL TEST GRIDS ===================
% Exercise different physical layer counts before committing to 20,000 steps.
% Reuse one solver directory to check that CPS/HVf do not reuse stale models.
check_par=par; check_par.mod.start.jitter=false;
check_model=b1_INITIATE_MODEL(check_par);
[born,birth]=b2_PERTURB_MODEL(check_model,par,struct( ...
    'type','birth','zone',2,'index',1,'split_f',.12,'delta_vs',.08));
assert(birth.valid,'thbi:test','Native preflight birth is outside prior support.');
[back,death]=b2_PERTURB_MODEL(born,par,struct('type','death','zone',2,'index',1));
assert(death.valid && isequal(back.K,check_model.K) && ...
    max(abs(back.matrix(:)-check_model.matrix(:)))<1e-10, ...
    'thbi:test','Native preflight merge did not restore the physical Earth.');
preflight_dir=fullfile(result_dir,'native_split_merge'); mkdir(preflight_dir);
check_models={check_model,born,back}; check_predictions=cell(3,1);
for j=1:3
    try
        [check_predictions{j},check_forward]=b3_FORWARD_MODEL(check_models{j},data,par,preflight_dir);
        for d=1:numel(par.inv.datatypes)
            name=par.inv.datatypes{d}; y=check_predictions{j}.(name);
            assert(isreal(y) && numel(y)==data.(name).N && all(isfinite(y)), ...
                'thbi:test','Invalid native preflight response for %s.',name);
            assert(isequal(check_forward.pred.(name).model_input,check_forward.model), ...
                'thbi:test','Native solvers used different Earths for %s.',name);
        end
    catch problem
        context=struct('stage','native split/merge preflight','iteration',0,'model_index',j);
        preserve_forward_failure(fullfile(result_dir,'preflight_failure'), ...
            check_models{j},data,par,problem,preflight_dir,context);
        rethrow(problem);
    end
end
native_check=struct();
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d};
    native_check.(name).split_change=max(abs(check_predictions{1}.(name)-check_predictions{2}.(name)));
    native_check.(name).roundtrip_error=max(abs(check_predictions{1}.(name)-check_predictions{3}.(name)));
end
save(fullfile(preflight_dir,'preflight.mat'),'check_models','check_predictions','native_check','-v7');
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d};
    assert(native_check.(name).split_change>1e-8 && native_check.(name).roundtrip_error<1e-6, ...
        'thbi:test','Native split/merge response check failed for %s; inspect preflight.mat.',name);
end
fprintf('PASS: all five native responses change after a split and return after the inverse merge.\n');

%% =================== a3: CONDITIONED PRIOR FOR POSTERIOR COMPARISON ===================
if par.res.nprior>0
    prior=a3_BUILD_EMPIRICAL_PRIOR(par,par.res.nprior);
    save(fullfile(result_dir,'prior.mat'),'prior','-v7');
end

%% =================== b1-b9: TWO INDEPENDENT TRANSDIMENSIONAL CHAINS ===================
% Each chain saves chain_XX/checkpoint.mat and chain_XX/chain.mat.
% The sampler preserves the failed model and native files if a solver fails.
chains=run_all_chains(par,data);

%% =================== CHECK COMPLETE HISTORIES AND ACTUAL COUNT TRANSITIONS ===================
assert(numel(chains)==par.inv.nchains,'thbi:test','Wrong number of completed chains.');
keep=(par.inv.burnin+1):par.inv.niter;
test_report=struct('nchains',par.inv.nchains,'iterations_per_chain',par.inv.niter, ...
    'burnin',par.inv.burnin,'retained_per_chain',numel(keep), ...
    'retained_total',numel(keep)*numel(chains),'implementation_checks_run',run_implementation_checks);
test_report.native_split_merge=native_check;
test_report.forward_failures=zeros(1,par.inv.nchains);
test_report.retained_births=zeros(par.inv.nchains,3);
test_report.retained_deaths=zeros(par.inv.nchains,3);
for c=1:numel(chains)
    chain=chains{c};
    assert(chain.last_iter==par.inv.niter && isequal(size(chain.K),[par.inv.niter 3]), ...
        'thbi:test','Chain %d is incomplete.',c);
    values=[chain.logL chain.logP chain.E2 chain.logL_indivdat chain.boundaries chain.halfspace chain.eta];
    assert(isreal(values) && all(isfinite(values(:))) && all(chain.E2(:)>=0), ...
        'thbi:test','Chain %d contains invalid states or likelihoods.',c);
    assert(chain.failures==sum(chain.forward_failed) && ...
        ~any(chain.accepted & (~chain.proposal_valid | chain.forward_failed)), ...
        'thbi:test','Chain %d has inconsistent failure/acceptance accounting.',c);
    test_report.forward_failures(c)=chain.failures;
    % A stage-two layered path contains BOTH moves x->y->z. Count each
    % split/merge along an accepted path, although y itself was not accepted.
    expected_dK=zeros(par.inv.niter,3);
    layered=any(strcmp(par.proposal.dr.mode,{'layered_path','targeted_layered'}));
    for z=1:3
        first_used=chain.accepted_stage==1 | (layered & chain.accepted_stage==2);
        births=double(first_used & chain.move==4 & chain.zone==z);
        deaths=double(first_used & chain.move==5 & chain.zone==z);
        births=births+double(chain.accepted_stage==2 & chain.dr_move==4 & chain.dr_zone==z);
        deaths=deaths+double(chain.accepted_stage==2 & chain.dr_move==5 & chain.dr_zone==z);
        expected_dK(:,z)=births-deaths;
        test_report.retained_births(c,z)=sum(births(keep));
        test_report.retained_deaths(c,z)=sum(deaths(keep));
    end
    assert(isequal(diff([chain.initial.model.K;chain.K],1,1),expected_dK), ...
        'thbi:test','Chain %d count changes disagree with its accepted one/two-leg paths.',c);
    for i=1:par.inv.niter
        model=model_from_chain(chain,i,par);
        [supported,why]=a1_TEST_CONDITIONS(model,par);
        assert(supported,'thbi:test','Chain %d iteration %d is outside prior support: %s',c,i,why);
        [Nfinite,K,Nall]=count_N_layers_in_mod(model,par);
        assert(isequal(K,chain.K(i,:)) && Nall==Nfinite+1, ...
            'thbi:test','Chain %d iteration %d has inconsistent physical counts.',c,i);
    end
    expected_iterations=zeros(0,1);
    if par.res.prediction_stride>0
        expected_iterations=(par.inv.burnin+1:par.res.prediction_stride:par.inv.niter)';
    end
    assert(isequal(chain.prediction_iterations,expected_iterations) && all(chain.prediction_saved), ...
        'thbi:test','Chain %d is missing retained responses.',c);
    for d=1:numel(par.inv.datatypes)
        name=par.inv.datatypes{d};
        if par.res.prediction_stride==0
            assert(~isfield(chain.predictions,name),'thbi:test','Response saving should be disabled.');
            continue;
        end
        values=chain.predictions.(name);
        assert(isequal(size(values),[numel(expected_iterations) data.(name).N]) && ...
            isreal(values) && all(isfinite(values(:))), ...
            'thbi:test','Chain %d has invalid saved %s responses.',c,name);
    end
    fprintf('PASS: chain %d completed, valid states/count transitions and %d response draws.\n', ...
        c,numel(expected_iterations));
end
assert(~isequaln(chains{1}.velocities,chains{2}.velocities) || ...
    ~isequaln(chains{1}.eta,chains{2}.eta),'thbi:test','The two full histories are identical; inspect seeds.');

%% =================== c1-c4: POSTERIOR, MODEL SUITE AND RETAINED FINAL EARTH ===================
summary=c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');
[posterior,retained]=c2_BUILD_POSTERIOR(chains,par);
suite_of_models=c3_BUILD_MODEL_SUITE(retained,par);
final_model=c4_FINAL_MODEL(posterior,suite_of_models,retained,par);
assert(summary.bestfit_chain==final_model.representative_chain && ...
    summary.bestfit_iteration==final_model.representative_iteration, ...
    'thbi:test','Summary fit and retained final-model selection disagree.');
assert(summary.retained_per_chain==numel(keep) && posterior.Nstored==test_report.retained_total && ...
    size(suite_of_models.VS,2)==test_report.retained_total, ...
    'thbi:test','Postprocessing dropped or duplicated retained states.');
assert(final_model.representative_iteration>par.inv.burnin, ...
    'thbi:test','Final model was selected from burn-in.');
save(fullfile(result_dir,'postprocessing.mat'),'posterior','suite_of_models','final_model','-v7.3');

%% =================== c5: FRESH FORWARD CHECK AGAINST THE SAVED CHAIN ===================
% Recompute, rather than trusting cached predictions of the selected Earth.
final_forward=c5_FINAL_FORWARD_MODEL(final_model,par,data);
c=final_model.representative_chain; i=final_model.representative_iteration;
saved_E2=chains{c}.E2(i,:)'; fresh_E2=final_forward.misfit.E2(:);
test_report.final_E2_relative_error=max(abs(fresh_E2-saved_E2)./(1+abs(saved_E2)));
test_report.final_logL_relative_error=abs(final_forward.logL-chains{c}.logL(i))/(1+abs(chains{c}.logL(i)));
assert(test_report.final_E2_relative_error<1e-6 && test_report.final_logL_relative_error<1e-6, ...
    'thbi:test','Fresh final-model predictions disagree with the saved chain; inspect final_forward/.');

%% =================== PLOT ALL FIVE FITS, VELOCITIES, COUNTS AND MIXING ===================
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary,posterior,suite_of_models,final_model);

%% =================== SAVE EXECUTION DIAGNOSTICS WITHOUT CLAIMING CONVERGENCE ===================
test_report.completed=true;
test_report.all_native_calls_succeeded=all(test_report.forward_failures==0);
test_report.retained_birth_and_death_per_chain= ...
    sum(test_report.retained_births,2)>0 & sum(test_report.retained_deaths,2)>0;
test_report.status='Execution checks passed; inspect convergence and recovery separately.';
if ~test_report.all_native_calls_succeeded || ~all(test_report.retained_birth_and_death_per_chain)
    test_report.status='Run completed with diagnostic issues; inspect solver failures and RJ movement.';
end
test_report.metric_names=summary.metric_names;
test_report.classical_split_rhat=summary.split_rhat;
test_report.note='No convergence, unique recovery or uncertainty-coverage claim follows from finishing this test.';
save(fullfile(result_dir,'test_report.mat'),'test_report','-v7');
fid=fopen(fullfile(result_dir,'test_report.txt'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',test_report.status);
fprintf(fid,'2 chains x 10000 iterations; burn-in 4000; 6000 retained/chain; 12000 total.\n');
fprintf(fid,'Native split/merge, state/count consistency, postprocessing and final-forward checks passed.\n');
for c=1:numel(chains)
    fprintf(fid,'Chain %d: native failures %d; retained births [sed crust mantle] [%d %d %d]; deaths [%d %d %d].\n', ...
        c,test_report.forward_failures(c),test_report.retained_births(c,:),test_report.retained_deaths(c,:));
end
fprintf(fid,'Final forward relative errors: E2 %.6g; logL %.6g.\n', ...
    test_report.final_E2_relative_error,test_report.final_logL_relative_error);
fprintf(fid,'See inversion_report.txt for fits, count masses, per-move acceptance and split R-hat.\n%s\n',test_report.note);
fclose(fid);
fid=fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',result_dir); fclose(fid);
if ~test_report.all_native_calls_succeeded
    warning('thbi:testDiagnostics','Native proposals failed: chain totals %s. Inspect first_forward_failure folders.', ...
        mat2str(test_report.forward_failures));
end
if ~all(test_report.retained_birth_and_death_per_chain)
    warning('thbi:testDiagnostics','At least one chain did not accept both birth and death moves after burn-in; inspect mixing.');
end
fprintf('\n%s\nResults and figures: %s\n',test_report.status,result_dir);
% After interruption, keep code/settings unchanged and call:
% addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES(); resume_thbi(result_dir);
% resume_thbi continues sampling and ordinary postprocessing. It does not run
% this launcher's extra validation/report blocks or the fresh c5 forward check.
