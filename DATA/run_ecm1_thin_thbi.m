% RUN_ECM1_THIN_THBI -- run manually in MATLAB when ready for synthetic testing.
% Three variable-layer zones, one physical Earth, five observable datasets.

%% =================== a0: PROJECT AND SOLVER PATHS ===================
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI'));
[repo,paths]=a0_STARTUP_BAYES();

%% =================== RUN SETTINGS AND ZONED PRIORS ===================
par=bayes_inv_parms('pilot');
par.inv.parallel=false;
par.res.show_figures=true;
run_final_forward=false; % true recomputes the selected retained Earth's five responses after MCMC
% Edit THBI/parms/bayes_inv_parms.m for counts, boundaries, moves and errors.
% par.inv.niter=30000; par.inv.burnin=10000; par.inv.nchains=4;
validate_configuration(par);
assert(~par.inv.prior_only,'Use DATA/check_thbi_implementation.m for prior-only checks.');

%% =================== CREATE A NEW OUTPUT DIRECTORY ===================
if ~exist(par.res.root,'dir'), mkdir(par.res.root); end
result_dir=fullfile(par.res.root,['ECM1_thin_THBI_' par.inv.profile '_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(result_dir,'dir'),'Result directory already exists.');
mkdir(result_dir); par.res.dir=result_dir;
fprintf('Zoned RJ-HBI output: %s\n',result_dir);

%% =================== a2: KNOWN THIN EARTH AND FIVE SYNTHETICS ===================
[data,TRUEmodel,forward]=a2_LOAD_DATA(par);
par.forc=forward.config;
runtime.version=version; runtime.computer=computer;
save(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','forward','runtime','-v7');

%% =================== a3: EMPIRICAL CONDITIONED PRIOR ===================
if par.res.nprior>0
    prior=a3_BUILD_EMPIRICAL_PRIOR(par,par.res.nprior);
    save(fullfile(result_dir,'prior.mat'),'prior','-v7');
end

%% =================== b1-b9: INITIALIZATION AND RJ CHAINS ===================
chains=run_all_chains(par,data);

%% =================== c1: COMMON-DEPTH AND COUNT SUMMARIES ===================
summary=c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');

%% =================== c2-c4: POSTERIOR, MODEL SUITE AND PHYSICAL FINAL MODEL ===================
[posterior,retained]=c2_BUILD_POSTERIOR(chains,par);
suite_of_models=c3_BUILD_MODEL_SUITE(retained,par);
final_model=c4_FINAL_MODEL(posterior,suite_of_models,retained,par);
save(fullfile(result_dir,'postprocessing.mat'),'posterior','suite_of_models','final_model','-v7.3');

%% =================== c2: FITS, COUNTS, SHALLOW STRUCTURE AND DIAGNOSTICS ===================
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary,posterior,suite_of_models,final_model);

%% =================== c5: OPTIONAL FRESH FORWARD CHECK OF THE RETAINED MODEL ===================
if run_final_forward
    final_forward=c5_FINAL_FORWARD_MODEL(final_model,par,data);
end
fid=fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',result_dir); fclose(fid);
fprintf('Completed. Inspect counts and mixing before interpreting recovery:\n%s\n',result_dir);
