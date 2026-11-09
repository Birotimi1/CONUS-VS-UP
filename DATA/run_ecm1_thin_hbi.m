% RUN_ECM1_THIN_HBI
% Run this script in MATLAB. One thin layered Earth, five observables.
% Generates noisy synthetics, runs fixed-layer HBI, saves and opens figures.

%% Paths
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'HBI'));
[repo,paths] = a0_STARTUP_BAYES();

%% Run settings
% 'pilot': 2 chains, 240 iterations, 120 burn-in (implementation check only).
% 'production': 4 chains, 30,000 iterations, 10,000 burn-in.
par = a1_DEFINE_PARAMETERS('pilot');
par.inv.parallel = false;       % true uses MATLAB process workers (PCT)
par.res.show_figures = true;

% Optional overrides; uncomment and edit when ready for a longer inversion.
% par.inv.niter = 30000;
% par.inv.burnin = 10000;
% par.inv.nchains = 4;

% Model bounds, initial model, data errors, RF settings, and proposal sizes
% are in HBI/a1_DEFINE_PARAMETERS.m. The HBI thin truth is defined directly
% in HBI/a2_LOAD_DATA.m. It is not used to initialize the chains.

%% =================== CHECK RUN LENGTHS ===================
assert(par.inv.niter > par.inv.burnin && par.inv.burnin >= 0 && ...
    par.inv.nchains >= 1 && ...
    all(mod([par.inv.niter par.inv.burnin par.inv.nchains],1)==0), ...
    'Invalid chain lengths.');

%% =================== CREATE THIS RUN'S OUTPUT DIRECTORY ===================
if ~exist(par.res.root,'dir'), mkdir(par.res.root); end
result_dir = fullfile(par.res.root,['ECM1_thin_HBI_' par.inv.profile '_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(result_dir,'dir'),'Result directory already exists.');
mkdir(result_dir);
par.res.dir = result_dir;
fprintf('Output: %s\n',result_dir);

%% =================== a2: TRUE MODEL AND FIVE SYNTHETIC DATASETS ===================
[data,TRUEmodel,forward] = a2_LOAD_DATA(par);
par.forc = forward.config;
runtime.version = version;
runtime.computer = computer;
save(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','forward','runtime','-v7');

%% =================== b1-b9: INITIALIZE AND RUN THE HBI CHAINS ===================
% run_one_chain calls the numbered condition/prior, proposal, forward,
% misfit, likelihood, acceptance and state-saving stages on each iteration.
chains = run_all_chains(par,data);

%% =================== c1: PROCESS ALL RETAINED CHAIN STATES ===================
summary = c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');

%% =================== c2: PLOT MODEL, OBSERVABLES AND DIAGNOSTICS ===================
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary);

%% =================== RECORD THE COMPLETED RUN ===================
fid = fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w');
assert(fid >= 0,'Cannot write latest-run pointer.');
fprintf(fid,'%s\n',result_dir);
fclose(fid);
fprintf('Saved results and figures in:\n%s\n',result_dir);
