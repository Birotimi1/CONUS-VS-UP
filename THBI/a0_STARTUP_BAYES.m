function [repo, paths] = a0_STARTUP_BAYES()
% Set the project, MATLAB, native-solver, data and output paths.
% Analogous to THBI a0_STARTUP_BAYES, restricted to dependencies used here.
% The optional paths output exposes each directory for inspection in MATLAB.

%% =================== PROJECT DIRECTORIES ===================
paths.THBI = fileparts(mfilename('fullpath'));
repo = fileparts(paths.THBI);
paths.project = repo;
paths.parms = fullfile(paths.THBI,'parms');
paths.functions = fullfile(repo,'functions');
paths.plotting = fullfile(paths.functions,'plotting');
paths.tests = fullfile(paths.THBI,'tests');

%% =================== MATLAB SEARCH PATH ===================
% Keep these explicit so a breakpoint shows exactly which code is loaded.
addpath(paths.plotting);
addpath(paths.functions);
addpath(paths.parms);
addpath(paths.tests);
addpath(paths.THBI); % numbered stages take precedence over other checkouts

%% =================== SHARED FORWARD-SOLVER PROJECT ===================
paths.shared_project = getenv('HBI_SHARED_PROJECT');
if isempty(paths.shared_project)
    paths.shared_project = fullfile(fileparts(fileparts(repo)),'HVSR-Joint-Inversion');
end

%% =================== CPS: DISPERSION AND RECEIVER FUNCTIONS ===================
paths.cps_bin = getenv('HBI_CPS_BIN');
if isempty(paths.cps_bin)
    candidates = [{fullfile(repo,'AGU_models','bin_v3.30')}, ...
        strsplit(getenv('PATH'),':'), ...
        {fullfile(paths.shared_project,'AGU_models','bin_v3.30')}];
    for j = 1:numel(candidates)
        if ~isempty(candidates{j}) && exist(fullfile(candidates{j},'surf96'),'file')
            paths.cps_bin = candidates{j};
            break;
        end
    end
end
paths.surf96 = fullfile(paths.cps_bin,'surf96');
paths.hrftn96 = fullfile(paths.cps_bin,'hrftn96');

%% =================== HV-DFA: HVSR AND ELLIPTICITY ===================
paths.hvf_bin = getenv('HBI_HVF_BIN');
if isempty(paths.hvf_bin)
    candidates = [{fullfile(repo,'AGU_models','bin_HV_DFA')}, ...
        strsplit(getenv('PATH'),':'), ...
        {fullfile(repo,'AGU_models','bin_v3.30'), ...
         fullfile(paths.shared_project,'AGU_models','bin_HV_DFA')}];
    for j = 1:numel(candidates)
        if ~isempty(candidates{j}) && exist(fullfile(candidates{j},'HVf'),'file')
            paths.hvf_bin = candidates{j};
            break;
        end
    end
end
paths.HVf = fullfile(paths.hvf_bin,'HVf');
% b3 checks the executables when forward calculation is requested, and also
% supports explicit par.forc overrides and executables on the system PATH.

%% =================== DATA, RESULTS AND LOCAL SCRATCH ===================
paths.data = fullfile(repo,'DATA');
paths.inversions = fullfile(paths.data,'transdimensional');
paths.synthetics = fullfile(paths.data,'synthetics');
paths.scratch = tempdir;
% Native solver files use a separate temporary directory for every chain.
% Detect accidental use of the fixed-dimensional or external THBI stages.
required = {'a1_TEST_CONDITIONS','a2_LOAD_DATA','a2_PREPARE_DATA', ...
    'a3_CALC_LOG_PRIOR','a3_BUILD_EMPIRICAL_PRIOR','b1_INITIATE_MODEL', ...
    'b2_PERTURB_MODEL','b3_FORWARD_MODEL','b4_CALC_MISFIT','b5_CALC_LIKELIHOOD', ...
    'b6_IFACCEPT','b9_SAVE_RESULT','run_one_chain','run_all_chains', ...
    'c1_PROCESS_RESULTS','c2_PLOT_RESULTS','c2_BUILD_POSTERIOR', ...
    'c3_BUILD_MODEL_SUITE','c4_FINAL_MODEL','c5_FINAL_FORWARD_MODEL'};
for j = 1:numel(required)
    assert(strcmp(fileparts(which(required{j})),paths.THBI), ...
        'thbi:path','Wrong implementation on MATLAB path: %s',required{j});
end
assert(strcmp(fileparts(which('bayes_inv_parms')),paths.parms), ...
    'thbi:path','Wrong bayes_inv_parms on MATLAB path.');
% Check every packaged helper, including file IO and plotting utilities.
% A current folder or another checkout can otherwise shadow a helper even
% when the numbered stages themselves resolve correctly.
helpers=[dir(fullfile(paths.functions,'*.m'));dir(fullfile(paths.plotting,'*.m'))];
assert(~isempty(helpers),'thbi:path','Missing project functions folder: %s',paths.functions);
for j=1:numel(helpers)
    [~,name]=fileparts(helpers(j).name);
    expected=fullfile(helpers(j).folder,helpers(j).name);
    assert(strcmp(which(name),expected),'thbi:path', ...
        'Wrong or missing helper on MATLAB path: %s. Expected %s.',name,expected);
end
% This startup stage does not launch solvers, clear variables, or delete files.
end
