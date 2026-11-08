%% =================== TARGETED TRANSDIMENSIONAL ACCEPTANCE CHECKS ===================
% Optional targeted comparison; the production profiles use random second moves.
% Run in MATLAB from this script. No CPS/HVf or production inversion.
% Also executable in Octave for the sampler kernel; MATLAB-only posterior
% summaries and exception-report formatting are tested only in MATLAB.
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();
is_octave=exist('OCTAVE_VERSION','builtin')~=0;

%% =================== PHYSICAL MAPS AND FULL ACCEPTANCE RULE ===================
test_thbi_rj();
test_thbi_layered_dr();                 % 1080 physical forward/reverse paths
balance_report=test_thbi_dr_stationarity('targeted_layered'); % exact transition matrix + negative controls

%% =================== THREE KERNELS: CACHE, COUNTS AND CHECKPOINTS ===================
test_thbi_chain('extended60','targeted_layered',~is_octave);
test_thbi_chain('extended60','layered_path',~is_octave);
test_thbi_chain('extended60','fixed_scalar',~is_octave);
if ~is_octave
    test_thbi_failures('extended60');
    test_thbi_results();
end

%% =================== INDEPENDENT CONTINUOUS PRIOR ENSEMBLE ===================
% Unlike one long chain, exact independent starts do not require mixing.
prior_ensemble_report=test_thbi_prior_invariance(6000,'targeted_layered');

%% =================== 250-KM PROFILE AND LONG-PERIOD RF PROCESSING ===================
test_thbi_extended60();
coverage_report=test_thbi_extended100();
test_thbi_chain('extended100','targeted_layered',~is_octave);

%% =================== OPTIONAL CONTINUOUS PRIOR-COUNT CHECK ===================
run_prior_counts=false; % 30,000 solver-free iterations; 5,000 burn-in
if run_prior_counts, test_thbi_prior_counts('extended60'); end
fprintf('Targeted RJ checks completed; native forward validation is separate.\n');
