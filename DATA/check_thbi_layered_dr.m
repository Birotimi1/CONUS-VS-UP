%% =================== LAYERED TWO-STAGE RJ CHECKS ONLY ===================
% Run this script yourself in MATLAB. No CPS/HVf or production inversion.
% Main calculation blocks remain in THBI/b2 and THBI/run_one_chain; this
% script calls the algebra, analytic-chain, resume, failure and result checks.
is_octave=exist('OCTAVE_VERSION','builtin')~=0;
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();

%% =================== SPLIT/MERGE MAPS AND BOTH ACCEPTANCE KERNELS ===================
test_thbi_rj();
test_thbi_layered_dr();
test_thbi_dr_stationarity();
test_thbi_extended60(); % also pure-MATLAB RF/filter checks; no native solvers

%% =================== CACHE, CHAIN HISTORY AND EXACT RESTART ===================
test_thbi_chain('extended60','layered_path',~is_octave);
test_thbi_chain('extended100','layered_path',~is_octave);
coverage_report=test_thbi_extended100();
prior_ensemble_report=test_thbi_prior_invariance();
if ~is_octave
    test_thbi_failures('extended60');
    test_thbi_failures('extended100');
    test_thbi_results();
end

%% =================== OPTIONAL LONGER PRIOR-ONLY COUNT EXPERIMENT ===================
run_prior_counts=false; % 30k cheap prior-only iterations, not native inversion
if run_prior_counts
    test_thbi_prior_counts('extended60');
end
fprintf('Layered DR implementation checks finished. Review outputs before native preflight.\n');
