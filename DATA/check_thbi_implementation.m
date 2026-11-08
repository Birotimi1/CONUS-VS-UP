% Run these checks manually in MATLAB before a long synthetic RJ inversion.
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();

%% =================== ALGEBRAIC AND CHEAP CHAIN CHECKS ===================
test_thbi_rj();
test_thbi_chain();
test_thbi_failures();
test_thbi_results();

%% =================== OPTIONAL PRIOR-ONLY COUNT EXPERIMENT ===================
run_prior_check=false; % 30,000 cheap iterations; no native forward calls
if run_prior_check, test_thbi_prior_counts(); end

%% =================== OPTIONAL CPS/HVf VARIABLE-LAYER CHECK ===================
run_native_check=false; % enables real forward calculations, no inversion
if run_native_check, test_thbi_forward(); end
