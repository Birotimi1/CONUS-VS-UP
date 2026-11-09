% GENERATE_ECM1_AK135_FIVE_OBSERVABLES -- forward-only MATLAB script.
% No MCMC. Produces the published thin ECM1 + labelled AK135 synthetic inputs.

%% =================== SETTINGS AND OUTPUT ===================
repo=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(repo,'THBI'));
a0_STARTUP_BAYES(); par=bayes_inv_parms('reference100_vti');
par.res.show_figures=true;
par.res.dir=fullfile(par.res.root,['ECM1_AK135_forward_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(par.res.dir,'dir'),'thbi:run','Output directory already exists.');
mkdir(par.res.dir); validate_configuration(par);
test_thbi_reference100();
test_thbi_vti();
vti_report=test_thbi_vti_native(par,fullfile(par.res.dir,'vti_preflight'));
assert(vti_report.passed);

%% =================== a2: TRUE STRUCTURE, CPS/HVf AND FIXED RF PRODUCTS ===================
[data,TRUEmodel,forward]=a2_LOAD_DATA(par); par.forc=forward.config;
save(fullfile(par.res.dir,'inputs.mat'),'par','data','TRUEmodel','forward','vti_report','-v7');

%% =================== CURVES, Ps WAVES/STACK AND Sp DEPTH STACK ===================
plot_reference_synthetics(TRUEmodel,data,par); plot_RF_BANDS(forward,par,'true');
fprintf('%s\nSaved synthetic inputs and figures: %s\n',TRUEmodel.description,par.res.dir);
disp(array2table(TRUEmodel.model_matrix,'VariableNames', ...
    {'Thickness_km','Vp_average_km_s','Vs_average_km_s','Density_g_cm3','xi','phi','elastic_eta'}));
