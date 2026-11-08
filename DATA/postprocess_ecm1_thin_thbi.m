% Rebuild the full THBI-style posterior/model/plot suite from saved chains.
% No MCMC. Native solvers run only if run_final_forward is set true below.
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();
result_dir=''; % enter an absolute completed run directory, or use latest
run_final_forward=false;
if isempty(result_dir)
    pointer=fullfile(repo,'DATA','transdimensional','LATEST_RUN.txt');
    assert(exist(pointer,'file')==2,'No completed transdimensional run yet.');
    result_dir=strtrim(fileread(pointer));
end
s=load(fullfile(result_dir,'results.mat'));
par=s.par; par.res.dir=result_dir; par.res.show_figures=true;
chains=s.chains; data=s.data; TRUEmodel=s.TRUEmodel;

%% =================== c1-c4: ALL RETAINED STATES AND COMMON-DEPTH SUMMARIES ===================
summary=c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'summary','-append');
[posterior,retained]=c2_BUILD_POSTERIOR(chains,par);
suite_of_models=c3_BUILD_MODEL_SUITE(retained,par);
final_model=c4_FINAL_MODEL(posterior,suite_of_models,retained,par);
save(fullfile(result_dir,'postprocessing.mat'),'posterior','suite_of_models','final_model','-v7.3');

%% =================== MODEL, DATA, COUNTS, TRADEOFFS AND SAMPLING FIGURES ===================
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary,posterior,suite_of_models,final_model);

%% =================== c5: OPTIONAL FIVE-OBSERVABLE RECOMPUTATION ===================
if run_final_forward
    final_forward=c5_FINAL_FORWARD_MODEL(final_model,par,data);
end
fprintf('Saved posterior/model suite and figures in:\n%s\n',result_dir);
