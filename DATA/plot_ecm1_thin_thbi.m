% Replot a completed transdimensional run without running MCMC or any solver.
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();
result_dir=''; % enter an absolute run folder, or leave empty for latest
if isempty(result_dir)
    pointer=fullfile(repo,'DATA','transdimensional','LATEST_RUN.txt');
    assert(exist(pointer,'file')==2,'No completed transdimensional run yet.');
    result_dir=strtrim(fileread(pointer));
end
s=load(fullfile(result_dir,'results.mat'));
s.par.res.dir=result_dir; s.par.res.show_figures=true;
if exist(fullfile(result_dir,'postprocessing.mat'),'file')
    p=load(fullfile(result_dir,'postprocessing.mat'));
    c2_PLOT_RESULTS(s.par,s.data,s.TRUEmodel,s.chains,s.summary,p.posterior,p.suite_of_models,p.final_model);
else
    c2_PLOT_RESULTS(s.par,s.data,s.TRUEmodel,s.chains,s.summary);
end
