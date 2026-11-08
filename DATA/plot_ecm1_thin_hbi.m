% Reopen plots from a completed inversion, without running any chain/solver.
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'HBI'));
a0_STARTUP_BAYES();

% Leave empty for the latest completed run, or enter an absolute run path.
result_dir = '';
if isempty(result_dir)
    pointer = fullfile(repo,'DATA','inversions','LATEST_RUN.txt');
    assert(exist(pointer,'file')==2,'No completed run yet. Run run_ecm1_thin_hbi.m first.');
    result_dir = strtrim(fileread(pointer));
end
saved = load(fullfile(result_dir,'results.mat'));
saved.par.res.dir = result_dir;
saved.par.res.show_figures = true;
c2_PLOT_RESULTS(saved.par,saved.data,saved.TRUEmodel,saved.chains,saved.summary);
