% One prescribed ECM1 thin Earth -> five clean forward observables.
% Run from MATLAB or Octave. Results stay under this repository's DATA folder.
repo = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo, 'functions'));
cfg = struct();
cfg.project_root = repo;
cfg.output_root = getenv('HBI_SYNTH_FIVE_OUTPUT');
if isempty(cfg.output_root), cfg.output_root = fullfile(repo, 'DATA', 'synthetics'); end
cfg.shared_project = getenv('HBI_SHARED_PROJECT');
if isempty(cfg.shared_project)
    % SharedData/US_TA/CONUS-VS-UP and SharedData/HVSR-Joint-Inversion.
    cfg.shared_project = fullfile(fileparts(fileparts(repo)), 'HVSR-Joint-Inversion');
end
% HBI_CPS_BIN and HBI_HVF_BIN can override executable discovery on any machine.
cfg.max_layer_thickness_km = 2;
cfg.rf = struct('rayp_s_km', 0.06, 'gaussian_alpha', 2.5, ...
    'dt', 0.025, 'npts', 4096, 'delay', 10, 'plot_window', [-2 25]);
manifest = run_synthetic_thin(cfg);
fprintf('\nPlot with: python3 "%s" "%s"\n', ...
    fullfile(repo, 'scripts', 'plot_synthetic_thin.py'), manifest.output_dir);
