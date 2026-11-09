function forward = forward_five_observables(model, grids, cfg, output_dir)
% Evaluate one layered Earth for all five observables (no likelihood here).
% grids.hvsr.f (Hz); grids.ellip.T/cph.T/ugr.T (s). RF grid is cfg.rf.
% Model columns: thickness km, Vp km/s, Vs km/s, density g/cm^3.
validate_layer_model(model);
if nargin < 3, cfg = struct(); end
if nargin < 4, output_dir = tempname; end
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
% Absolute solver paths are necessary because each executable changes cwd.
[ok, attrs] = fileattrib(output_dir);
if ~ok, error('fiveobs:io', 'Cannot access output directory.'); end
output_dir = attrs.Name;
cfg = merge_defaults(cfg, struct('cps_bin', '', 'hvf_bin', '', ...
    'project_root', fileparts(fileparts(mfilename('fullpath'))), ...
    'shared_project', '', 'max_layer_thickness_km', 2, 'verbose', true));
if ~isfield(cfg, 'hvsr'), cfg.hvsr = struct(); end
cfg.hvsr = merge_defaults(cfg.hvsr, struct('nmr', 10, 'nml', 10, 'nks', 1000, ...
    'prec', 1e-4, 'apsv', 1e-3, 'ash', 1e-3));
if ~isfield(cfg, 'rf'), cfg.rf = struct(); end
cfg.rf = merge_defaults(cfg.rf, struct('rayp_s_km', 0.06, 'gaussian_alpha', 2.5, ...
    'dt', 0.025, 'npts', 4096, 'delay', 10, 'plot_window', [-2 25]));
cfg.cps_bin = find_bin(cfg.cps_bin, 'HBI_CPS_BIN', ...
    {fullfile(cfg.project_root, 'AGU_models', 'bin_v3.30'), ...
     fullfile(cfg.shared_project, 'AGU_models', 'bin_v3.30')}, 'surf96');
cfg.hvf_bin = find_bin(cfg.hvf_bin, 'HBI_HVF_BIN', ...
    {fullfile(cfg.project_root, 'AGU_models', 'bin_HV_DFA'), ...
     fullfile(cfg.project_root, 'AGU_models', 'bin_v3.30'), ...
     fullfile(cfg.shared_project, 'AGU_models', 'bin_HV_DFA')}, 'HVf');
if ~exist(fullfile(cfg.cps_bin, 'hrftn96'), 'file')
    error('fiveobs:dependency', 'hrftn96 is missing from %s.', cfg.cps_bin);
end
old_path = getenv('PATH');
restore_path = onCleanup(@() setenv('PATH', old_path)); %#ok<NASGU>
setenv('PATH', [cfg.cps_bin ':' cfg.hvf_bin ':' old_path]);
source_model = model;
model = refine_layer_model(model, cfg.max_layer_thickness_km);
forward = struct('source_model', source_model, 'model', model, ...
    'config', cfg, 'output_dir', output_dir, 'pred', struct());
forward.geometry = 'Flat, isotropic, elastic; identical layers for all five observables';
names = {'hvsr', 'ellip', 'cph', 'ugr'};
coords = {'f', 'T', 'T', 'T'};
% Validate all requested grids before invoking any executable.
for k = 1:4
    dt = names{k}; coord = coords{k};
    if ~isfield(grids, dt) || ~isfield(grids.(dt), coord)
        error('fiveobs:grid', 'Missing grids.%s.%s.', dt, coord);
    end
    x = grids.(dt).(coord);
    if ~isnumeric(x) || ~isreal(x) || ~isvector(x) || isempty(x) || any(~isfinite(x) | x <= 0)
        error('fiveobs:grid', '%s coordinates must be finite positive vectors.', dt);
    end
    if any(strcmp(dt, {'cph', 'ugr'})) && any(round(x * 1e4) <= 0)
        error('fiveobs:grid', 'Periods must remain positive at CPS exchange precision.');
    end
end
started = tic;
for k = 1:4
    dt = names{k}; x = grids.(dt).(coords{k})(:);
    if cfg.verbose
        fprintf('Calculating %s from the shared %d-layer model...\n', dt, size(model, 1));
    end
    workdir = fullfile(output_dir, dt);
    if ~exist(workdir, 'dir'), mkdir(workdir); end
    try
        if any(strcmp(dt, {'hvsr', 'ellip'}))
            prediction = hv_prediction(model, x, dt, cfg, workdir);
        else
            prediction = dispersion_prediction(model, x, dt, cfg, workdir);
        end
        if numel(prediction.val) ~= numel(x) || any(~isfinite(prediction.val) | prediction.val <= 0)
            error('fiveobs:incomplete', '%s has missing or nonphysical predictions; see %s.', dt, workdir);
        end
        prediction.status = 'ok';
    catch problem
        prediction = struct('val', NaN(size(x)), 'status', 'failed', ...
            'message', problem.message, 'model_input', model);
        warning('fiveobs:solver', '%s: %s', dt, problem.message);
    end
    prediction.x = x;
    forward.pred.(dt) = prediction;
end
if cfg.verbose, fprintf('Calculating synthetic P-to-S receiver function...\n'); end
try
    forward.pred.rf = forward_receiver_function_cps(model, cfg.rf, cfg.cps_bin, fullfile(output_dir, 'rf'));
catch problem
    forward.pred.rf = struct('t', [], 'val', [], 'status', 'failed', ...
        'message', problem.message, 'model_input', model, 'observed', false);
    warning('fiveobs:solver', 'Receiver function: %s', problem.message);
end
forward.all_forward_succeeded = strcmp(forward.pred.rf.status, 'ok');
for k = 1:4
    forward.all_forward_succeeded = forward.all_forward_succeeded && strcmp(forward.pred.(names{k}).status, 'ok');
end
forward.elapsed_seconds = toc(started);
end
