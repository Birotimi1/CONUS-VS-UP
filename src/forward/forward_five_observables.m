function forward = forward_five_observables(model, grids, cfg, output_dir)
% Forward-only evaluation: one layered Earth, five observables, no data fit.
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
    'project_root', fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
    'shared_project', '', 'max_layer_thickness_km', 2));
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
    fprintf('Calculating %s from the shared %d-layer model...\n', dt, size(model, 1));
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
fprintf('Calculating synthetic P-to-S receiver function...\n');
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

function a = merge_defaults(a, b)
names = fieldnames(b);
for k = 1:numel(names)
    if ~isfield(a, names{k}), a.(names{k}) = b.(names{k}); end
end
end

function directory = find_bin(explicit, env_name, candidates, executable)
directory = explicit;
if isempty(directory), directory = getenv(env_name); end
if ~isempty(directory)
    if exist(fullfile(directory, executable), 'file'), return; end
    error('fiveobs:dependency', '%s is not in configured directory %s.', executable, directory);
end
% Prefer project-local binaries, then PATH, then the supplied shared project.
search = [candidates(1) strsplit(getenv('PATH'), ':') candidates(2:end)];
for k = 1:numel(search)
    if ~isempty(search{k}) && exist(fullfile(search{k}, executable), 'file')
        directory = search{k}; return;
    end
end
error('fiveobs:dependency', 'Cannot find %s. Set %s or its cfg bin directory.', executable, env_name);
end

function prediction = hv_prediction(model, x, dtype, cfg, workdir)
previous = pwd; cleanup = onCleanup(@() cd(previous)); %#ok<NASGU>
cd(workdir);
fid = fopen('model.txt', 'w');
if fid < 0, error('fiveobs:io', 'Cannot write HV-DFA model.'); end
fprintf(fid, '%d\n', size(model, 1));
fprintf(fid, '%.12g %.12g %.12g %.12g\n', (1000 * model)');
fclose(fid);
if strcmp(dtype, 'hvsr')
    frequencies = x;
    nmr = cfg.hvsr.nmr; nml = cfg.hvsr.nml; nks = cfg.hvsr.nks;
    apsv = cfg.hvsr.apsv; ash = cfg.hvsr.ash;
else
    frequencies = 1 ./ x;
    nmr = 1; nml = 0; nks = 0; apsv = 0; ash = 0;
end
[requested, ~, map] = unique(frequencies);
fid = fopen('frequencies.txt', 'w');
fprintf(fid, '%.15g\n', requested); fclose(fid);
cmd = sprintf('%s -f model.txt -ff frequencies.txt -nmr %d -nml %d -nks %d -prec %.10g -apsv %.10g -ash %.10g -rep -hv', ...
    quote_shell_path(fullfile(cfg.hvf_bin, 'HVf')), nmr, nml, nks, cfg.hvsr.prec, apsv, ash);
[status, message] = system([cmd ' > HV.dat 2> HVf.log']);
if status ~= 0, error('fiveobs:hvf', 'HVf failed (%d): %s. See %s', status, message, workdir); end
fid = fopen('HV.dat', 'r');
if fid < 0, error('fiveobs:hvf', 'HVf did not produce output.'); end
values = fscanf(fid, '%f'); fclose(fid);
if numel(values) ~= 2 * numel(requested)
    error('fiveobs:hvf', 'HVf returned an incomplete frequency grid; see %s', fullfile(workdir, 'HV.dat'));
end
values = reshape(values, 2, [])';
if any(abs(values(:, 1) - requested) > 1e-6 * max(1, abs(requested)))
    error('fiveobs:hvf', 'HVf returned a different frequency grid.');
end
prediction = struct('val', values(map, 2), 'model_input', model, 'command', cmd, 'status', 'computed');
end

function prediction = dispersion_prediction(model, periods, dtype, cfg, workdir)
previous = pwd; cleanup = onCleanup(@() cd(previous)); %#ok<NASGU>
cd(workdir);
write_elastic_cps_model('start.mod', model);
if strcmp(dtype, 'cph'), velocity_type = 'C'; else, velocity_type = 'U'; end
% CPS exchange files print periods to four decimals. Match by returned period,
% never by position or by truncating the requested observations.
[requested, ~, map] = unique(round(periods * 1e4) / 1e4);
fid = fopen('disp_obs.dsp', 'w');
for k = 1:numel(requested)
    fprintf(fid, 'SURF96 R %s X 0 %.4f 3.0 0.1\n', velocity_type, requested(k));
end
fclose(fid);
fid = fopen('sobs.d', 'w');
fprintf(fid, '0.005 0.005 0.0 0.005 0.0\n0 0 0 0 0 1 1 0 1 0\nstart.mod\ndisp_obs.dsp\n');
fclose(fid);
if exist('predicted.dsp', 'file'), delete('predicted.dsp'); end
cmd = [quote_shell_path(fullfile(cfg.cps_bin, 'surf96')) ' 1 27 predicted.dsp'];
[status, message] = system([cmd ' > surf96.log 2>&1']);
if status ~= 0 || ~exist('predicted.dsp', 'file')
    error('fiveobs:cps', 'surf96 failed (%d): %s. See %s', status, message, fullfile(workdir, 'surf96.log'));
end
fid = fopen('predicted.dsp', 'r');
columns = textscan(fid, '%s %s %s %s %f %f %f %f'); fclose(fid);
if isempty(columns{6}) || any(~strcmp(columns{2}, 'R')) || ...
        any(~strcmp(columns{3}, velocity_type)) || any(columns{5} ~= 0)
    error('fiveobs:cps', 'Unexpected CPS wave/mode/velocity type or empty output.');
end
predicted = NaN(size(requested));
for k = 1:numel(requested)
    rows = find(abs(columns{6} - requested(k)) < 6e-5);
    if numel(rows) == 1, predicted(k) = columns{7}(rows); end
end
prediction = struct('val', predicted(map), 'model_input', model, 'command', cmd, 'status', 'computed');
end
