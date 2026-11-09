function prediction = dispersion_prediction(model, periods, dtype, cfg, workdir)
assert(size(model,2)==4,'fiveobs:anisotropy', ...
    'Legacy surf96 wrapper is isotropic. Use THBI/b3_FORWARD_MODEL for VTI.');
previous = pwd; cleanup = onCleanup(@() cd(previous)); %#ok<NASGU>
cd(workdir);
% surf96 reuses tmpsrfi.* (including the old model) even after start.mod is
% replaced. Equivalent to `surf96 39`: clear only CPS's generated cache.
% Essential when MCMC evaluates successive models in the same directory.
cache = dir('tmpsrfi.*');
for k = 1:numel(cache)
    if ~cache(k).isdir, delete(cache(k).name); end
end
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
