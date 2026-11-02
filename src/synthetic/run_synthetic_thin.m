function manifest = run_synthetic_thin(cfg)
% Compute all five observables from one known layered Earth; no inversion.
true_model = ecm1_thin_true_model();
grids.hvsr.f = logspace(log10(.2), log10(25), 80)';
grids.ellip.T = linspace(3, 25, 50)';
grids.cph.T = linspace(9, 35, 35)';
grids.ugr.T = linspace(3, 25, 40)';
if ~exist(cfg.output_root, 'dir'), mkdir(cfg.output_root); end
run_dir = fullfile(cfg.output_root, ['ECM1_thin_' datestr(now, 'yyyymmdd_HHMMSS')]);
if exist(run_dir, 'dir'), run_dir = tempname(cfg.output_root); end
mkdir(run_dir);
[ok, attrs] = fileattrib(run_dir);
if ~ok, error('fiveobs:io', 'Cannot access output directory.'); end
run_dir = attrs.Name;
fprintf('%s\nOutput: %s\n', true_model.description, run_dir);
forward = forward_five_observables(true_model.model_matrix, grids, cfg, run_dir);
if ~forward.all_forward_succeeded
    save(fullfile(run_dir, 'failed_forward.mat'), 'forward', 'true_model', '-v7');
    error('fiveobs:incomplete', 'One or more observables failed; diagnostics retained in %s.', run_dir);
end
true_model.solver_model = forward.model;
true_model.geometry = forward.geometry;
syn_data = grids;
syn_data.rf.t = forward.pred.rf.t;
syn_data.rf.rayp_s_km = forward.config.rf.rayp_s_km;
syn_data.rf.gaussian_alpha = forward.config.rf.gaussian_alpha;
names = {'hvsr', 'ellip', 'cph', 'ugr', 'rf'};
coords = {'f', 'T', 'T', 'T', 't'};
samples = zeros(1,5);
for k = 1:5
    name = names{k};
    assert(isequal(forward.pred.(name).model_input, forward.model), 'Solvers used different Earth models.');
    clean = forward.pred.(name).val(:);
    assert(~isempty(clean) && all(isfinite(clean)), 'Missing/nonfinite predictions.');
    true_model.clean.(name) = clean;
    syn_data.(name).clean = clean;
    syn_data.(name).obs = clean;
    syn_data.(name).sig = zeros(size(clean)); % Forward-only output, not likelihood inputs.
    samples(k) = numel(clean);
    write_csv(fullfile(run_dir, [name '.csv']), [coords{k} ',clean'], ...
        [syn_data.(name).(coords{k}), clean]);
    fprintf('%s: %d finite samples\n', name, samples(k));
end
model = true_model.model_matrix;
p = forward.config.rf.rayp_s_km;
i = true_model.moho_interface;
true_model.moho_primary_ps_delay_s = sum(model(1:i,1) .* ...
    (sqrt(model(1:i,3).^-2-p^2)-sqrt(model(1:i,2).^-2-p^2)));
write_csv(fullfile(run_dir, 'true_model.csv'), 'thickness_km,vp_km_s,vs_km_s,density_g_cm3', model);
write_csv(fullfile(run_dir, 'solver_model.csv'), 'thickness_km,vp_km_s,vs_km_s,density_g_cm3', forward.model);
save(fullfile(run_dir, 'synthetic_thin.mat'), 'true_model', 'syn_data', 'forward', '-v7');
manifest = struct('output_dir', run_dir, 'case_name', 'ecm1_thin', ...
    'all_forward_succeeded', true, 'samples', samples, ...
    'source_layers', size(model,1), 'numerical_layers', size(forward.model,1), ...
    'moho_depth_km', true_model.moho_depth, ...
    'moho_primary_ps_delay_s', true_model.moho_primary_ps_delay_s, ...
    'elapsed_seconds', forward.elapsed_seconds, 'noise_added', false, ...
    'cps_bin', forward.config.cps_bin, 'hvf_bin', forward.config.hvf_bin);
save(fullfile(run_dir, 'manifest.mat'), 'manifest', '-v7');
write_text(fullfile(run_dir, 'manifest.json'), jsonencode(manifest));
fid = fopen(fullfile(run_dir, 'forward_report.txt'), 'w');
if fid < 0, error('fiveobs:io', 'Cannot write report.'); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n%s\n', true_model.description, forward.geometry);
fprintf(fid, 'Truth: prescribed h and Vs; Vp and density derived with original HBI equations.\n');
fprintf(fid, 'All five model-input matrices are identical: %d physical rows / %d numerical rows.\n', size(model,1), size(forward.model,1));
fprintf(fid, 'HVSR 80 samples, 0.2-25 Hz; ellipticity 50, 3-25 s; phase 35, 9-35 s; group 40, 3-25 s.\n');
fprintf(fid, 'P RF: p=%.4g s/km; alpha=%.4g; dt=%.4g s; %d samples; %.4g s pre-arrival.\n', ...
    forward.config.rf.rayp_s_km, forward.config.rf.gaussian_alpha, ...
    forward.config.rf.dt, forward.config.rf.npts, forward.config.rf.delay);
fprintf(fid, 'Moho primary Ps ray-theory lag: %.6f s (not a picked waveform peak).\n', true_model.moho_primary_ps_delay_s);
fprintf(fid, 'No noise, inversion, data fitting, or parameter recovery. Zero sig arrays must not be used as likelihood errors.\n');
fprintf(fid, 'HVf: %s\nCPS: %s\n', forward.config.hvf_bin, forward.config.cps_bin);
[~, basename] = fileparts(run_dir);
write_text(fullfile(cfg.output_root, 'LATEST_RUN.txt'), basename);
fprintf('All five observables completed in %.2f s.\n', forward.elapsed_seconds);
end

function write_csv(filename, header, values)
fid = fopen(filename, 'w');
if fid < 0, error('fiveobs:io', 'Cannot write %s.', filename); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', header);
fprintf(fid, [repmat('%.12g,',1,size(values,2)-1) '%.12g\n'], values');
end

function write_text(filename, value)
fid = fopen(filename, 'w');
if fid < 0, error('fiveobs:io', 'Cannot write %s.', filename); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', value);
end
