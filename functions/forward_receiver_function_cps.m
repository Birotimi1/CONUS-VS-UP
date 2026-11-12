function rf = forward_receiver_function_cps(model, cfg, cps_bin, workdir)
% Synthetic P receiver function from the shared layered model.
% hrftn96 computes the radial/vertical frequency-domain response directly,
% including multiples, then applies its Gaussian filter. No observed RF,
% arbitrary amplitude normalization, or fabricated RF misfit is used.
validate_layer_model(model);
rf_program='hrftn96';
if size(model,2)==7, rf_program='trftn96'; end
rf_sac=[rf_program '.sac'];
vti_vertical_slowness(model,cfg.rayp_s_km);
if ~isscalar(cfg.rayp_s_km) || ~isfinite(cfg.rayp_s_km) || cfg.rayp_s_km <= 0
    error('fiveobs:rayp', 'Require 0 < rayp < 1/Vp in every layer (rayp in s/km).');
end
if ~isscalar(cfg.npts) || cfg.npts < 64 || cfg.npts > 4096 ...
        || log2(cfg.npts) ~= floor(log2(cfg.npts))
    error('fiveobs:rf', 'npts must be a power of two, 64..4096; -2 uses twice this length.');
end
if ~isscalar(cfg.dt) || ~isfinite(cfg.dt) || cfg.dt <= 0 || ...
        ~isscalar(cfg.delay) || ~isfinite(cfg.delay) || cfg.delay <= 0 || ...
        cfg.delay >= cfg.npts * cfg.dt || ...
        ~isscalar(cfg.gaussian_alpha) || ~isfinite(cfg.gaussian_alpha) || cfg.gaussian_alpha <= 0
    error('fiveobs:rf', 'Invalid RF sample interval, delay, or Gaussian alpha.');
end
if ~exist(workdir, 'dir'), mkdir(workdir); end
previous = pwd;
cleanup = onCleanup(@() cd(previous)); %#ok<NASGU>
cd(workdir);
write_elastic_cps_model('rf_model.mod', model);
cmd = sprintf('%s -P -2 -M rf_model.mod -RAYP %.10g -ALP %.10g -DT %.10g -NSAMP %d -D %.10g', ...
    quote_shell_path(fullfile(cps_bin,rf_program)), cfg.rayp_s_km, ...
    cfg.gaussian_alpha, cfg.dt, cfg.npts, cfg.delay);
if exist(rf_sac, 'file'), delete(rf_sac); end
[status, message] = run_native_command([cmd ' > hrftn96.log 2>&1'],cfg,'Ps/CPS');
if status ~= 0 || ~exist(rf_sac, 'file')
    error('fiveobs:rf', 'hrftn96 failed (%d): %s. See %s', status, message, fullfile(workdir, 'hrftn96.log'));
end
rf = read_cps_sac(rf_sac);
if numel(rf.val) ~= cfg.npts || abs(rf.dt - cfg.dt) > 1e-6 || ...
        abs(rf.t(1) + cfg.delay) > 1e-5 || ...
        abs(rf.rayp_s_km - cfg.rayp_s_km) > 1e-6 || ...
        abs(rf.gaussian_alpha - cfg.gaussian_alpha) > 1e-5
    error('fiveobs:rf', 'RF SAC header does not match the requested settings.');
end
rf.kind = 'Synthetic P-to-S RF: Gaussian-filtered radial/vertical transfer response';
rf.model_input = model;
rf.command = cmd;
rf.observed = false;
rf.status = 'ok';
% Kinematic primary Ps delays to each interface, not full waveform picks.
[qp,qs]=vti_vertical_slowness(model,cfg.rayp_s_km);
qp=qp(1:end-1); qs=qs(1:end-1);
rf.interface_depth_km = cumsum(model(1:end-1, 1));
rf.primary_ps_delay_s = cumsum(model(1:end-1, 1) .* (qs - qp));
rf.interface_contrast = any(abs(diff(model(:, 2:end), 1, 1)) > 1e-10, 2);
end
