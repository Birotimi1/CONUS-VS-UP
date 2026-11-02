function trace = read_cps_sac(filename)
% Read a single evenly sampled SAC binary time series, either byte order.
orders = {'ieee-le', 'ieee-be'};
for k = 1:2
    fid = fopen(filename, 'rb', orders{k});
    if fid < 0, error('fiveobs:sac', 'Cannot open %s', filename); end
    floats = fread(fid, 70, 'float32=>double');
    ints = fread(fid, 40, 'int32=>double');
    fread(fid, 192, 'uint8');
    valid = numel(floats) == 70 && numel(ints) == 40;
    if valid
        valid = any(ints(7) == [6 7]) && ints(10) > 0 && ints(10) <= 8192 ...
            && ints(16) == 1 && ints(36) == 1 && floats(1) > 0;
    end
    if valid
        values = fread(fid, ints(10), 'float32=>double');
        fclose(fid);
        if numel(values) ~= ints(10) || any(~isfinite(values))
            error('fiveobs:sac', 'Truncated or nonfinite SAC data: %s', filename);
        end
        trace.val = values;
        trace.dt = floats(1);
        trace.t = floats(6) + (0:ints(10)-1)' * trace.dt;
        trace.gaussian_alpha = floats(41); % USER0
        trace.rayp_s_km = floats(45);       % USER4
        return;
    end
    fclose(fid);
end
error('fiveobs:sac', 'Unsupported or invalid SAC header: %s', filename);
end
