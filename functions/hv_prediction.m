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
