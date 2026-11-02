function write_elastic_cps_model(filename, model)
% Same elastic, flat-Earth structure used by HV-DFA and hrftn96.
validate_layer_model(model);
fid = fopen(filename, 'w');
if fid < 0, error('fiveobs:io', 'Cannot write %s', filename); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, ['MODEL.01\nShared five-observable model\nISOTROPIC\nKGS\n' ...
    'FLAT EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\n' ...
    'H VP VS RHO QP QS ETAP ETAS FREFP FREFS\n']);
for k = 1:size(model, 1)
    fprintf(fid, '%.12g %.12g %.12g %.12g 0 0 0 0 1 1\n', model(k, :));
end
end
