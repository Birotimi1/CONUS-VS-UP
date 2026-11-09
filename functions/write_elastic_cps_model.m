function write_elastic_cps_model(filename, model)
% Same elastic, flat-Earth structure used by HV-DFA and hrftn96.
validate_layer_model(model);
fid = fopen(filename, 'w');
if fid < 0, error('fiveobs:io', 'Cannot write %s', filename); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
if size(model,2)==7
    e=vti_elastic(model);
    assert(all(e.F>=0),'fiveobs:model','CPS model96 VPF requires nonnegative F.');
    fprintf(fid, ['MODEL.01\nShared VTI five-observable model\nTRANSVERSE ISOTROPIC\nKGS\n' ...
        'FLAT EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\n' ...
        'H VPV VSV RHO QP QS ETAP ETAS FREFP FREFS\nVPH VSH VPF\n']);
    for k=1:size(model,1)
        fprintf(fid,'%.12g %.12g %.12g %.12g 0 0 0 0 1 1\n', ...
            model(k,1),e.Vpv(k),e.Vsv(k),e.rho(k));
        fprintf(fid,'%.12g %.12g %.12g\n',e.Vph(k),e.Vsh(k),sqrt(e.F(k)/e.rho(k)));
    end
    return;
end
fprintf(fid, ['MODEL.01\nShared five-observable model\nISOTROPIC\nKGS\n' ...
    'FLAT EARTH\n1-D\nCONSTANT VELOCITY\nLINE08\nLINE09\nLINE10\nLINE11\n' ...
    'H VP VS RHO QP QS ETAP ETAS FREFP FREFS\n']);
for k = 1:size(model, 1)
    fprintf(fid, '%.12g %.12g %.12g %.12g 0 0 0 0 1 1\n', model(k, :));
end
end
