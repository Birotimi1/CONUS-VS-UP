function model = make_mod_from_parms(model, par)
% THBI-style model assembly for piecewise-constant layers instead of splines.
% f contains K-1 internal fractional depths; VS has K physical layer values.
boundaries = [0 model.zsed model.zmoh par.mod.maxz];
h = []; vs = []; model.K = zeros(1,3);
for z = 1:3
    field = par.mod.zone_fields{z};
    f = model.(field).f(:); v = model.(field).VS(:);
    model.(field).f = f; model.(field).VS = v;
    model.(field).Nlay = numel(v); model.K(z) = numel(v);
    h = [h;diff([0;f;1])*(boundaries(z+1)-boundaries(z))]; %#ok<AGROW>
    vs = [vs;v]; %#ok<AGROW>
end
vs = [vs;model.halfspaceVS];
vp = .9409 + 2.0947*vs - .8206*vs.^2 + .2683*vs.^3 - .0251*vs.^4;
rho = 1.6612*vp - .4721*vp.^2 + .0671*vp.^3 - .0043*vp.^4 + .000106*vp.^5;
vp = max(vp,1.5*vs); rho = max(rho,1.5);
model.matrix = [[h;0],vp,vs,rho];
model.VS = vs; model.VP = vp; model.rho = rho;
model.zlayt = [0;cumsum(h)];
model.zmax = par.mod.maxz;
moho_pair=[model.crustmparm.VS(end),model.mantmparm.VS(1)];
model.dVSmoh=diff(moho_pair);                       % km/s
model.fdVSmoh=100*model.dVSmoh/mean(moho_pair);     % percent, THBI convention
model.M = 2*sum(model.K)+numel(par.inv.datatypes);
model.datahparm.eta = model.datahparm.eta(:);
model.datahparm.gamma = exp(model.datahparm.eta);
end
