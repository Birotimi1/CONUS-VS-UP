function [ifpass, violation] = a1_TEST_CONDITIONS(model, par, ifverbose)
% Explicit validity blocks for the zoned, variable-layer model.
if nargin < 3, ifverbose = false; end
ifpass = false; violation = '';

%% =================== SEDIMENT BASE, MOHO AND HALF-SPACE ===================
if ~isscalar(model.zsed) || ~isscalar(model.zmoh) || ...
        ~isreal(model.zsed) || ~isreal(model.zmoh)
    violation = 'Zone boundaries must be real scalars.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end
b = [model.zsed model.zmoh];
if any(~isfinite(b)) || any(b<par.mod.boundary_lower | b>par.mod.boundary_upper) || ...
        model.zsed<=0 || model.zmoh<=model.zsed || model.zmoh>=par.mod.maxz
    violation = 'Zone boundaries are outside their ordered prior bounds.';
elseif ~isscalar(model.halfspaceVS) || ~isreal(model.halfspaceVS) || ~isfinite(model.halfspaceVS) || ...
        model.halfspaceVS<par.mod.halfspace_bounds(1) || model.halfspaceVS>par.mod.halfspace_bounds(2)
    violation = 'Half-space Vs is outside its prior bounds.';
end
if ~isempty(violation)
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== COUNTS, INTERFACES AND VELOCITIES IN EACH ZONE ===================
for z = 1:3
    zone = model.(par.mod.zone_fields{z});
    v = zone.VS(:); f = zone.f(:); K = numel(v);
    if K<par.mod.kmin(z) || K>par.mod.kmax(z) || numel(f)~=K-1
        violation = sprintf('%s layer/interface counts are invalid.',par.mod.zone_names{z});
    elseif ~isreal(v) || ~isreal(f) || any(~isfinite([v;f]))
        violation = sprintf('%s has nonfinite or complex parameters.',par.mod.zone_names{z});
    elseif any(diff([0;f;1])<par.mod.min_fraction(z))
        violation = sprintf('%s interfaces violate ordering/minimum separation.',par.mod.zone_names{z});
    elseif any(v<par.mod.vsmin(z) | v>par.mod.vsmax(z))
        violation = sprintf('%s Vs is outside its prior bounds.',par.mod.zone_names{z});
    elseif par.conditions.monotonic(z) && any(diff(v)<0)
        violation = sprintf('%s Vs decreases while monotonicity is enabled.',par.mod.zone_names{z});
    end
    if ~isempty(violation)
        if ifverbose, fprintf('Failed: %s\n',violation); end
        return;
    end
end

%% =================== BASEMENT AND MOHO VELOCITY CONTRASTS ===================
moho_pair=[model.crustmparm.VS(end),model.mantmparm.VS(1)];
moho_jump=diff(moho_pair);
moho_fraction=moho_jump/mean(moho_pair);
if par.conditions.pos_sed2basement && model.crustmparm.VS(1)<model.sedmparm.VS(end)
    violation = 'Negative sediment-to-basement velocity jump.';
elseif par.conditions.pos_moho && model.mantmparm.VS(1)<model.crustmparm.VS(end)
    violation = 'Negative Moho velocity jump.';
elseif moho_fraction > par.conditions.max_moho_fraction
    violation = 'Moho velocity increase exceeds its allowed fraction.';
elseif moho_jump > par.conditions.max_moho_jump
    violation = 'Moho velocity increase exceeds its allowed absolute jump.';
end
if ~isempty(violation)
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== HIERARCHICAL NOISE MULTIPLIERS ===================
eta = model.datahparm.eta(:);
if numel(eta)~=numel(par.inv.datatypes) || ~isreal(eta) || any(~isfinite(eta)) || ...
        any(eta<par.mod.noise.lower | eta>par.mod.noise.upper)
    violation = 'Log-noise parameters are invalid or outside their bounds.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end
ifpass = true;
end
