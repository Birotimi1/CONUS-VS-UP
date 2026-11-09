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
if isfield(par.mod,'vpvs_bounds') && (~isfield(model,'vpvs') || ...
        ~isscalar(model.vpvs) || ~isreal(model.vpvs) || ~isfinite(model.vpvs) || ...
        model.vpvs<par.mod.vpvs_bounds(1) || model.vpvs>par.mod.vpvs_bounds(2))
    violation='Crustal Vp/Vs is outside its uniform prior.'; return;
end
% Layered counterpart of reference near-Moho gradient constraints: secants
% over a FIXED physical distance on each side, excluding the Moho jump.
% This permits LVZs elsewhere; do not differentiate delta-function interfaces.
if isfield(par.conditions,'moho_gradient_window_km')
    w=par.conditions.moho_gradient_window_km;
    edges=[0 model.zsed model.zmoh par.mod.maxz]; gradients=zeros(1,2);
    for z=2:3
        zone=model.(par.mod.zone_fields{z}); width=min(w,diff(edges(z:z+1)));
        f=zone.f(:); v=zone.VS(:);
        if z==2
            j=1+sum(f<1-width/diff(edges(z:z+1)));
            gradients(1)=(v(end)-v(j))/width;
        else
            j=1+sum(f<width/diff(edges(z:z+1)));
            gradients(2)=(v(j)-v(1))/width;
        end
    end
    limits=par.conditions.moho_gradient_bounds;
    if any(gradients<limits(1) | gradients>limits(2))
        violation='Near-Moho layered secant gradient exceeds its bounds.'; return;
    end
end
eta = model.datahparm.eta(:);
if numel(eta)~=numel(par.inv.datatypes) || ~isreal(eta) || any(~isfinite(eta)) || ...
        any(eta<par.mod.noise.lower | eta>par.mod.noise.upper)
    violation = 'Log-noise parameters are invalid or outside their bounds.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end
if isfield(par.mod,'xi_bounds')
    if ~isfield(model,'xi') || ~isequal(size(model.xi),[1 3]) || ...
            ~isreal(model.xi) || any(~isfinite(model.xi)) || ...
            any(model.xi'<par.mod.xi_bounds(:,1) | model.xi'>par.mod.xi_bounds(:,2))
        violation='Zone xi is outside its prior (including fixed values).'; return;
    end
    % Tensor stability depends jointly on xi, kappa and Vs. Test the proposed
    % parameters, never a stale pre-proposal matrix; this is part of support.
    try
        make_mod_from_parms(model,par);
    catch problem
        if ~strcmp(problem.identifier,'fiveobs:model'), rethrow(problem); end
        violation=problem.message; return;
    end
end
ifpass = true;
end
