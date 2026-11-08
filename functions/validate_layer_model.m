function validate_layer_model(model)
% Solid, isotropic layers [h(km), Vp(km/s), Vs(km/s), rho(g/cm^3)].
if ~isnumeric(model) || size(model, 2) ~= 4 || isempty(model) || ...
        ~isreal(model) || any(~isfinite(model(:)))
    error('fiveobs:model', 'Model must be a finite real N-by-4 matrix.');
end
if size(model, 1) > 200
    error('fiveobs:model', 'hrftn96 supports at most 200 layers.');
end
if any(model(1:end-1, 1) <= 0) || model(end, 1) ~= 0
    error('fiveobs:model', 'Finite thicknesses must be positive; last thickness must be 0 (half-space).');
end
if any(model(:, 2) <= 0) || any(model(:, 3) <= 0) || any(model(:, 4) <= 0) || ...
        any(model(:, 2).^2 <= (4/3) * model(:, 3).^2)
    error('fiveobs:model', 'Require positive Vp/Vs/density and positive isotropic bulk modulus.');
end
end
