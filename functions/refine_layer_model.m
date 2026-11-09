function refined = refine_layer_model(model, max_thickness_km)
% Change numerical subdivision only, preserving every physical interface.
% Merge exactly equal adjacent materials, including redundant mantle layers
% above an identical half-space; then split thick homogeneous finite layers.
validate_layer_model(model);
if ~isscalar(max_thickness_km) || ~isfinite(max_thickness_km) || max_thickness_km <= 0
    error('fiveobs:layers', 'max_thickness_km must be finite and positive.');
end
merged = model(1, :);
for k = 2:size(model, 1)
    if isequal(model(k, 2:end), merged(end, 2:end))
        if model(k, 1) == 0
            merged(end, 1) = 0;
        else
            merged(end, 1) = merged(end, 1) + model(k, 1);
        end
    else
        merged(end+1, :) = model(k, :); %#ok<AGROW>
    end
end
refined = zeros(0, size(model,2));
for k = 1:size(merged, 1)
    n = max(1, ceil(merged(k, 1) / max_thickness_km));
    rows = repmat(merged(k, :), n, 1);
    rows(:, 1) = merged(k, 1) / n;
    refined = [refined; rows]; %#ok<AGROW>
end
validate_layer_model(refined);
end
