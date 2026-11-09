function weight = ccp_depth_weights(z, cfg)
% Zach/Brennan Sp depth weighting, Gaussian smoothing over 80 km, RMS=1.
z=z(:); knots=cfg.weight_depth_val;
weight=interp1(knots(:,1),knots(:,2),z,'linear');
assert(all(isfinite(weight) & weight>0),'thbi:ccp','Depth weights do not cover CCP grid.');
% Explicit Gaussian moving window; truncate and renormalize at the ends.
% MATLAB smoothdata uses a Gaussian window with std=window_length/5.
dz=median(diff(z)); n=max(1,round(cfg.weight_smooth_km/dz)+1);
if exist('smoothdata','file') || exist('smoothdata','builtin')
    weight=smoothdata(weight,'gaussian',n);
    weight=weight/sqrt(mean(weight.^2));
    return;
end
% Octave solver-free fallback; native MATLAB uses the reference call above.
offset=(1:n)'-((n+1)/2); kernel=exp(-.5*(offset/(n/5)).^2);
weight=conv(weight,kernel,'same')./conv(ones(size(weight)),kernel,'same');
weight=weight/sqrt(mean(weight.^2));
end
