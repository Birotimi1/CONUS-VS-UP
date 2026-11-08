function values = sample_layer_model(matrix, depth_km)
% Evaluate [Vp, Vs, rho] on a depth grid without smoothing physical interfaces.
% A point exactly on an interface belongs to the layer immediately below it.
validate_layer_model(matrix);
assert(isnumeric(depth_km) && isreal(depth_km) && isvector(depth_km) && ...
    all(isfinite(depth_km)) && all(depth_km>=0),'thbi:depth','Invalid depth grid.');
interfaces=cumsum(matrix(1:end-1,1));
layer=1+sum(bsxfun(@ge,depth_km(:),interfaces'),2);
values=matrix(layer,2:4);
end
