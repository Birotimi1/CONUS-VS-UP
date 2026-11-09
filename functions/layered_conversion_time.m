function lag = layered_conversion_time(matrix, depth_km, rayp_s_km, phase)
% Exact flat-layer vertical-slowness integral; half-space extends below base.
validate_layer_model(matrix); z=depth_km(:);
assert(all(isfinite(z)) && isscalar(rayp_s_km) && rayp_s_km>0 && ...
    all(rayp_s_km*matrix(:,2)<1),'thbi:migration','Invalid depths/ray parameter.');
h=matrix(:,1); tops=[0;cumsum(h(1:end-1))]; h(end)=Inf;
lengths=max(0,min(bsxfun(@minus,abs(z),tops'),h'));
q=sqrt(matrix(:,3).^-2-rayp_s_km^2)-sqrt(matrix(:,2).^-2-rayp_s_km^2);
lag=lengths*q; lag(z<0)=-lag(z<0);
if strcmp(phase,'Sp'), lag=-lag;
else, assert(strcmp(phase,'Ps'),'thbi:migration','Choose Ps or Sp.'); end
end
