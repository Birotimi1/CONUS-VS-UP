function [E2,rms_error] = ccp_cross_convolution_misfit(observed, predicted, z, cfg)
% Reference b4: depth-weight daughters, FULL cross convolution, integral RMS.
% E2 is the reference normalized error (not an ordinary squared waveform norm).
assert(isequal(size(observed),size(predicted)) && size(observed,2)==2 && ...
    size(observed,1)==numel(z) && isreal(observed) && isreal(predicted) && ...
    all(isfinite([observed(:);predicted(:)])) && numel(z)>=3 && ...
    all(diff(z)>0) && max(abs(diff(z)-median(diff(z))))<1e-8, ...
    'thbi:ccp','CCP requires matching finite [P SV] component pairs.');
weight=ccp_depth_weights(z,cfg);
observed(:,1)=observed(:,1).*weight; predicted(:,1)=predicted(:,1).*weight;
[~,error,first,second]=xconv_misfit(observed(:,1),observed(:,2),predicted(:,1),predicted(:,2));
dz=median(diff(z)); x=(0:numel(error)-1)'*dz;
integral_rms=@(v) sqrt(trapz(x,v.^2)/(x(end)-x(1)));
denominator=sqrt(integral_rms(first)*integral_rms(second));
assert(isfinite(denominator) && denominator>0,'thbi:ccp','Zero cross-convolution normalization.');
E2=integral_rms(error)/denominator;
rms_error=sqrt(mean((observed(:,1)-predicted(:,1)).^2));
end
