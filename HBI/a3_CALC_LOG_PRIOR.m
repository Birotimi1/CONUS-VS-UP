function logP = a3_CALC_LOG_PRIOR(model, par)
% Density with respect to d(theta) d(eta), eta=ln(gamma).
% Uniform physical bounds + monotone Vs; truncated Gaussian in eta.
% The -ln(gamma) in a lognormal density cancels the transform Jacobian.
%% Reject models outside the structural and log-noise support
if ~a1_TEST_CONDITIONS(model,par)
    logP = -Inf; return;
end
%% Truncated Gaussian prior in eta = ln(gamma)
e = model.datahparm.eta(:);
n = par.mod.noise;
logP = -.5*sum(((e-n.log_center)./n.log_width).^2);
% Fixed normalizing constants cancel in all MH ratios.
end
