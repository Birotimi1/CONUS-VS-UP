function [logL, misfit] = b5_CALC_LIKELIHOOD(misfit, data, datahparm, par)
% Full Gaussian log likelihood, summed over five independent data types.
% C_d = exp(2 eta_d) C0_d. Retain the N*eta normalization when inferring noise.

%% Evaluate each normalized Gaussian likelihood with its current gamma
nd = numel(par.inv.datatypes);
misfit.logL_indivdat = zeros(nd,1); misfit.chi2 = zeros(nd,1);
for k = 1:nd
    d = par.inv.datatypes{k}; e = datahparm.eta(k);
    misfit.chi2(k) = misfit.E2(k)*exp(-2*e);
    misfit.logL_indivdat(k) = -.5*(data.(d).N*log(2*pi) + ...
        data.(d).logdetC0 + 2*data.(d).N*e + misfit.chi2(k));
end
%% Joint likelihood: independent data types in this controlled synthetic
logL = sum(misfit.logL_indivdat);
end
