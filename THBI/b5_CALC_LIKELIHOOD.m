function [logL, misfit] = b5_CALC_LIKELIHOOD(misfit, data, datahparm, par)
% Full Gaussian log likelihood, summed over the configured independent data blocks.
% C_d = exp(2 eta_d) C0_d. Retain the N*eta normalization when inferring noise.

%% Evaluate each normalized Gaussian likelihood with its current gamma
nd = numel(par.inv.datatypes);
assert(numel(datahparm.eta)==nd && isreal(datahparm.eta) && all(isfinite(datahparm.eta)) && ...
    numel(misfit.E2)==nd && isreal(misfit.E2) && all(isfinite(misfit.E2) & misfit.E2>=0), ...
    'hbi:likelihood','Invalid log-noise parameters or cached residual energies.');
misfit.logL_indivdat = zeros(nd,1); misfit.chi2 = zeros(nd,1);
for k = 1:nd
    d = par.inv.datatypes{k}; e = datahparm.eta(k);
    misfit.chi2(k) = misfit.E2(k)*exp(-2*e);
    misfit.logL_indivdat(k) = -.5*(data.(d).N*log(2*pi) + ...
        data.(d).logdetC0 + 2*data.(d).N*e + misfit.chi2(k));
end
%% Joint likelihood: independent data types in this controlled synthetic
logL = sum(misfit.logL_indivdat);
assert(isfinite(logL),'hbi:likelihood','Nonfinite likelihood; inspect errors, covariance and noise scales.');
end
