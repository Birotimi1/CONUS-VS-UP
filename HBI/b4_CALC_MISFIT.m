function misfit = b4_CALC_MISFIT(data, predata, par)
% Whiten residuals with fixed C0; cache these energies for noise-only moves.

%% Allocate one energy and one amplitude RMS for each observable
nd = numel(par.inv.datatypes);
misfit.E2 = zeros(nd,1); misfit.rms = zeros(nd,1);
%% Residuals -> covariance whitening -> squared misfit energy
for k = 1:nd
    d = par.inv.datatypes{k}; v = predata.(d)(:);
    if numel(v)~=data.(d).N || any(~isfinite(v))
        error('hbi:prediction','Incomplete prediction for %s.',d);
    end
    r = data.(d).obs-v;            % residual in the observable's own units
    w = data.(d).L0\r;             % dimensionless whitened residual
    misfit.E2(k) = w'*w;
    misfit.rms(k) = sqrt(mean(r.^2));
end
end
