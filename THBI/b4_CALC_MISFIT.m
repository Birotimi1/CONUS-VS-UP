function misfit = b4_CALC_MISFIT(data, predata, par)
% Whiten residuals with fixed C0; cache these energies for noise-only moves.

%% Allocate one energy and one amplitude RMS for each observable
nd = numel(par.inv.datatypes);
misfit.E2 = zeros(nd,1); misfit.rms = zeros(nd,1);
%% Residuals -> covariance whitening -> squared misfit energy
for k = 1:nd
    d = par.inv.datatypes{k}; v = predata.(d)(:);
    if ~isreal(v) || numel(v)~=data.(d).N || any(~isfinite(v))
        error('hbi:prediction','Incomplete prediction for %s.',d);
    end
    if isfield(par,'rf') && strcmp(par.rf.method,'hk_ccp') && strcmp(d,'rf')
        % Reference hk_forward_model: normalized energy deficit + error floor.
        % It is not a covariance-whitened waveform squared residual.
        assert(v<=1+1e-10,'thbi:hk','Stack energy exceeds its maximum.');
        misfit.E2(k)=par.rf.hk.error_scale*(1-v+par.rf.hk.min_error);
        misfit.rms(k)=abs(1-v);
    elseif isfield(par,'rf') && strcmp(par.rf.method,'hk_ccp') && strcmp(d,'rf_sp')
        [misfit.E2(k),misfit.rms(k)]=ccp_cross_convolution_misfit( ...
            [data.(d).obs,data.(d).parent],[v,predata.rf_sp_parent],data.(d).x,par.rf.ccp);
    else
        r = data.(d).obs-v;        % residual in the observable's own units
        w = data.(d).L0\r;         % dimensionless whitened residual
        misfit.E2(k) = w'*w;
        misfit.rms(k) = sqrt(mean(r.^2));
    end
    assert(isfinite(misfit.E2(k)) && isfinite(misfit.rms(k)), ...
        'hbi:misfit','Nonfinite residual energy for %s; check data and covariance scale.',d);
end
end
