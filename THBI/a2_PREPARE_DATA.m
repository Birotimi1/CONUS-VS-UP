function data = a2_PREPARE_DATA(data, par)
% Fix coordinates, uncertainties, and covariance BEFORE starting any chain.
% No proposal-dependent masks, truncation, normalization, or dropped data.

%% Validate each dataset and build its fixed whitening matrix
for k = 1:numel(par.inv.datatypes)
    d = par.inv.datatypes{k};
    if ~isfield(data,d), error('hbi:data','Missing %s data.',d); end
    if isfield(par,'rf') && strcmp(par.rf.method,'hk_ccp') && k>=5
        v=data.(d);
        assert(v.N==numel(v.x) && v.N==numel(v.obs) && ...
            isreal(v.obs) && all(isfinite(v.obs)) && v.dof>0 && v.nrecords==1, ...
            'hbi:data','Invalid reference RF data product.');
        if k==5
            w=v.waves;
            assert(v.N==1 && size(w.rf,1)==numel(w.t) && size(w.rf,2)==numel(w.rayp_s_km) && ...
                isreal(w.rf) && all(isfinite(w.rf(:))) && all(isfinite(w.t)) && ...
                all(diff(w.t)>0) && all(isfinite(w.rayp_s_km) & w.rayp_s_km>0), ...
                'hbi:data','Invalid observed Ps waves for H-kappa.');
        else
            assert(isequal(v.x,par.rf.ccp.depth_km) && numel(v.parent)==v.N && ...
                isreal(v.parent) && all(isfinite(v.parent)),'hbi:data','Invalid Sp depth grid/parent.');
        end
        continue;
    end
    v = data.(d); v.x = v.x(:); v.obs = v.obs(:); v.sig = v.sig(:);
    N = numel(v.obs);
    if ~isreal(v.x) || ~isreal(v.obs) || ~isreal(v.sig) || ...
            N < 2 || numel(v.x) ~= N || numel(v.sig) ~= N || ...
            any(~isfinite([v.x;v.obs;v.sig])) || any(v.sig <= 0) || ...
            any(diff(v.x) <= 0)
        error('hbi:data','%s requires matching finite arrays, positive sig, increasing x.',d);
    end
    if ~any(strcmp(d,{'rf','rf_sp'})) && any(v.x <= 0)
        error('hbi:data','Positive frequency/period required.');
    end
    if ~isfield(v,'R'), v.R = eye(N); end
    if ~isreal(v.R) || ~isequal(size(v.R),[N N]) || any(~isfinite(v.R(:))) || ...
            norm(v.R-v.R','fro') > 1e-10 || any(abs(diag(v.R)-1)>1e-10)
        error('hbi:covariance','%s R must be a finite symmetric correlation matrix.',d);
    end
    % C0 = diag(sig)*R*diag(sig); L0*L0' = C0.
    [L,flag] = chol(v.R,'lower');
    if flag, error('hbi:covariance','%s R must be positive definite.',d); end
    v.N = N;
    v.L0 = bsxfun(@times,v.sig,L); % C0 = diag(sig)*R*diag(sig)
    v.logdetC0 = 2*sum(log(diag(v.L0)));
    assert(all(isfinite(v.L0(:))) && all(diag(v.L0)>0) && isfinite(v.logdetC0), ...
        'hbi:covariance','%s covariance scaling is not numerically representable.',d);
    data.(d) = v;
end
end
