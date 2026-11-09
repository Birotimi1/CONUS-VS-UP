function [logP, terms] = a3_CALC_LOG_PRIOR(model, par)
% Density with respect to zone-boundary depths, ordered interface FRACTIONS,
% layer Vs, half-space Vs, and eta=ln(gamma), summed over discrete counts K.
% Dimension-dependent normalizations MUST remain in an RJ target.
terms = struct('counts',0,'interfaces',0,'velocity',0,'boundaries',0,'noise',0);
if ~a1_TEST_CONDITIONS(model,par), logP=-Inf; return; end

%% =================== DISCRETE LAYER-COUNT MASS ===================
for z = 1:3
    K = numel(model.(par.mod.zone_fields{z}).VS);
    counts = par.mod.kmin(z):par.mod.kmax(z);
    switch par.mod.count_prior
        case 'jeffreys'
            terms.counts = terms.counts-log(K)-log(sum(1./counts));
        case 'uniform'
            terms.counts = terms.counts-log(numel(counts));
        otherwise
            error('thbi:prior','Unknown layer-count prior.');
    end

    %% ORDERED INTERFACES WITH A MINIMUM FRACTIONAL THICKNESS
    % f_i-i*epsilon are K-1 ordered points in [0,1-K*epsilon].
    S = 1-K*par.mod.min_fraction(z);
    terms.interfaces = terms.interfaces+gammaln(K)-(K-1)*log(S);

    %% INDEPENDENT OR ORDERED UNIFORM VELOCITIES WITHIN THIS ZONE
    width = par.mod.vsmax(z)-par.mod.vsmin(z);
    terms.velocity = terms.velocity-K*log(width);
    if par.conditions.monotonic(z)
        terms.velocity = terms.velocity+gammaln(K+1); % ordered-uniform K!
    end
end

%% =================== BOUNDARIES AND HALF-SPACE VELOCITY ===================
terms.boundaries = -sum(log(par.mod.boundary_upper-par.mod.boundary_lower));
terms.velocity = terms.velocity-log(diff(par.mod.halfspace_bounds));

%% =================== TRUNCATED NORMAL IN LOG-NOISE COORDINATES ===================
n = par.mod.noise; eta = model.datahparm.eta(:);
PhiUpper = .5*(1+erf((n.upper-n.log_center)./(sqrt(2)*n.log_width)));
PhiLower = .5*(1+erf((n.lower-n.log_center)./(sqrt(2)*n.log_width)));
terms.noise = sum(-.5*((eta-n.log_center)./n.log_width).^2 ...
    -log(n.log_width)-.5*log(2*pi)-log(PhiUpper-PhiLower));
logP = terms.counts+terms.interfaces+terms.velocity+terms.boundaries+terms.noise;
% Cross-zone conditions in a1 multiply this joint density by an indicator.
% Their overall normalization cancels in MH; their induced K marginal is
% not asserted to equal the base count mass above.
end
