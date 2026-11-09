function model = b1_INITIATE_MODEL(par, theta, eta)
% Build the physical layered Earth and initialize a chain when requested.
% b1_INITIATE_MODEL(par): use the reference start and independently jitter it.
% b1_INITIATE_MODEL(par,theta,eta): build the supplied model WITHOUT jitter.
% b2 uses the second form so Vp/density have a single implementation here.

%% Choose the structural parameters and hierarchical noise parameters
initialize_chain = nargin < 2;
if initialize_chain
    theta = par.mod.start;
end
if nargin < 3
    eta = par.mod.noise.log_center;
end
theta = theta(:);
eta = eta(:);

%% =================== INDEPENDENT STARTING-MODEL PERTURBATIONS ===================
if initialize_chain
    initial.theta = theta;
    initial.datahparm.eta = eta;
    [ifpass,violation] = a1_TEST_CONDITIONS(initial,par);
    assert(ifpass,'Initial reference model violates the prior: %s',violation);
    steps = par.mod.step*par.mod.start_jitter;
    order = randperm(numel(theta));
    for j = order
        % Conditional bounds are fixed by the unchanged neighboring Vs values.
        lower = par.mod.lower(j);
        upper = par.mod.upper(j);
        loc = find(par.mod.idx.vs == j,1);
        if par.mod.monotonic && ~isempty(loc)
            if loc > 1
                lower = max(lower,theta(par.mod.idx.vs(loc-1)));
            end
            if loc < numel(par.mod.idx.vs)
                upper = min(upper,theta(par.mod.idx.vs(loc+1)));
            end
        end

        % Gaussian proposal, then repeated reflection at the two boundaries.
        proposed = theta(j)+steps(j)*randn;
        if lower == upper
            theta(j) = lower;
        else
            assert(isfinite(lower) && isfinite(upper) && lower < upper, ...
                'Invalid starting-model proposal interval.');
            width = upper-lower;
            distance = mod(proposed-lower,2*width);
            if distance > width
                distance = 2*width-distance;
            end
            theta(j) = min(upper,max(lower,lower+distance));
        end
    end

end

%% =================== FINAL LAYER PROPERTIES AND MODEL STRUCTURE ===================
h = theta(par.mod.idx.h);             % finite-layer thicknesses, km
vs = theta(par.mod.idx.vs);           % Vs including the half-space, km/s
%% Derive Vp and density using the original HBI equations and floor order
vp = .9409 + 2.0947*vs - .8206*vs.^2 + .2683*vs.^3 - .0251*vs.^4;
rho = 1.6612*vp - .4721*vp.^2 + .0671*vp.^3 - .0043*vp.^4 + .000106*vp.^5;
vp = max(vp,1.5*vs);
rho = max(rho,1.5);

%% Assemble the shared model matrix and data-noise parameters
model.theta = theta;
model.matrix = [[h;0],vp,vs,rho];      % [h, Vp, Vs, rho]; last h=0 is half-space
model.datahparm.eta = eta;
model.datahparm.gamma = exp(eta);

end
