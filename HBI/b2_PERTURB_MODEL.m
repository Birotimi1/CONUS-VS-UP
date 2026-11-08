function [trial, move] = b2_PERTURB_MODEL(model, par, steps, j)
% Random-scan, one-coordinate Gaussian proposals. Reflection is symmetric.
% Adapted from HBI c5_perturb_theta; noise-only branch follows THBI separation.

%% Select a structural parameter or a log-noise multiplier
nt = numel(model.theta); nd = numel(par.inv.datatypes);
if nargin < 4
    if rand < par.inv.noise_move_probability, j = nt+randi(nd);
    else, j = randi(nt); end
end
move.index = j; move.structure = j <= nt;
t = model.theta; e = model.datahparm.eta;
if move.structure
    %% Structural proposal: bounds, monotonic neighbors, Gaussian reflection
    lo = par.mod.lower(j); hi = par.mod.upper(j);
    loc = find(par.mod.idx.vs == j,1);
    if par.mod.monotonic && ~isempty(loc)
        if loc>1, lo=max(lo,t(par.mod.idx.vs(loc-1))); end
        if loc<numel(par.mod.idx.vs), hi=min(hi,t(par.mod.idx.vs(loc+1))); end
    end
    t(j) = reflect_to_interval(t(j)+steps(j)*randn,lo,hi);
    % Rebuild Vp/density using the actual equations in b1, without jitter.
    trial = b1_INITIATE_MODEL(par,t,e);
else
    %% Noise-only proposal: keep the Earth unchanged and update eta/gamma
    q = j-nt;
    e(q) = reflect_to_interval(e(q)+steps(j)*randn, ...
        par.mod.noise.lower(q),par.mod.noise.upper(q));
    trial = model; trial.datahparm.eta = e; trial.datahparm.gamma = exp(e);
end
end
