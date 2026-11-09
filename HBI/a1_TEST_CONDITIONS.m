function [ifpass, violation] = a1_TEST_CONDITIONS(model, par, ifverbose)
% Check the layered model against its prior conditions, as in THBI a1.
% These are the existing fixed-layer HBI constraints, expressed in explicit
% blocks. Spline, anisotropy and variable-layer-count tests do not apply here.
if nargin < 3, ifverbose = false; end
ifpass = true;
violation = '';
t = model.theta(:);
eta = model.datahparm.eta(:);

%% =================== PARAMETER COUNTS ===================
if numel(t) ~= numel(par.mod.lower) || numel(eta) ~= numel(par.mod.noise.lower)
    ifpass = false;
    violation = 'Structural/noise parameter count does not match the prior.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== FINITE MODEL AND NOISE PARAMETERS ===================
if any(~isfinite(t)) || any(~isfinite(eta))
    ifpass = false;
    violation = 'Model or log-noise parameter is NaN/Inf.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== FINITE-LAYER THICKNESS BOUNDS ===================
ih = par.mod.idx.h;
if any(t(ih) < par.mod.lower(ih) | t(ih) > par.mod.upper(ih))
    ifpass = false;
    violation = 'A finite-layer thickness is outside its prior bounds.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== SHEAR-VELOCITY BOUNDS ===================
iv = par.mod.idx.vs;
if any(t(iv) < par.mod.lower(iv) | t(iv) > par.mod.upper(iv))
    ifpass = false;
    violation = 'A shear velocity is outside its prior bounds.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== MONOTONIC VELOCITY CONDITION ===================
if par.mod.monotonic && any(diff(t(iv)) < 0)
    ifpass = false;
    violation = 'Vs decreases with depth while monotonicity is enabled.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end

%% =================== HIERARCHICAL NOISE BOUNDS ===================
if any(eta < par.mod.noise.lower | eta > par.mod.noise.upper)
    ifpass = false;
    violation = 'A log-noise multiplier is outside its prior bounds.';
    if ifverbose, fprintf('Failed: %s\n',violation); end
    return;
end
end
