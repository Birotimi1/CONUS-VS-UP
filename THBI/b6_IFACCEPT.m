function [accept, logalpha] = b6_IFACCEPT(current, trial, ptb, draw_acceptance)
% Reversible-jump MH at temperature 1. Every ratio term is explicit.
% current/trial contain logL and logP, where logP includes count-dependent
% prior normalizations. q terms are recorded in the proposal from b2.
% draw_acceptance=false evaluates this SAME first-stage rule without RNG,
% including for the hypothetical reverse first stage inside delayed rejection.
if nargin<4, draw_acceptance=true; end
assert(islogical(draw_acceptance) && isscalar(draw_acceptance), ...
    'thbi:accept','draw_acceptance must be a logical scalar.');

%% =================== REJECT INVALID SUPPORT OR SOLVER OUTPUT ===================
accept=false; logalpha=-Inf;
if ~ptb.valid || ~isfinite(trial.logL+trial.logP), return; end

%% =================== LIKELIHOOD, PRIOR, PROPOSAL AND JACOBIAN ===================
delta_logL=trial.logL-current.logL;
delta_logP=trial.logP-current.logP;
log_q_ratio=ptb.log_q_reverse-ptb.log_q_forward;
logalpha=delta_logL+delta_logP+log_q_ratio+ptb.log_abs_jacobian;
assert(isreal(logalpha) && ~isnan(logalpha),'thbi:accept','Undefined MH acceptance ratio.');
if draw_acceptance, accept=log(rand)<min(0,logalpha); end
end
