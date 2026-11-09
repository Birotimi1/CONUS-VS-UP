function [accept, logalpha] = delayed_rejection_damp(current, first, second, ptb)
% Correct two-stage delayed rejection for ONE unbounded Gaussian coordinate.
% Adapts the HBI second-try role, not its likelihood-only acceptance equation.
% x -> y rejected; z is proposed from x using a smaller Gaussian on the SAME
% selected coordinate. All three models have the same dimension. The smaller
% q2 is symmetric and independent of y, so its forward/reverse terms cancel.
% q1(y|z)/q1(y|x) and BOTH first-stage rejection factors must remain.
% Green & Mira (2001), doi:10.1093/biomet/88.4.1035.
accept=false; logalpha=-Inf;
assert(ptb.dr_eligible && ~any(strcmp(ptb.type,{'birth','death'})), ...
    'thbi:dr','Delayed rejection here is restricted to a fixed scalar coordinate.');
lx=current.logL+current.logP; ly=first.logL+first.logP; lz=second.logL+second.logP;
if ~isfinite(lx) || ~isfinite(lz) || isnan(ly), return; end
ax=min(0,ly-lx); az=min(0,ly-lz);
if ax>=0 || az>=0, return; end % reverse first stage would accept with certainty
x=ptb.value_from; y=ptb.value_to; z=second.ptb.value_to;
assert(isequal(current.model.K,second.model.K) && isequal(current.model.K,first.model.K) && ...
    strcmp(ptb.type,second.ptb.type) && ptb.zone==second.ptb.zone && ...
    ptb.index==second.ptb.index && ptb.value_from==second.ptb.value_from && ...
    ptb.sigma>0 && isfinite(ptb.sigma) && second.ptb.sigma>0 && ...
    second.ptb.sigma<ptb.sigma && all(isfinite([x y z])), ...
    'thbi:dr','Second-stage dimension/scale is invalid.');
% Factoring the difference avoids subtracting two large squared distances.
log_q1_ratio=.5*((z-x)/ptb.sigma)*((y-z)/ptb.sigma+(y-x)/ptb.sigma);
% expm1 is stable when a first-stage acceptance probability is nearly one.
logalpha=lz-lx+log_q1_ratio+log(-expm1(az))-log(-expm1(ax));
assert(isreal(logalpha) && ~isnan(logalpha),'thbi:dr','Undefined delayed-rejection acceptance ratio.');
accept=log(rand)<min(0,logalpha);
end
