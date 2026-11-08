function [accept, logalpha, terms] = delayed_rejection_layered(current, first, second, ptb, par)
% Full two-stage RJ path: x -> y (rejected) -> z, with second proposal FROM y.
% Q1 is the ordinary b2 kernel. Q2 may target the first move's zone, retaining
% all-zone support. The reverse path is z -> y (rejected) -> x.
% Both alpha1 evaluations call b6_IFACCEPT(...,false): no duplicated first-
% stage rule and no random numbers consumed by hypothetical reverse tests.
% See docs/THBI_TARGETED_DR.md for the augmented-path Jacobian and derivation.
accept=false; logalpha=-Inf; terms=struct();
mode=par.proposal.dr.mode;
assert(par.proposal.dr.enabled && any(strcmp(mode,{'layered_path','targeted_layered'})), ...
    'thbi:dr','Layered DR requires a supported full-path proposal mode.');
if ~ptb.valid || ~second.ptb.valid, return; end
p2=second.ptb;
assert(ptb.step_scale==1 && p2.step_scale==1 && ...
    strcmp(ptb.dr_mode,mode) && strcmp(p2.dr_mode,mode) && ...
    isnan(ptb.previous_zone) && ...
    isequal(ptb.K_before,current.model.K) && isequal(ptb.K_after,first.model.K) && ...
    isequal(p2.K_before,first.model.K) && isequal(p2.K_after,second.model.K), ...
    'thbi:dr','Second proposal must start at the rejected Earth; first proposal must use Q1.');

%% =================== FINITE TARGETS AND BASE RJ PROPOSAL DENSITIES ===================
lx=current.logL+current.logP; ly=first.logL+first.logP; lz=second.logL+second.logP;
assert(isreal([lx ly lz]) && ~any(isnan([lx ly lz])), ...
    'thbi:dr','Undefined path target density.');
if ~all(isfinite([lx ly lz])), return; end
values=[ptb.log_q_forward ptb.log_q_reverse ptb.log_abs_jacobian ...
    p2.log_q_forward p2.log_q_reverse p2.log_abs_jacobian];
assert(isreal(values) && all(isfinite(values)),'thbi:dr','Undefined path proposal density.');

%% =================== SAME FIRST-STAGE RULE IN BOTH DIRECTIONS ===================
[~,loga_xy]=b6_IFACCEPT(current,first,ptb,false);
reverse=p2;
reverse.log_q_forward=p2.log_q_reverse;
reverse.log_q_reverse=p2.log_q_forward;
reverse.log_abs_jacobian=-p2.log_abs_jacobian;
[~,loga_zy]=b6_IFACCEPT(second,first,reverse,false);
terms.logalpha_xy=min(0,loga_xy); terms.logalpha_zy=min(0,loga_zy);

%% =================== ACTUAL Q2: TARGET ZONE AND REVERSE TARGET ZONE ===================
% Q2(y->z | x,y) / Q1(y->z) is a zone-weight ratio only. Scalar step
% distributions, birth contrasts, eligible layers and move probabilities stay
% unchanged. Global moves (zone=0) have multiplier one in both directions.
log_tf=0; log_tr=0;
if strcmp(mode,'targeted_layered')
    assert(isequal(p2.previous_zone,ptb.zone), ...
        'thbi:dr','Second-stage targeting does not match the first move.');
    w=par.proposal.dr.same_zone_weight; base=par.proposal.zone_probability;
    assert(isvector(w) && any(numel(w)==[1 3]) && isreal(w) && ...
        all(isfinite(w) & w>=0 & w<1),'thbi:dr','Invalid targeting weight.');
    if isscalar(w), w=repmat(w,1,3); end
    if ptb.zone>0 && p2.zone>0
        forward_zone=(1-w(ptb.zone))*base;
        forward_zone(ptb.zone)=forward_zone(ptb.zone)+w(ptb.zone);
        reverse_zone=(1-w(p2.zone))*base;
        reverse_zone(p2.zone)=reverse_zone(p2.zone)+w(p2.zone);
        log_tf=log(forward_zone(p2.zone)/base(p2.zone));
        log_tr=log(reverse_zone(ptb.zone)/base(ptb.zone));
    end
else
    assert(isnan(p2.previous_zone),'thbi:dr','Random path mode cannot use targeted draws.');
end
assert(abs(p2.log_targeting_forward-log_tf)<1e-12, ...
    'thbi:dr','Recorded second-stage proposal density disagrees with targeting.');
terms.log_targeting_ratio=log_tr-log_tf;
terms.log_q1_xy=ptb.log_q_forward; terms.log_q1_zy=p2.log_q_reverse;
terms.log_q2_xyz=p2.log_q_forward+log_tf;
terms.log_q2_zyx=ptb.log_q_reverse+log_tr;
terms.log_abs_jacobian=ptb.log_abs_jacobian+p2.log_abs_jacobian;

%% =================== BOTH REJECTION FACTORS AND FULL PATH ACCEPTANCE ===================
if terms.logalpha_xy>=0 || terms.logalpha_zy>=0, return; end
terms.log_reject_xy=log(-expm1(terms.logalpha_xy));
terms.log_reject_zy=log(-expm1(terms.logalpha_zy));
terms.log_target_ratio=lz-lx;
terms.log_proposal_jacobian=terms.log_q1_zy-terms.log_q1_xy + ...
    terms.log_q2_zyx-terms.log_q2_xyz+terms.log_abs_jacobian;
logalpha=terms.log_target_ratio+terms.log_proposal_jacobian + ...
    terms.log_reject_zy-terms.log_reject_xy;
assert(isreal(logalpha) && ~isnan(logalpha),'thbi:dr','Undefined targeted RJ acceptance ratio.');
accept=log(rand)<min(0,logalpha);
end
