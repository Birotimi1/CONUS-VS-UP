function test_thbi_chain(profile,dr_mode,include_summaries)
% Cheap analytic likelihood checks, cache checks and exact checkpoint resume.
if nargin<1, profile='pilot'; end
if nargin<3, include_summaries=true; end
% Octave kernel runs may omit MATLAB-only histcounts summary processing.
par=bayes_inv_parms(profile); par.mod.start.jitter=false;
if nargin>1
    par.proposal.dr.mode=dr_mode;
    % Explicit optional comparison; the production profile uses random moves.
    if strcmp(dr_mode,'targeted_layered'), par.proposal.dr.same_zone_weight=.5; end
end
nd=numel(par.inv.datatypes);
if par.proposal.dr.enabled, par.proposal.noise_step=.7; end
par.inv.niter=240; par.inv.burnin=80; par.inv.print_every=1000; par.inv.checkpoint_every=37;
root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
model=b1_INITIATE_MODEL(par);
for d=1:nd
    name=par.inv.datatypes{d};
    data.(name)=struct('x',(1:3)','obs',[1;2;3],'sig',[.5;.8;1.1], ...
        'R',[1 .2 .1;.2 1 .3;.1 .3 1]);
end
if nd==6, data.rf_sp.x=[-2;0;2]; end
data=a2_PREPARE_DATA(data,par); predicted=toy_forward(model,data,par,'');
misfit=b4_CALC_MISFIT(data,predicted,par);
[L,~]=b5_CALC_LIKELIHOOD(misfit,data,model.datahparm,par); expected=0;
for d=1:nd
    name=par.inv.datatypes{d}; r=data.(name).obs-predicted.(name);
    C=diag(data.(name).sig)*data.(name).R*diag(data.(name).sig);
    expected=expected-.5*(3*log(2*pi)+log(det(C))+r'*(C\r));
end
assert(abs(L-expected)<1e-10,'Correlated Gaussian likelihood mismatch.');
noise=model.datahparm; noise.eta(5)=log(2);
[L2,~]=b5_CALC_LIKELIHOOD(misfit,data,noise,par);
assert(abs(L2-L-(-3*log(2)+.375*misfit.E2(5)))<1e-10,'Noise normalization mismatch.');

%% =================== NOISE MOVES DO NOT CALL FORWARD SOLVERS ===================
original_p=par.proposal.probability;
par.proposal.probability=[0 0 0 0 0 1]; par.res.dir=fullfile(root,'noise');
c=run_one_chain(par,data,1,[],@toy_forward);
assert(c.nforward==1 && all(all(c.K==repmat(model.K,par.inv.niter,1))));
assert(size(c.K,1)==par.inv.niter && any(~c.accepted),'Rejected states were not retained.');

%% =================== UNINTERRUPTED VERSUS RESUMED RJ CHAIN ===================
par.proposal.probability=original_p; par.res.dir=fullfile(root,'continuous');
a=run_one_chain(par,data,1,[],@toy_forward);
par.res.dir=fullfile(root,'resumed');
run_one_chain(par,data,1,91,@toy_forward);
b=run_one_chain(par,data,1,[],@toy_forward);
for field={'fractions','velocities','K','boundaries','halfspace','eta','logL','logP','E2', ...
        'bestfit_retained','bestfit_retained_iteration','bestfit_iteration', ...
        'accepted_stage','dr_attempted','dr_proposal_valid','dr_logalpha','first_logalpha', ...
        'dr_move','dr_zone','dr_K_before','dr_K_after','first_K_after', ...
        'dr_log_q_ratio','dr_log_jacobian','dr_log_targeting_ratio','dr_reverse_first_logalpha','first_forward_called','dr_forward_called', ...
        'move','accepted','logalpha','predictions','prediction_iterations','prediction_saved','forward_failed'}
    assert(isequaln(a.(field{1}),b.(field{1})),'Resumed history mismatch: %s',field{1});
end
assert(b.nforward==1+sum(b.first_forward_called)+sum(b.dr_forward_called), ...
    'Forward-cache accounting mismatch.');
assert(isequal(b.first_forward_called,b.first_proposal_valid & b.move~=6) && ...
    isequal(b.dr_forward_called,b.dr_proposal_valid & b.dr_move~=6));
assert(isequal(b.accepted,b.accepted_stage>0),'Invalid DR stage accounting.');
if strcmp(par.proposal.dr.mode,'fixed_scalar')
    assert(~any(b.dr_attempted & ismember(b.move,[4 5])));
else
    rows=b.dr_attempted;
    assert(isequal(b.dr_K_before(rows,:),b.first_K_after(rows,:)));
end
if par.proposal.dr.enabled
    assert(any(b.dr_attempted) && any(b.accepted_stage==2),'DR was not exercised.');
else
    assert(~any(b.dr_attempted),'Disabled DR was attempted.');
end
for i=1:par.inv.niter
    state=model_from_chain(b,i,par);
    assert(isequal(state.K,b.K(i,:)) && a1_TEST_CONDITIONS(state,par));
    % Check EVERY state's likelihood, including first-structure/second-noise
    % paths: second noise must reuse rejected y's response, not accepted x's.
    prediction=toy_forward(state,data,par,'');
    check=b4_CALC_MISFIT(data,prediction,par);
    [logL,~]=b5_CALC_LIKELIHOOD(check,data,state.datahparm,par);
    assert(max(abs(check.E2'-b.E2(i,:)))<1e-10 && abs(logL-b.logL(i))<1e-10);
    priorK=b.initial.model.K; if i>1, priorK=b.K(i-1,:); end
    expected=zeros(1,3);
    use_first=b.accepted_stage(i)==1 || ...
        (b.accepted_stage(i)==2 && any(strcmp(par.proposal.dr.mode,{'layered_path','targeted_layered'})));
    if use_first && any(b.move(i)==[4 5])
        expected(b.zone(i))=double(b.move(i)==4)-double(b.move(i)==5);
    end
    if b.accepted_stage(i)==2 && any(b.dr_move(i)==[4 5])
        z=b.dr_zone(i); expected(z)=expected(z)+double(b.dr_move(i)==4)-double(b.dr_move(i)==5);
    end
    assert(isequal(state.K-priorK,expected),'Accepted path layer counts disagree.');
end
c=run_one_chain(par,data,2,[],@toy_forward);
assert(~isequaln(b.eta,c.eta),'Independent chain seeds produced identical histories.');

%% =================== SAVED RESPONSES MATCH THEIR EXACT SAVED EARTH ===================
for row=1:numel(b.prediction_iterations)
    i=b.prediction_iterations(row);
    predicted=toy_forward(model_from_chain(b,i,par),data,par,'');
    for d=1:nd
        name=par.inv.datatypes{d};
        assert(max(abs(b.predictions.(name)(row,:)'-predicted.(name)))<1e-12);
    end
end
[~,j]=min(sum(b.E2(par.inv.burnin+1:end,:),2));
assert(b.bestfit_retained_iteration==par.inv.burnin+j);
retained_prediction=toy_forward(b.bestfit_retained.model,data,par,'');
assert(isequal(retained_prediction,b.bestfit_retained.predata), ...
    'thbi:test','Retained-best cache does not match its physical model.');
if include_summaries
    summary=c1_PROCESS_RESULTS({b;c},data,par);
    keep=(par.inv.burnin+1):par.inv.niter;
    assert(sum(summary.move_accepts(1,:))==sum(b.accepted_stage(keep)==1) && ...
        sum(summary.dr_move_accepts(1,:))==sum(b.accepted_stage(keep)==2));
    assert(sum(summary.dr_move_attempts(1,:))==sum(b.dr_attempted(keep)));
    assert(isequal(size(summary.vs_quantiles),[3 numel(par.res.depth_km)]));
    assert(numel(summary.vs_split_rhat)==numel(par.res.depth_km));
    assert(summary.nresponse_per_chain==numel(b.prediction_iterations));
    for d=1:nd
        name=par.inv.datatypes{d};
        q=summary.response_quantiles.(name);
        assert(isequal(size(q),[3 data.(name).N]) && all(all(diff(q,1,1)>=0)));
    end
end
% A deterministic trending sequence has this analytic classical split R-hat.
assert(abs(split_rhat(repmat((1:8)',1,1,2))-sqrt(79/20))<1e-12);
assert(isnan(split_rhat(ones(8,1,2))),'Identical constant chains have undefined R-hat.');
assert(isinf(split_rhat(cat(3,ones(8,1),2*ones(8,1)))),'Distinct constant chains need infinite R-hat.');

%% =================== INVALID INPUTS AND INCOMPATIBLE RESUME ===================
bad=data; bad.rf.obs(1)=bad.rf.obs(1)+1i;
expect_error(@() a2_PREPARE_DATA(bad,par),'hbi:data');
bad=data; bad.rf.R(1,2)=.2+.1i; bad.rf.R(2,1)=.2-.1i;
expect_error(@() a2_PREPARE_DATA(bad,par),'hbi:covariance');
predicted=toy_forward(model,data,par,''); predicted.rf(1)=predicted.rf(1)+1i;
expect_error(@() b4_CALC_MISFIT(data,predicted,par),'hbi:prediction');
changed=data; changed.rf.obs(1)=changed.rf.obs(1)+.1;
try
    run_one_chain(par,changed,1,[],@toy_forward);
    error('thbi:test','Changed data were incorrectly accepted for resume.');
catch problem
    assert(strcmp(problem.identifier,'thbi:resume'),'Unexpected resume error.');
end
expect_error(@() run_one_chain(par,data,1,1.5,@toy_forward),'thbi:chain');
expect_error(@() run_one_chain(par,data,0,[],@toy_forward),'thbi:chain');
bad=par; bad.proposal.dr.same_zone_weight=1; bad.proposal.dr.mode='targeted_layered';
expect_error(@() validate_configuration(bad),'thbi:config');
bad=par; bad.proposal.dr.mode='unknown'; expect_error(@() validate_configuration(bad),'thbi:config');
bad=par; bad.res.nprior=-1; expect_error(@() validate_configuration(bad),'thbi:config');
bad=par; bad.synth.seed=-1; expect_error(@() validate_configuration(bad),'thbi:config');
bad=misfit; bad.E2(1)=Inf;
expect_error(@() b5_CALC_LIKELIHOOD(bad,data,model.datahparm,par),'hbi:likelihood');
bad=toy_forward(model,data,par,''); bad.rf(:)=1e200;
expect_error(@() b4_CALC_MISFIT(data,bad,par),'hbi:misfit');

%% =================== DISABLED OPTIONAL OUTPUTS KEEP VALID RETAINED FITS ===================
disabled=par; disabled.inv.niter=40; disabled.inv.burnin=10;
disabled.res.nprior=0; disabled.res.prediction_stride=0;
disabled.res.dir=fullfile(root,'no_response_history');
minimal=run_one_chain(disabled,data,1,[],@toy_forward);
assert(isempty(minimal.prediction_iterations) && isempty(fieldnames(minimal.predictions)));
assert(minimal.bestfit_retained_iteration>disabled.inv.burnin);
if include_summaries
    summary=c1_PROCESS_RESULTS({minimal},data,disabled);
    assert(summary.nresponse_per_chain==0 && isempty(fieldnames(summary.response_quantiles)) && ...
        summary.bestfit_iteration>disabled.inv.burnin);
end
fprintf('PASS: Gaussian target, noise/cache, repeated states, saved responses, exact resume, seeds and input checks.\n');
end

function expect_error(callback, identifier)
try
    callback();
    error('thbi:test','Expected error %s was not raised.',identifier);
catch problem
    assert(strcmp(problem.identifier,identifier),'Unexpected error: %s',problem.identifier);
end
end

function predicted=toy_forward(model,data,par,~)
% Deliberately inexpensive physical-model functional; no external executable.
v=sum(model.matrix(:,1).*model.matrix(:,3))/par.mod.maxz+.1*model.halfspaceVS;
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; predicted.(name)=v*(1:data.(name).N)'/10;
end
end
