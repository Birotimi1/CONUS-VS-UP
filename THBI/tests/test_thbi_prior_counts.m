function report=test_thbi_prior_counts(profile,proposal_mode,niter,seed,driver)
% Optional stochastic prior-only check. No CPS/HVf calls. This diagnostic
% is not a proof of posterior convergence and can be repeated across seeds.
if nargin<1, profile='pilot'; end
if nargin<2, proposal_mode='production_steps'; end
if nargin<3, niter=30000; end
if nargin<4, seed=42000; end
if nargin<5, driver='chain'; end
par=bayes_inv_parms(profile);
par.inv.prior_only=true; par.inv.niter=niter; par.inv.burnin=floor(niter/6); par.inv.seed=seed;
par.inv.print_every=5000; par.inv.checkpoint_every=5000;
par.mod.kmax=[4 5 4]; par.mod.start.method='prior';
par.conditions.pos_moho=false; par.conditions.pos_sed2basement=false;
par.conditions.max_moho_fraction=Inf; % removes count-dependent cross-zone conditioning
par.conditions.max_moho_jump=Inf;
par.proposal.birth_vs_sigma=[.4 .3 .3];
if strcmp(proposal_mode,'wide_prior_steps')
    % Controlled mixing diagnostic over the broad PRIOR, not production
    % tuning. Tiny posterior-oriented Vs steps can give misleading count
    % frequencies in a short prior-only random walk. Leave defaults intact.
    par.proposal.vs_step=(par.mod.vsmax-par.mod.vsmin)/3;
    par.proposal.interface_step=[.1 .1 .1];
elseif ~strcmp(proposal_mode,'production_steps')
    error('thbi:test','Choose production_steps or wide_prior_steps.');
end
fprintf('Prior-count diagnostic proposal mode: %s\n',proposal_mode);
root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
par.res.dir=root;
if strcmp(driver,'chain')
    chain=run_one_chain(par,struct(),1);
elseif strcmp(driver,'kernel')
    % Long prior diagnostics need only counts, not every Earth/response and
    % checkpoint. First require an exact match to the production driver.
    pilot=par; pilot.inv.niter=min(500,niter); pilot.inv.burnin=0;
    pilot.res.dir=fullfile(root,'driver_comparison');
    actual=run_one_chain(pilot,struct(),1);
    compact=prior_kernel(pilot);
    assert(isequal(actual.K,compact.K) && isequal(actual.accepted_stage,compact.accepted_stage) && ...
        isequal(actual.logP,compact.logP),'thbi:test','Compact diagnostic diverged from run_one_chain.');
    chain=prior_kernel(par);
else
    error('thbi:test','Choose chain or kernel diagnostic driver.');
end
assert(chain.nforward==0,'Prior-only sampler invoked a forward model.');
report=struct('proposal_mode',proposal_mode,'niter',niter,'seed',seed,'driver',driver,'maximum_count_error',zeros(1,3), ...
    'target',{{}},'observed',{{}});
for z=1:3
    support=par.mod.kmin(z):par.mod.kmax(z);
    target=(1./support)/sum(1./support);
    retained=chain.K(par.inv.burnin+1:end,z);
    observed=arrayfun(@(k) mean(retained==k),support);
    fprintf('%s: K / analytic base prior / sampled frequency\n',par.mod.zone_names{z});
    disp([support' target' observed']);
    report.target{z}=target; report.observed{z}=observed;
    report.maximum_count_error(z)=max(abs(target-observed));
end
assert(all(report.maximum_count_error<.08), ...
    'Prior-count discrepancy: inspect proposal ratios and mixing; repeat with longer runs.');
fprintf('PASS: empirical count frequencies close to analytic 1/K masses with cross-zone conditions disabled.\n');
end

function result=prior_kernel(par)
% Count-only execution of the production kernels; parity checked above.
assert(par.proposal.dr.enabled && any(strcmp(par.proposal.dr.mode,{'layered_path','targeted_layered'})));
saved=rng; restore=onCleanup(@() rng(saved)); %#ok<NASGU>
rng(par.inv.seed+1009,'twister'); validate_configuration(par);
x.model=b1_INITIATE_MODEL(par); x.logP=a3_CALC_LOG_PRIOR(x.model,par); x.logL=0;
result=struct('K',zeros(par.inv.niter,3),'logP',zeros(par.inv.niter,1), ...
    'accepted_stage',zeros(par.inv.niter,1),'nforward',0);
for i=1:par.inv.niter
    y=struct('logP',-Inf,'logL',-Inf);
    [y.model,p1]=b2_PERTURB_MODEL(x.model,par);
    if p1.valid, y.logP=a3_CALC_LOG_PRIOR(y.model,par); y.logL=0; end
    [accepted,~]=b6_IFACCEPT(x,y,p1);
    if accepted
        x=y; result.accepted_stage(i)=1;
    elseif p1.dr_eligible && isfinite(y.logP+y.logL)
        z=struct('logP',-Inf,'logL',-Inf);
        request=struct();
        if strcmp(par.proposal.dr.mode,'targeted_layered'), request.previous_zone=p1.zone; end
        [z.model,z.ptb]=b2_PERTURB_MODEL(y.model,par,request);
        if z.ptb.valid, z.logP=a3_CALC_LOG_PRIOR(z.model,par); z.logL=0; end
        [accepted,~]=delayed_rejection_layered(x,y,z,p1,par);
        if accepted, x=z; result.accepted_stage(i)=2; end
    end
    result.K(i,:)=x.model.K; result.logP(i)=x.logP;
    if mod(i,10000)==0, fprintf('Prior kernel: %d/%d\n',i,par.inv.niter); end
end
end
