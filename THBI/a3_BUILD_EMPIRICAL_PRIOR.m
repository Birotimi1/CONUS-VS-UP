function prior = a3_BUILD_EMPIRICAL_PRIOR(par, nsamples)
% THBI-style empirical prior: independent draws from the specified layered
% joint prior, including physical conditions; no likelihood or forward solver.
assert(isscalar(nsamples) && nsamples>0 && nsamples==fix(nsamples),'Invalid prior sample count.');
validate_configuration(par);
saved_rng=rng; restore_rng=onCleanup(@() rng(saved_rng)); %#ok<NASGU>
rng(par.inv.seed-1,'twister');
drawpar=par; drawpar.mod.start.method='prior'; drawpar.mod.start.jitter=false;

%% =================== INDEPENDENT MODELS AND COMMON-DEPTH PROFILES ===================
K=zeros(nsamples,3); profiles=zeros(nsamples,numel(par.res.depth_km));
boundaries=zeros(nsamples,2);
for i=1:nsamples
    % b1's prior mode draws the base normalized priors and rejects models
    % violating a1, exactly matching the cross-zone-conditioned joint target.
    model=b1_INITIATE_MODEL(drawpar);
    K(i,:)=model.K; boundaries(i,:)=[model.zsed model.zmoh];
    depth=cumsum(model.matrix(1:end-1,1));
    layer=1+sum(bsxfun(@ge,par.res.depth_km(:),depth'),2);
    profiles(i,:)=model.matrix(layer,3)';
end

%% =================== INDUCED COUNT MASSES AND VELOCITY ENVELOPES ===================
prior.nsamples=nsamples; prior.seed=par.inv.seed-1;
prior.depth_km=par.res.depth_km(:);
prior.vs_quantiles=sample_quantiles(profiles,[.05 .5 .95]);
prior.boundary_quantiles=sample_quantiles(boundaries,[.05 .5 .95]);
prior.count_values=cell(1,3); prior.count_probability=cell(1,3);
for z=1:3
    support=par.mod.kmin(z):par.mod.kmax(z);
    prior.count_values{z}=support;
    prior.count_probability{z}=arrayfun(@(k) mean(K(:,z)==k),support);
end
prior.note='Independent samples from the joint prior after cross-zone conditions, not raw 1/K weights.';
end
