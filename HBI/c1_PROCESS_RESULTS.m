function summary = c1_PROCESS_RESULTS(chains, data, par)
% Parameter/depth summaries from ALL retained states, not just accepted ones.

%% Collect the same retained iteration range from every chain
keep=(par.inv.burnin+1):par.inv.niter;
samples=[]; cube=[];
best=chains{1}.bestfit;
for k=1:numel(chains)
    c=chains{k}; assert(c.last_iter==par.inv.niter,'Incomplete chain.');
    x=[c.theta(keep,:),exp(c.eta(keep,:))];
    samples=[samples;x]; %#ok<AGROW>
    cube(:,:,k)=x; %#ok<AGROW>
    if sum(c.bestfit.misfit.E2)<sum(best.misfit.E2), best=c.bestfit; end
end
%% Parameter quantiles and between/within-chain mixing
summary.parameter_names=[par.mod.names,cellfun(@(s)['gamma_' s],par.inv.datatypes,'UniformOutput',false)];
summary.parameter_quantiles=sample_quantiles(samples,[.05 .5 .95]);
summary.split_rhat=split_rhat(cube);
summary.rhat_kind='Classical split R-hat (not rank-normalized); necessary diagnostic, not a convergence certificate';
summary.retained_per_chain=numel(keep);
summary.bestfit=best;
summary.bestfit_selection='Minimum sum of residual energies whitened by C0, among visited states including burn-in';
%% Map each sampled layer model to the common depth grid
summary.depth_km=par.res.depth_km;
profiles=zeros(size(samples,1),numel(par.res.depth_km));
for j=1:size(samples,1)
    z=cumsum(samples(j,par.mod.idx.h));
    ii=1+sum(bsxfun(@ge,par.res.depth_km(:),z),2);
    vs=samples(j,par.mod.idx.vs); profiles(j,:)=vs(ii);
end
summary.vs_quantiles=sample_quantiles(profiles,[.05 .5 .95]);
summary.warning='Exploratory finite-chain summaries. Check mixing, multiple starts, effective sample size and longer runs before interpreting uncertainty.';
summary.base_reduced_chi2=zeros(5,1);
for k=1:5, summary.base_reduced_chi2(k)=best.misfit.E2(k)/data.(par.inv.datatypes{k}).N; end
end
