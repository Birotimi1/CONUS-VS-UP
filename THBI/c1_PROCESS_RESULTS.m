function summary = c1_PROCESS_RESULTS(chains, data, par)
% Summarize RJ states on a common depth grid and as zone-count probabilities.
% Layer j has no common meaning across different K; never average those slots.

%% =================== RETAINED STATES AND COMMON METRICS ===================
keep=(par.inv.burnin+1):par.inv.niter;
n=numel(keep); nc=numel(chains); nz=numel(par.res.depth_km);
profiles=zeros(n*nc,nz); counts=zeros(n*nc,3);
nd=numel(par.inv.datatypes); nmetrics=6+nd+isfield(par.mod,'vpvs_bounds');
metrics=zeros(n,nmetrics,nc);
retained_best=all(cellfun(@(q) isfield(q,'bestfit_retained') && ...
    ~isempty(q.bestfit_retained),chains));
best_field='bestfit'; summary.bestfit_label='Best visited fit (legacy)';
summary.bestfit_selection='Legacy minimum base-covariance energy among visited states, possibly including burn-in. Cached retained-best responses were not stored in this format.';
if retained_best
    best_field='bestfit_retained'; summary.bestfit_label='Best retained fit';
    summary.bestfit_selection='Minimum fixed-reference-noise fit score among ALL retained states, including repeats; same selection as c4, not a posterior mean or cross-dimension MAP.';
end
best=chains{1}.(best_field); summary.bestfit_chain=1;
summary.bestfit_iteration=NaN;
interface_counts=zeros(1,numel(par.res.interface_bins_km)-1);
summary.move_attempts=zeros(nc,numel(par.proposal.names)); summary.move_accepts=zeros(nc,numel(par.proposal.names));
summary.move_valid=zeros(nc,numel(par.proposal.names));
summary.zone_move_attempts=zeros(nc,numel(par.proposal.names),3); summary.zone_move_accepts=zeros(nc,numel(par.proposal.names),3);
summary.retained_forward_failures=zeros(nc,1);
summary.dr_attempts=zeros(nc,1); summary.dr_accepts=zeros(nc,1);
summary.dr_move_attempts=zeros(nc,numel(par.proposal.names)); summary.dr_move_valid=zeros(nc,numel(par.proposal.names));
summary.dr_move_accepts=zeros(nc,numel(par.proposal.names));
summary.dr_zone_move_attempts=zeros(nc,numel(par.proposal.names),3); summary.dr_zone_move_accepts=zeros(nc,numel(par.proposal.names),3);
summary.move_note='First-stage acceptance is separate from accepted second-stage paths.';
for c=1:nc
    chain=chains{c}; assert(chain.last_iter==par.inv.niter,'thbi:results','Incomplete chain.');
    metrics(:,1:6+nd,c)=[chain.boundaries(keep,:),chain.halfspace(keep),exp(chain.eta(keep,:)),chain.K(keep,:)];
    if isfield(chain,'vpvs'), metrics(:,end,c)=chain.vpvs(keep); end
    counts((c-1)*n+(1:n),:)=chain.K(keep,:);
    candidate=chain.(best_field);
    if retained_best
        [energy,first_minimum]=min(model_fit_score(chain.E2(keep,:),par));
        assert(chain.bestfit_retained_iteration==keep(first_minimum) && ...
            abs(model_fit_score(candidate.misfit.E2',par)-energy)<=1e-10*(1+abs(energy)), ...
            'thbi:results','Retained-best cache disagrees with the stored chain.');
    end
    if model_fit_score(candidate.misfit.E2',par)<model_fit_score(best.misfit.E2',par)
        best=candidate; summary.bestfit_chain=c;
    end
    first_accepted=chain.accepted; first_valid=chain.proposal_valid;
    if isfield(chain,'accepted_stage'), first_accepted=chain.accepted_stage==1; end
    if isfield(chain,'first_proposal_valid'), first_valid=chain.first_proposal_valid; end
    for m=1:numel(par.proposal.names)
        summary.move_attempts(c,m)=sum(chain.move(keep)==m);
        summary.move_accepts(c,m)=sum(first_accepted(keep) & chain.move(keep)==m);
        summary.move_valid(c,m)=sum(first_valid(keep) & chain.move(keep)==m);
        for z=1:3
            selected=chain.move(keep)==m & chain.zone(keep)==z;
            summary.zone_move_attempts(c,m,z)=sum(selected);
            summary.zone_move_accepts(c,m,z)=sum(first_accepted(keep) & selected);
        end
    end
    summary.retained_forward_failures(c)=sum(chain.forward_failed(keep));
    if isfield(chain,'dr_attempted')
        summary.dr_attempts(c)=sum(chain.dr_attempted(keep));
        summary.dr_accepts(c)=sum(chain.accepted_stage(keep)==2);
        % New files preserve the independently selected second family/zone.
        % Old scalar-DR files used the same first/second family and zone.
        dm=chain.move; dz=chain.zone;
        if isfield(chain,'dr_move'), dm=chain.dr_move; dz=chain.dr_zone; end
        for m=1:numel(par.proposal.names)
            selected=chain.dr_attempted(keep) & dm(keep)==m;
            summary.dr_move_attempts(c,m)=sum(selected);
            summary.dr_move_valid(c,m)=sum(selected & chain.dr_proposal_valid(keep));
            summary.dr_move_accepts(c,m)=sum(selected & chain.accepted_stage(keep)==2);
            for z=1:3
                selected_zone=selected & dz(keep)==z;
                summary.dr_zone_move_attempts(c,m,z)=sum(selected_zone);
                summary.dr_zone_move_accepts(c,m,z)=sum(selected_zone & chain.accepted_stage(keep)==2);
            end
        end
    end

    %% =================== DEPTH PROFILE AND INTERFACES FOR EACH DRAW ===================
    for j=1:n
        model=model_from_chain(chain,keep(j),par);
        depth=cumsum(model.matrix(1:end-1,1));
        layer=1+sum(bsxfun(@ge,par.res.depth_km(:),depth'),2);
        profiles((c-1)*n+j,:)=model.matrix(layer,3)';
        % The fixed model base is not a sampled physical interface.
        interface_counts=interface_counts+histcounts(depth(1:end-1),par.res.interface_bins_km);
    end
end

%% =================== POSTERIOR VELOCITIES AND LAYER-COUNT MASSES ===================
summary.depth_km=par.res.depth_km(:);
summary.vs_quantiles=sample_quantiles(profiles,[.05 .5 .95]);
% Evaluate mixing at fixed depths as well as for K and boundary scalars.
% Work in blocks to avoid duplicating the entire production profile array.
summary.vs_split_rhat=NaN(1,nz);
summary.vs_rank_rhat=NaN(1,nz); summary.vs_bulk_ess=NaN(1,nz); summary.vs_tail_ess=NaN(1,nz);
for first=1:32:nz
    columns=first:min(first+31,nz);
    x=permute(reshape(profiles(:,columns),n,nc,numel(columns)),[1 3 2]);
    summary.vs_split_rhat(columns)=split_rhat(x);
    diagnostic=mcmc_rank_diagnostics(x);
    summary.vs_rank_rhat(columns)=diagnostic.rhat;
    summary.vs_bulk_ess(columns)=diagnostic.bulk_ess;
    summary.vs_tail_ess(columns)=diagnostic.tail_ess;
end
summary.count_values=cell(1,3); summary.count_probability=cell(1,3);
for z=1:3
    support=par.mod.kmin(z):par.mod.kmax(z);
    summary.count_values{z}=support;
    summary.count_probability{z}=arrayfun(@(k) mean(counts(:,z)==k),support);
end
summary.interface_bins_km=par.res.interface_bins_km;
summary.expected_interfaces_per_bin=interface_counts/(n*nc);
summary.retained_per_chain=n;
summary.prior=[];
if exist(fullfile(par.res.dir,'prior.mat'),'file')
    saved=load(fullfile(par.res.dir,'prior.mat'),'prior');
    summary.prior=saved.prior;
end
summary.bestfit=best;
if retained_best
    summary.bestfit_iteration=chains{summary.bestfit_chain}.bestfit_retained_iteration;
end

%% =================== SAMPLED RESPONSE DISTRIBUTIONS ===================
summary.response_quantiles=struct(); summary.nresponse_per_chain=0;
summary.response_note='Pointwise 5/50/95 percentiles of saved model responses; observational noise is not added.';
if ~par.inv.prior_only && par.res.prediction_stride>0
    iterations=chains{1}.prediction_iterations;
    summary.nresponse_per_chain=numel(iterations);
    for c=1:nc
        assert(isequal(chains{c}.prediction_iterations,iterations) && all(chains{c}.prediction_saved), ...
            'thbi:results','Incomplete or inconsistent saved response history.');
    end
    for d=1:numel(par.inv.datatypes)
        name=par.inv.datatypes{d};
        responses=cellfun(@(chain) chain.predictions.(name),chains,'UniformOutput',false);
        responses=vertcat(responses{:});
        assert(size(responses,2)==data.(name).N && all(isfinite(responses(:))), ...
            'thbi:results','Invalid saved predictions for %s.',name);
        summary.response_quantiles.(name)=sample_quantiles(responses,[.05 .5 .95]);
    end
end

%% =================== FIXED-METRIC DIAGNOSTICS ===================
summary.metric_names=[{'zsed','zmoh','Vs_halfspace'}, ...
    cellfun(@(d) ['gamma_' d],par.inv.datatypes,'UniformOutput',false), ...
    {'K_sed','K_crust','K_mantle'}];
if isfield(par.mod,'vpvs_bounds')
    summary.metric_names{end+1}='crust_vpvs';
    summary.metric_names{8}='sigma_Ps_HK'; summary.metric_names{9}='sigma_Sp_CCP';
end
summary.metric_quantiles=sample_quantiles(reshape(permute(metrics,[1 3 2]),n*nc,nmetrics),[.05 .5 .95]);
summary.split_rhat=split_rhat(metrics);
summary.rank_diagnostics=mcmc_rank_diagnostics(metrics);
summary.rhat_kind='Classical split and rank-normalized/folded R-hat; bulk and tail ESS. Diagnostics are estimates, not recovery guarantees.';
summary.warning='Pilot summaries do not establish mixing, count recovery, identifiability or uncertainty coverage.';
summary.base_reduced_chi2=NaN(nd,1);
if ~par.inv.prior_only
    for d=1:nd
        v=data.(par.inv.datatypes{d});
        dof=v.N; if isfield(v,'dof'), dof=v.dof; end
        summary.base_reduced_chi2(d)=best.misfit.E2(d)/dof;
    end
end
end
