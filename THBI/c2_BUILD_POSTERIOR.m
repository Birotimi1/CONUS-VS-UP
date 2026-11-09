function [posterior, retained] = c2_BUILD_POSTERIOR(chains, par)
% THBI-style posterior parameter collection for variable physical layers.
% Keep all post-burn-in states. Never average layer slots across different K.

%% =================== COLLATE STATES WITH CHAIN/ITERATION IDENTITIES ===================
retained=collate_allmodels_perchain(chains,par);
n=retained.Nstored;
posterior.Nstored=n;
posterior.chain=retained.chain; posterior.iteration=retained.iteration;
posterior.K=retained.K; posterior.zsed=retained.boundaries(:,1);
posterior.zmoh=retained.boundaries(:,2); posterior.halfspaceVS=retained.halfspace;
posterior.gamma=exp(retained.eta); posterior.logL=retained.logL;
posterior.base_energy=model_fit_score(retained.E2,par);
if isfield(retained,'vpvs'), posterior.vpvs=retained.vpvs; end
if isfield(retained,'xi'), posterior.xi=retained.xi; end
posterior.Nfinite=zeros(n,1); posterior.Nincluding_halfspace=zeros(n,1);
posterior.zone_mean_vs=zeros(n,3); posterior.zone_top_vs=zeros(n,3);
posterior.zone_bottom_vs=zeros(n,3); posterior.fdVSsed=zeros(n,1);
posterior.fdVSmoh=zeros(n,1);

%% =================== PHYSICAL COUNTS, VELOCITIES AND DISCONTINUITIES ===================
for i=1:n
    model=model_from_chain(retained,i,par);
    [posterior.Nfinite(i),K,posterior.Nincluding_halfspace(i)]=count_N_layers_in_mod(model,par);
    assert(isequal(K,retained.K(i,:)),'thbi:results','Stored and reconstructed counts disagree.');
    for z=1:3
        zone=model.(par.mod.zone_fields{z});
        posterior.zone_top_vs(i,z)=zone.VS(1);
        posterior.zone_bottom_vs(i,z)=zone.VS(end);
        posterior.zone_mean_vs(i,z)=sum(diff([0;zone.f(:);1]).*zone.VS(:));
    end
    pair=[model.sedmparm.VS(end),model.crustmparm.VS(1)];
    posterior.fdVSsed(i)=100*diff(pair)/mean(pair);
    posterior.fdVSmoh(i)=model.fdVSmoh;
end

%% =================== DISCRETE COUNT MASSES ===================
posterior.count_values=cell(1,3); posterior.count_probability=cell(1,3);
for z=1:3
    support=par.mod.kmin(z):par.mod.kmax(z);
    posterior.count_values{z}=support;
    posterior.count_probability{z}=arrayfun(@(k) mean(posterior.K(:,z)==k),support);
end
posterior.total_count_values=sum(par.mod.kmin):sum(par.mod.kmax);
posterior.total_count_probability=arrayfun(@(k) mean(posterior.Nfinite==k),posterior.total_count_values);
posterior.note=retained.note;
end
