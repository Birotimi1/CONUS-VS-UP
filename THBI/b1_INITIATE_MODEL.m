function model = b1_INITIATE_MODEL(par)
% Initialize zoned layer counts, boundaries, velocities and data errors.
% The actual draw/jitter blocks are here; assembly is in functions/.

%% =================== REFERENCE START OR DISPERSED PRIOR DRAW ===================
for attempt = 1:par.mod.start.max_attempts
    model = struct();
    if strcmp(par.mod.start.method,'reference')
        model.zsed = par.mod.start.zsed;
        model.zmoh = par.mod.start.zmoh;
        model.halfspaceVS = par.mod.start.halfspace;
    elseif strcmp(par.mod.start.method,'prior')
        b = par.mod.boundary_lower+rand(1,2).*(par.mod.boundary_upper-par.mod.boundary_lower);
        model.zsed=b(1); model.zmoh=b(2);
        model.halfspaceVS=par.mod.halfspace_bounds(1)+rand*diff(par.mod.halfspace_bounds);
    else
        error('thbi:start','Starting method must be reference or prior.');
    end

    %% =================== SEDIMENT, CRUST AND MANTLE PARAMETERS ===================
    for z = 1:3
        field = par.mod.zone_fields{z};
        if strcmp(par.mod.start.method,'reference')
            f=par.mod.start.f{z}(:); v=par.mod.start.vs{z}(:);
        else
            counts=par.mod.kmin(z):par.mod.kmax(z);
            if strcmp(par.mod.count_prior,'jeffreys'), p=1./counts; else, p=ones(size(counts)); end
            cumulative=cumsum(p/sum(p)); cumulative(end)=1;
            K=counts(find(rand<cumulative,1,'first'));
            spacing=par.mod.min_fraction(z); S=1-K*spacing;
            f=sort(rand(K-1,1))*S+(1:K-1)'*spacing;
            v=par.mod.vsmin(z)+rand(K,1)*(par.mod.vsmax(z)-par.mod.vsmin(z));
            if par.conditions.monotonic(z), v=sort(v); end
        end
        model.(field)=struct('f',f,'VS',v);
    end

    %% =================== REFERENCE-START JITTER ===================
    if strcmp(par.mod.start.method,'reference') && par.mod.start.jitter
        b=[model.zsed model.zmoh];
        for j=1:2
            b(j)=reflect_to_interval(b(j)+par.proposal.boundary_step(j)*randn, ...
                par.mod.boundary_lower(j),par.mod.boundary_upper(j));
        end
        model.zsed=b(1); model.zmoh=b(2);
        for z=1:3
            field=par.mod.zone_fields{z}; v=model.(field).VS; f=model.(field).f;
            for j=randperm(numel(v))
                lo=par.mod.vsmin(z); hi=par.mod.vsmax(z);
                if par.conditions.monotonic(z)
                    if j>1, lo=max(lo,v(j-1)); end
                    if j<numel(v), hi=min(hi,v(j+1)); end
                end
                v(j)=reflect_to_interval(v(j)+par.proposal.vs_step(z)*randn,lo,hi);
            end
            for j=randperm(numel(f))
                e=[0;f;1]; gap=par.mod.min_fraction(z);
                f(j)=reflect_to_interval(f(j)+par.proposal.interface_step(z)*randn,e(j)+gap,e(j+2)-gap);
            end
            model.(field).VS=v; model.(field).f=f;
        end
    end

    %% =================== CRUSTAL VP/VS AND HIERARCHICAL DATA PARAMETERS ===================
    if isfield(par.mod,'vpvs_bounds')
        model.vpvs=par.mod.start.vpvs;
        if strcmp(par.mod.start.method,'prior')
            model.vpvs=par.mod.vpvs_bounds(1)+rand*diff(par.mod.vpvs_bounds);
        end
    end
    if isfield(par.mod,'xi_bounds')
        model.xi=par.mod.start.xi;
        if strcmp(par.mod.start.method,'prior')
            model.xi=(par.mod.xi_bounds(:,1)+rand(3,1).*diff(par.mod.xi_bounds,1,2))';
        end
    end
    model.datahparm.eta=par.mod.noise.log_center;
    if strcmp(par.mod.start.method,'prior')
        for d=1:numel(par.inv.datatypes)
            if isfield(par.mod.noise,'log_uniform') && par.mod.noise.log_uniform(d)
                model.datahparm.eta(d)=par.mod.noise.lower(d)+rand*(par.mod.noise.upper(d)-par.mod.noise.lower(d));
                continue;
            end
            eta=Inf;
            while eta<par.mod.noise.lower(d) || eta>par.mod.noise.upper(d)
                eta=par.mod.noise.log_center(d)+par.mod.noise.log_width(d)*randn;
            end
            model.datahparm.eta(d)=eta;
        end
    end

    %% =================== CONDITIONS AND PHYSICAL MODEL ===================
    [ok,why]=a1_TEST_CONDITIONS(model,par);
    if ok
        model=make_mod_from_parms(model,par);
        return;
    end
end
error('thbi:start','No valid start after %d attempts: %s',attempt,why);
end
