function [model, ptb] = b2_PERTURB_MODEL(model, par, request)
% Propose velocity, interface, zone boundary, birth, death or log-noise move.
% All RJ transformations and densities are written in this main body.
% Optional request.type/zone/index/split_f/delta_vs/value_to are for reproducible
% debugging/tests; normal MCMC calls this function with only model and par.
if nargin<3, request=struct(); end
dr_enabled=isfield(par.proposal,'dr') && par.proposal.dr.enabled;
dr_mode='fixed_scalar';
if isfield(par.proposal,'dr') && isfield(par.proposal.dr,'mode')
    dr_mode=par.proposal.dr.mode;
end
% Stage one always uses the base zone weights. A targeted second proposal
% may pass previous_zone; all probability corrections remain explicit.
previous_zone=NaN; zone_draw=par.proposal.zone_probability;
if isfield(request,'previous_zone')
    previous_zone=request.previous_zone;
    assert(isscalar(previous_zone) && any(previous_zone==0:3), ...
        'thbi:move','Invalid first-proposal zone.');
    assert(strcmp(dr_mode,'targeted_layered'),'thbi:move','Targeting requires targeted_layered mode.');
    if previous_zone>0
        weight=par.proposal.dr.same_zone_weight;
        if numel(weight)==3, weight=weight(previous_zone); end
        zone_draw=(1-weight)*zone_draw;
        zone_draw(previous_zone)=zone_draw(previous_zone)+weight;
    end
end
step_scale=1;
if isfield(request,'step_scale'), step_scale=request.step_scale; end
assert(isscalar(step_scale) && isfinite(step_scale) && step_scale>0,'thbi:move','Invalid step scale.');
names=par.proposal.names; probabilities=par.proposal.probability;
if isfield(request,'type')
    itype=find(strcmp(request.type,names),1);
    assert(~isempty(itype),'thbi:move','Unknown move type.');
else
    cumulative=cumsum(probabilities); cumulative(end)=1;
    itype=find(rand<cumulative,1,'first');
end
ptb=struct('type',names{itype},'type_index',itype,'zone',0,'index',0, ...
    'structure',~strcmp(names{itype},'noise'),'valid',false,'reason','', ...
    'log_q_forward',0,'log_q_reverse',0,'log_abs_jacobian',0, ...
    'log_hastings',0,'split_f',NaN,'delta_vs',NaN, ...
    'K_before',model.K,'K_after',model.K,'dr_eligible',false, ...
    'value_from',NaN,'value_to',NaN,'sigma',NaN, ...
    'step_scale',step_scale,'dr_mode',dr_mode, ...
    'previous_zone',previous_zone,'log_targeting_forward',0);
if probabilities(itype)==0, ptb.reason='Move disabled by configuration.'; return; end
if any(strcmp(ptb.type,{'velocity','interface','birth','death'}))
    if isfield(request,'zone')
        z=request.zone;
    else
        cumulative=cumsum(zone_draw); cumulative(end)=1;
        z=find(rand<cumulative,1,'first');
    end
    assert(isscalar(z) && any(z==1:3),'thbi:move','Invalid zone.');
    field=par.mod.zone_fields{z}; zone=model.(field);
    v=zone.VS(:); f=zone.f(:); K=numel(v); edges=[0;f;1];
    gap=par.mod.min_fraction(z); ptb.zone=z;
end

switch ptb.type
    %% =================== WITHIN-DIMENSION VELOCITY MOVE ===================
    case 'velocity'
        halfspace=isfield(request,'halfspace') && request.halfspace;
        if ~isfield(request,'index') && ~isfield(request,'halfspace')
            halfspace=rand<par.proposal.halfspace_probability;
        end
        if halfspace
            ptb.value_from=model.halfspaceVS; ptb.sigma=step_scale*par.proposal.halfspace_step;
            value=ptb.value_from+ptb.sigma*randn;
            if isfield(request,'value_to'), value=request.value_to; end
            if ~dr_enabled, value=reflect_to_interval(value,par.mod.halfspace_bounds(1),par.mod.halfspace_bounds(2)); end
            model.halfspaceVS=value; ptb.value_to=value; ptb.dr_eligible=dr_enabled;
            ptb.zone=0;
        else
            if isfield(request,'index'), j=request.index; else, j=randi(K); end
            assert(isscalar(j) && j==fix(j) && j>=1 && j<=K,'thbi:move','Invalid velocity index.');
            lo=par.mod.vsmin(z); hi=par.mod.vsmax(z);
            if par.conditions.monotonic(z)
                if j>1, lo=max(lo,v(j-1)); end
                if j<K, hi=min(hi,v(j+1)); end
            end
            ptb.value_from=v(j); ptb.sigma=step_scale*par.proposal.vs_step(z);
            value=v(j)+ptb.sigma*randn;
            if isfield(request,'value_to'), value=request.value_to; end
            if ~dr_enabled, value=reflect_to_interval(value,lo,hi); end
            v(j)=value; ptb.value_to=value; ptb.dr_eligible=dr_enabled;
            model.(field).VS=v; ptb.index=j;
        end

    %% =================== WITHIN-ZONE INTERFACE MOVE ===================
    case 'interface'
        if K==1, ptb.reason='No internal interface in this zone.'; return; end
        if isfield(request,'index'), j=request.index; else, j=randi(K-1); end
        assert(isscalar(j) && j==fix(j) && j>=1 && j<K,'thbi:move','Invalid interface index.');
        ptb.value_from=f(j); ptb.sigma=step_scale*par.proposal.interface_step(z);
        value=f(j)+ptb.sigma*randn;
        if isfield(request,'value_to'), value=request.value_to; end
        if ~dr_enabled, value=reflect_to_interval(value,edges(j)+gap,edges(j+2)-gap); end
        f(j)=value; ptb.value_to=value; ptb.dr_eligible=dr_enabled;
        model.(field).f=f; ptb.index=j;

    %% =================== BASEMENT OR MOHO DEPTH MOVE ===================
    case 'boundary'
        b=[model.zsed model.zmoh];
        if isfield(request,'index'), j=request.index; else, j=randi(2); end
        assert(isscalar(j) && any(j==[1 2]),'thbi:move','Invalid boundary index.');
        ptb.value_from=b(j); ptb.sigma=step_scale*par.proposal.boundary_step(j);
        value=b(j)+ptb.sigma*randn;
        if isfield(request,'value_to'), value=request.value_to; end
        if ~dr_enabled, value=reflect_to_interval(value,par.mod.boundary_lower(j),par.mod.boundary_upper(j)); end
        b(j)=value; ptb.value_to=value; ptb.dr_eligible=dr_enabled;
        model.zsed=b(1); model.zmoh=b(2); ptb.index=j;
        % Fractional interfaces remain fixed. In this parameter measure the
        % boundary random walk is symmetric and no rescaling Jacobian enters.

    %% =================== BIRTH: SPLIT ONE PHYSICAL LAYER ===================
    case 'birth'
        if K>=par.mod.kmax(z), ptb.reason='Maximum layer count reached.'; return; end
        eligible=find(diff(edges)>2*gap); E=numel(eligible);
        if E==0, ptb.reason='No layer can be split with the required separation.'; return; end
        if isfield(request,'index'), j=request.index; else, j=eligible(randi(E)); end
        assert(isscalar(j) && any(j==eligible),'thbi:move','Requested layer is not splittable.');
        a=edges(j); b=edges(j+1); width=b-a-2*gap;
        if isfield(request,'split_f'), s=request.split_f;
        else, s=a+gap+rand*width; end
        if ~isscalar(s) || ~isfinite(s) || s<a+gap || s>b-gap
            ptb.reason='Split point is outside the eligible interval.'; return;
        end
        sigma=par.proposal.birth_vs_sigma(z);
        if isfield(par.proposal,'birth_broad_probability') && ...
                par.proposal.birth_broad_probability>0 && ~isfield(request,'delta_vs') && ...
                rand<par.proposal.birth_broad_probability
            sigma=par.proposal.birth_vs_sigma_broad(z);
        end
        if isfield(request,'delta_vs'), delta=request.delta_vs;
        else, delta=sigma*randn; end
        assert(isscalar(delta) && isfinite(delta),'thbi:move','Invalid split contrast.');

        % Invertible map, preserving the thickness-weighted mean velocity:
        % v_left = v-(1-r)*delta, v_right = v+r*delta, r=(s-a)/(b-a).
        % Auxiliary coordinates are s (a fraction) and delta (km/s).
        % The full map (old parameters,s,delta)->new parameters has |J|=1.
        r=(s-a)/(b-a);
        left=v(j)-(1-r)*delta; right=v(j)+r*delta;
        model.(field).f=[f(1:j-1);s;f(j:end)];
        model.(field).VS=[v(1:j-1);left;right;v(j+1:end)];

        % Forward: select birth, zone, eligible layer, uniform split, contrast.
        % Reverse: select death, same zone, one of K new internal interfaces.
        logphi=birth_contrast_logpdf(delta,par,z);
        pb=probabilities(strcmp(names,'birth')); pd=probabilities(strcmp(names,'death'));
        pz=par.proposal.zone_probability(z);
        ptb.log_q_forward=log(pb)+log(pz)-log(E)-log(width)+logphi;
        ptb.log_q_reverse=log(pd)+log(pz)-log(K);
        ptb.log_abs_jacobian=0;
        ptb.index=j; ptb.split_f=s; ptb.delta_vs=delta;
        ptb.K_after(z)=K+1;

    %% =================== DEATH: MERGE ADJACENT PHYSICAL LAYERS ===================
    case 'death'
        if K<=par.mod.kmin(z), ptb.reason='Minimum layer count reached.'; return; end
        if isfield(request,'index'), j=request.index; else, j=randi(K-1); end
        assert(isscalar(j) && j==fix(j) && j>=1 && j<K,'thbi:move','Invalid merge interface.');
        a=edges(j); s=edges(j+1); b=edges(j+2);
        r=(s-a)/(b-a); delta=v(j+1)-v(j);
        merged=r*v(j)+(1-r)*v(j+1);
        remaining_f=f; remaining_f(j)=[];
        model.(field).f=remaining_f(:); % keep a 0-by-1 vector when K becomes 1
        model.(field).VS=[v(1:j-1);merged;v(j+2:end)];

        % Recover exactly the auxiliary variables that the reverse birth uses.
        reverse_edges=[0;model.(field).f;1];
        E=sum(diff(reverse_edges)>2*gap); width=b-a-2*gap;
        if width<=0 || E==0, ptb.reason='Reverse split has zero support.'; return; end
        logphi=birth_contrast_logpdf(delta,par,z);
        pb=probabilities(strcmp(names,'birth')); pd=probabilities(strcmp(names,'death'));
        pz=par.proposal.zone_probability(z);
        ptb.log_q_forward=log(pd)+log(pz)-log(K-1);
        ptb.log_q_reverse=log(pb)+log(pz)-log(E)-log(width)+logphi;
        ptb.log_abs_jacobian=0; % inverse of the unit-Jacobian birth map
        ptb.index=j; ptb.split_f=s; ptb.delta_vs=delta;
        ptb.K_after(z)=K-1;

    %% =================== HIERARCHICAL DATA-ERROR MOVE ===================
    case 'noise'
        nd=numel(par.inv.datatypes);
        if isfield(request,'index'), j=request.index; else, j=randi(nd); end
        assert(isscalar(j) && j==fix(j) && j>=1 && j<=nd,'thbi:move','Invalid noise index.');
        eta=model.datahparm.eta;
        ptb.value_from=eta(j); ptb.sigma=step_scale*par.proposal.noise_step;
        value=eta(j)+ptb.sigma*randn;
        if isfield(request,'value_to'), value=request.value_to; end
        if ~dr_enabled, value=reflect_to_interval(value,par.mod.noise.lower(j),par.mod.noise.upper(j)); end
        eta(j)=value; ptb.value_to=value; ptb.dr_eligible=dr_enabled;
        model.datahparm.eta=eta; model.datahparm.gamma=exp(eta); ptb.index=j;
end

%% =================== SUPPORT, PROPOSAL RATIO AND MODEL ASSEMBLY ===================
% Actual scalar Gaussian densities when DR is enabled. They are symmetric,
% but storing their absolute values lets us audit mixed scalar/RJ paths.
% Disabled-DR reflected moves keep cancelled zero log-ratios as before.
if dr_enabled && ~any(strcmp(ptb.type,{'birth','death'}))
    switch ptb.type
        case 'velocity'
            if ptb.zone==0
                selection=probabilities(itype)*par.proposal.halfspace_probability;
            else
                selection=probabilities(itype)*(1-par.proposal.halfspace_probability)* ...
                    par.proposal.zone_probability(z)/K;
            end
        case 'interface'
            selection=probabilities(itype)*par.proposal.zone_probability(z)/(K-1);
        case 'boundary'
            selection=probabilities(itype)/2;
        case 'noise'
            selection=probabilities(itype)/numel(par.inv.datatypes);
    end
    ptb.log_q_forward=log(selection)-log(sqrt(2*pi)*ptb.sigma) - ...
        .5*((ptb.value_to-ptb.value_from)/ptb.sigma)^2;
    ptb.log_q_reverse=ptb.log_q_forward;
end
% q_forward/reverse always describe the BASE first-stage kernel Q1.
% Actual Q2(y->z) = Q1(y->z)*exp(log_targeting_forward). Global moves
% (noise, boundaries, half-space) have zone=0 and no targeting factor.
if ptb.zone>0
    ptb.log_targeting_forward=log(zone_draw(ptb.zone)/par.proposal.zone_probability(ptb.zone));
end
[ptb.valid,ptb.reason]=a1_TEST_CONDITIONS(model,par);
if any(strcmp(dr_mode,{'layered_path','targeted_layered'}))
    % The reverse path also passes through this SAME valid middle Earth.
    % Out-of-support first trials are ordinary self-transitions, not retried.
    ptb.dr_eligible=dr_enabled && ptb.valid;
end
ptb.log_hastings=ptb.log_q_reverse-ptb.log_q_forward+ptb.log_abs_jacobian;
if ~ptb.valid, return; end % Reject outside support; never clip/redraw a birth.
if ptb.structure, model=make_mod_from_parms(model,par); end
% Symmetric fixed-dimensional moves have zero log proposal ratio. Their
% unchanged move/zone/index selection probabilities cancel in both directions.
end
