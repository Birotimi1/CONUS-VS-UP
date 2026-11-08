function test_thbi_layered_dr()
% User-run checks for the FULL x->y(rejected)->z layered RJ path.
% No native solvers. All 36 pairs of move families, including two births,
% two deaths, mixed scalar/RJ paths and different zones, are checked below.
par=bayes_inv_parms('extended60'); par.mod.start.jitter=false;
par.proposal.probability=[.20 .10 .10 .25 .10 .25];
par.proposal.zone_probability=[.2 .3 .5];
validate_configuration(par);
saved=rng; cleanup=onCleanup(@() rng(saved)); %#ok<NASGU>
rng(4207,'twister'); base=b1_INITIATE_MODEL(par);

%% =================== TWO-LEG MAPS AND INDEPENDENT REVERSE PATHS ===================
npaths=0;
for setting=1:5
    par.proposal.dr.mode='layered_path';
    if setting>1
        par.proposal.dr.mode='targeted_layered';
        weights={0,.5,.9,[.2 .5 .8]}; par.proposal.dr.same_zone_weight=weights{setting-1};
    end
for first_zone=[1 2]
    for i=1:6
        for j=1:6
            for second_zone=1:3
                [middle,p1]=b2_PERTURB_MODEL(base,par,request(base,par,i,first_zone));
                second_request=request(middle,par,j,second_zone);
                if setting>1, second_request.previous_zone=p1.zone; end
                [last,p2]=b2_PERTURB_MODEL(middle,par,second_request);
                assert(p1.valid && p2.valid && p1.dr_eligible && p2.dr_eligible);
                npaths=npaths+1;
                % Actually generate z->y and y->x with b2, not just negate ratios.
                [reverse_middle,r1]=b2_PERTURB_MODEL(last,par,inverse(p2));
                reverse_request=inverse(p1);
                if setting>1, reverse_request.previous_zone=r1.zone; end
                [reverse_last,r2]=b2_PERTURB_MODEL(reverse_middle,par,reverse_request);
                assert(r1.valid && r2.valid);
                same_model(middle,reverse_middle); same_model(base,reverse_last);
                assert(abs(p2.log_hastings+r1.log_hastings)<1e-9 && ...
                    abs(p1.log_hastings+r2.log_hastings)<1e-9);
                assert(abs(p2.log_q_forward-r1.log_q_reverse)<1e-9 && ...
                    abs(p1.log_q_reverse-r2.log_q_forward)<1e-9);

                %% =================== FULL PATH PROBABILITY, INCLUDING REJECTIONS ===================
                % Non-negligible rejection corrections: alpha_xy=exp(-.4),
                % alpha_zy=exp(-.7). Using an almost-zero alpha would fail to
                % detect accidentally omitted rejection factors at this tolerance.
                lx=0; ly=-.4-p1.log_hastings; lz=ly-p2.log_hastings+.7;
                same_endpoint=isequal(base.K,last.K) && ...
                    max(abs(base.matrix(:)-last.matrix(:)))<1e-10 && ...
                    max(abs(base.datahparm.eta-last.datahparm.eta))<1e-12;
                if same_endpoint, lz=lx; end % a return to x must have pi(z)=pi(x)
                x=target(base,lx,par); y=target(middle,ly,par); z=target(last,lz,par);
                z.ptb=p2;
                [~,a2,terms]=delayed_rejection_layered(x,y,z,p1,par);
                ry=target(reverse_middle,ly,par); rx=target(reverse_last,lx,par); rx.ptb=r2;
                [~,a2reverse]=delayed_rejection_layered(z,ry,rx,r1,par);
                % Independent flux expression from pi, absolute q and alpha1.
                forward=lx+p1.log_q_forward+p2.log_q_forward+p2.log_targeting_forward + ...
                    log1p(-min(1,exp(ly-lx+p1.log_hastings)))+min(0,a2);
                backward=lz+r1.log_q_forward+r2.log_q_forward+r2.log_targeting_forward + ...
                    log1p(-min(1,exp(ly-lz+r1.log_hastings)))+min(0,a2reverse);
                assert(abs(forward-backward-p1.log_abs_jacobian-p2.log_abs_jacobian)<1e-9, ...
                    'thbi:test','Layered DR path balance failed for %s then %s.',p1.type,p2.type);
                assert(abs(a2+a2reverse)<1e-9 && isfinite(terms.log_reject_zy));

                %% =================== ZERO REVERSE-REJECTION SUPPORT ===================
                if ~same_endpoint
                    x=target(base,p1.log_hastings+2,par); y=target(middle,0,par);
                    z=target(last,-p2.log_hastings-1,par); z.ptb=p2;
                    [accepted,ratio]=delayed_rejection_layered(x,y,z,p1,par);
                    assert(~accepted && ratio==-Inf,'thbi:test','Impossible reverse rejection was accepted.');
                end
            end
        end
    end
end
end

%% =================== SCALAR DENSITY AND INVALID PATH GUARDS ===================
[middle,p1]=b2_PERTURB_MODEL(base,par,request(base,par,1,2));
sigma=par.proposal.vs_step(2); delta=p1.value_to-p1.value_from;
density=par.proposal.probability(1)*(1-par.proposal.halfspace_probability)* ...
    par.proposal.zone_probability(2)/base.K(2) * exp(-.5*(delta/sigma)^2)/(sqrt(2*pi)*sigma);
assert(abs(p1.log_q_forward-log(density))<1e-12);
req=request(middle,par,4,2); req.previous_zone=p1.zone;
[last,p2]=b2_PERTURB_MODEL(middle,par,req);
x=target(base,0,par); y=target(middle,-40,par); z=target(last,-1,par); z.ptb=p2;
bad=z; bad.ptb.step_scale=.25;
expect_error(@() delayed_rejection_layered(x,y,bad,p1,par),'thbi:dr');
bad=z; bad.ptb.K_before=base.K+[1 0 0];
expect_error(@() delayed_rejection_layered(x,y,bad,p1,par),'thbi:dr');
bad=z; bad.ptb.log_q_reverse=NaN;
expect_error(@() delayed_rejection_layered(x,y,bad,p1,par),'thbi:dr');
bad=z; bad.logL=-Inf;
[accepted,ratio]=delayed_rejection_layered(x,y,bad,p1,par);
assert(~accepted && ratio==-Inf);
limited=par; limited.mod.kmax=base.K;
[~,invalid]=b2_PERTURB_MODEL(base,limited,struct('type','birth','zone',2));
assert(~invalid.valid && ~invalid.dr_eligible,'thbi:test','Impossible first split must not trigger DR.');
% Reverse first-stage evaluations must not change RNG state.
before=rng; [~,first_ratio]=b6_IFACCEPT(x,y,p1,false); after=rng;
assert(isequal(before,after));
[~,~,t]=delayed_rejection_layered(x,y,z,p1,par);
assert(t.logalpha_xy==min(0,first_ratio));
assert(npaths==1080);
fprintf('PASS: %d layered two-stage paths, inverse maps, absolute densities, path balance and support guards.\n',npaths);
end

function r=request(model,par,type,z)
r=struct('type',par.proposal.names{type},'zone',z,'index',1);
zone=model.(par.mod.zone_fields{z}); edges=[0;zone.f;1];
switch r.type
    case 'velocity'
        r.halfspace=false; r.value_to=zone.VS(1)+.005;
    case 'interface'
        r.value_to=zone.f(1)+.001;
    case 'boundary'
        r.index=2; r.value_to=model.zmoh+.01;
    case 'birth'
        r.split_f=edges(2)/2; r.delta_vs=.02;
    case 'noise'
        r.value_to=model.datahparm.eta(1)+.02;
end
end

function r=inverse(p)
r=struct('type',p.type,'zone',p.zone,'index',p.index);
switch p.type
    case 'birth'
        r.type='death';
    case 'death'
        r.type='birth'; r.split_f=p.split_f; r.delta_vs=p.delta_vs;
    otherwise
        r.value_to=p.value_from;
        if strcmp(p.type,'velocity')
            r.halfspace=p.zone==0;
            if r.halfspace, r.zone=1; end
        end
end
end

function state=target(model,logtarget,par)
prior=a3_CALC_LOG_PRIOR(model,par);
state=struct('model',model,'logP',prior,'logL',logtarget-prior);
end

function same_model(a,b)
assert(isequal(a.K,b.K) && max(abs(a.matrix(:)-b.matrix(:)))<1e-9 && ...
    max(abs(a.datahparm.eta-b.datahparm.eta))<1e-12,'thbi:test','Reverse path did not recover the Earth.');
end

function expect_error(callback,id)
try
    callback(); error('thbi:test','Expected error was not raised.');
catch problem
    assert(strcmp(problem.identifier,id),'thbi:test','Unexpected error: %s',problem.message);
end
end
