function test_thbi_reference100()
% Solver-free checks of the reference-data formulation and extended RJ moves.
par=bayes_inv_parms('reference100'); validate_configuration(par);
saved=rng; cleanup=onCleanup(@() rng(saved)); %#ok<NASGU>
rng(9127,'twister'); [~,truth,forward]=a2_LOAD_DATA(par,true);
assert(strcmp(forward.status,'not_computed'));
base=truth.zoned_parameters;
assert(isequal(base.K,[1 3 5]) && abs(base.zsed-.01)<1e-10 && ...
    abs(base.zmoh-43.7731)<1e-9 && abs(base.halfspaceVS-4.5912)<1e-10);
assert(all(diff(base.mantmparm.VS)>0) && a1_TEST_CONDITIONS(base,par));
assert(max(abs(base.crustmparm.VS-truth.published_ecm1_layers(2:4,4)))<1e-12);
assert(all(abs(base.matrix(2:4,2)./base.matrix(2:4,3)-base.vpvs)<1e-12));
lvz=base; lvz.crustmparm.VS(2)=3.35; lvz.mantmparm.VS(3)=4.1;
assert(a1_TEST_CONDITIONS(lvz,par),'Crustal and mantle LVZs should be supported.');
bad=base; bad.vpvs=2.2; assert(~a1_TEST_CONDITIONS(bad,par));
for j=1:20
    m=b1_INITIATE_MODEL(par); assert(a1_TEST_CONDITIONS(m,par));
end

%% =================== ALL 100 PAIRS, INCLUDING NEW MOVES AND BIRTH/DEATH ===================
% Two consecutive crustal deaths in the actual truth would violate the
% 0.8 km/s Moho jump bound. Check rejection, then use a supported path fixture.
[merged,~]=b2_PERTURB_MODEL(base,par,struct('type','death','zone',2,'index',1));
[~,invalid]=b2_PERTURB_MODEL(merged,par,struct('type','death','zone',2,'index',1));
assert(~invalid.valid && ~invalid.dr_eligible);
reference_model=base; base.crustmparm.VS=[3.6;3.7;3.9]; base=make_mod_from_parms(base,par);
count=0;
for i=1:10
    for j=1:10
        [middle,p1]=b2_PERTURB_MODEL(base,par,request(base,par,i));
        [last,p2]=b2_PERTURB_MODEL(middle,par,request(middle,par,j));
        assert(p1.valid && p2.valid && p1.dr_eligible && p2.dr_eligible, ...
            'thbi:test','Moves %s/%s invalid: %s / %s',p1.type,p2.type,p1.reason,p2.reason);
        [reverse_middle,r1]=b2_PERTURB_MODEL(last,par,inverse(p2));
        [back,r2]=b2_PERTURB_MODEL(reverse_middle,par,inverse(p1));
        assert(r1.valid && r2.valid);
        assert(max(abs(back.matrix(:)-base.matrix(:)))<1e-9 && abs(back.vpvs-base.vpvs)<1e-12);
        assert(abs(p2.log_q_forward-r1.log_q_reverse)<1e-9 && ...
            abs(p1.log_q_reverse-r2.log_q_forward)<1e-9);
        lx=0; ly=-.4-p1.log_hastings; lz=ly-p2.log_hastings+.7;
        if isequal(base.K,last.K) && max(abs(base.matrix(:)-last.matrix(:)))<1e-10 && ...
                max(abs(base.datahparm.eta-last.datahparm.eta))<1e-12
            lz=lx;
        end
        x=target(base,lx,par); y=target(middle,ly,par); z=target(last,lz,par); z.ptb=p2;
        [~,a]=delayed_rejection_layered(x,y,z,p1,par);
        ry=target(reverse_middle,ly,par); rx=target(back,lx,par); rx.ptb=r2;
        [~,b]=delayed_rejection_layered(z,ry,rx,r1,par);
        flux1=lx+p1.log_q_forward+p2.log_q_forward+log1p(-min(1,exp(ly-lx+p1.log_hastings)))+min(0,a);
        flux2=lz+r1.log_q_forward+r2.log_q_forward+log1p(-min(1,exp(ly-lz+r1.log_hastings)))+min(0,b);
        assert(abs(flux1-flux2)<1e-9 && abs(a+b)<1e-9);
        count=count+1;
    end
end
base=reference_model;

%% =================== INDEPENDENT H-KAPPA ARRIVAL-TIME FIXTURE ===================
H=base.zmoh; K=base.vpvs; h=base.matrix(1:4,1); vs=H/sum(h./base.matrix(1:4,3));
t=(-100:.2:500)'; p=par.rf.hk.rayp_s_km; rf=zeros(numel(t),numel(p));
for j=1:numel(p)
    qs=sqrt(1/vs^2-p(j)^2); qp=sqrt(1/(K*vs)^2-p(j)^2);
    arrivals=H*[qs-qp qs+qp 2*qs]; amplitudes=[1 .6 -.4];
    for k=1:3, rf(:,j)=rf(:,j)+amplitudes(k)*exp(-.5*((t-arrivals(k))/.35).^2); end
end
observed.waves=struct('t',t,'rf',rf,'rayp_s_km',p);
stack=b3_FORWARD_MODEL_HK(base,observed,par);
assert(stack.E_by_Emax>.99,'H-kappa phase signs/times do not recover the fixture.');
wrong=base; wrong.zmoh=H+3; wrong=make_mod_from_parms(wrong,par);
wrong_stack=b3_FORWARD_MODEL_HK(wrong,observed,par);
assert(wrong_stack.E_by_Emax<stack.E_by_Emax-.1);

%% =================== SIGNED Sp ROTATION, MIGRATION AND RESIDUAL ===================
cfg=par.forc.rf_sp; time=(0:cfg.npts-1)'*cfg.dt-cfg.delay;
lag=layered_conversion_time(base.matrix,100,cfg.rayp_s_km,'Sp');
assert(lag<0);
sigma=par.rf.ccp.source_period_s/6;
daughter=.08*exp(-.5*((time-lag)/sigma).^2); parent=exp(-.5*(time/sigma).^2);
[P,SV]=Rotate_XZ_to_PSV([1;0],[0;1],base.VP(1),base.VS(1),cfg.rayp_s_km);
components=[P SV]'\[daughter parent]';
r=struct('t',time,'dt',cfg.dt,'val',components(1,:)');
z=r; z.val=-components(2,:)'; % convert Z-down fixture to CPS Z-up
ccp=b3_FORWARD_MODEL_RF_CCP(base.matrix,{r,z},par);
[~,peak]=max(ccp.val); assert(abs(ccp.x(peak)-100)<=1);
zero=ccp_cross_convolution_misfit(ccp.PSV,ccp.PSV,ccp.x,par.rf.ccp);
assert(abs(zero)<1e-12);
scaled=ccp.PSV; scaled(:,1)=2*scaled(:,1);
energy=ccp_cross_convolution_misfit(ccp.PSV,scaled,ccp.x,par.rf.ccp);
assert(abs(energy-1/sqrt(2))<1e-12,'Reference RMS ratio was changed to a squared norm.');

%% =================== REFERENCE RF LIKELIHOOD AND NOISE NORMALIZATION ===================
for d=1:4
    name=par.inv.datatypes{d}; data.(name)=struct('x',[1;2;3],'obs',[1;2;3], ...
        'sig',ones(3,1)); pre.(name)=[1;2;3];
end
data.rf=struct('x',1,'obs',1,'N',1,'waves',observed.waves,'dof',15,'nrecords',1);
data.rf_sp=struct('x',ccp.x,'obs',ccp.val,'parent',ccp.parent,'N',numel(ccp.x),'dof',43,'nrecords',1);
data=a2_PREPARE_DATA(data,par); pre.rf=.8; pre.rf_sp=scaled(:,1); pre.rf_sp_parent=scaled(:,2);
misfit=b4_CALC_MISFIT(data,pre,par);
assert(abs(misfit.E2(5)-.202)<1e-12 && abs(misfit.E2(6)-1/sqrt(2))<1e-12);
[L,m]=b5_CALC_LIKELIHOOD(misfit,data,base.datahparm,par);
e=base.datahparm.eta(5); expected=-15*e-.5*misfit.E2(5)*exp(-2*e);
assert(abs(m.logL_indivdat(5)-expected)<1e-12);
noise=base.datahparm; noise.eta(5)=noise.eta(5)+log(2);
L2=b5_CALC_LIKELIHOOD(misfit,data,noise,par);
assert(abs(L2-L-(-15*log(2)+.375*misfit.E2(5)*exp(-2*e)))<1e-10);

%% =================== DIAGNOSTICS: IID, STICKY AND DISAGREEING CHAINS ===================
iid=randn(2000,1,4); good=mcmc_rank_diagnostics(iid);
sticky=filter(1,[1 -.97],iid,[],1); poor=mcmc_rank_diagnostics(sticky);
shifted=iid; shifted(:,:,4)=shifted(:,:,4)+4; separate=mcmc_rank_diagnostics(shifted);
assert(good.rhat<1.02 && good.bulk_ess>4000 && good.tail_ess>3000);
assert(poor.bulk_ess<good.bulk_ess/10 && separate.rhat>1.2);
fprintf('PASS: ECM1/AK135 provenance, LVZ support, %d full DR paths, HK timings, Sp signs/migration, cross-convolution, RF likelihood and rank/ESS controls.\n',count);
end

function r=request(m,p,i)
r=struct('type',p.proposal.names{i},'zone',2,'index',1);
switch r.type
    case 'velocity', r.halfspace=false; r.value_to=m.crustmparm.VS(1)+.002;
    case 'interface', r.value_to=m.crustmparm.f(1)+.001;
    case 'boundary', r.index=2; r.value_to=m.zmoh+.01;
    case 'birth', r.split_f=m.crustmparm.f(1)/2; r.delta_vs=.005;
    case 'noise', r.index=5; r.value_to=m.datahparm.eta(5)+.02;
    case 'vpvs', r.value_to=m.vpvs+.002;
    case 'boundary_mean', r.index=2; r.value_to=mean([m.crustmparm.VS(end),m.mantmparm.VS(1)])+.002;
    case 'boundary_jump', r.index=2; r.value_to=m.mantmparm.VS(1)-m.crustmparm.VS(end)+.002;
    case 'moho_hk', r.value_to=[m.zmoh+.01,m.vpvs+.001];
end
end
function r=inverse(p)
r=struct('type',p.type,'zone',p.zone,'index',p.index);
if strcmp(p.type,'birth'), r.type='death';
elseif strcmp(p.type,'death'), r.type='birth'; r.split_f=p.split_f; r.delta_vs=p.delta_vs;
else, r.value_to=p.value_from; if strcmp(p.type,'velocity'), r.halfspace=false; end
end
end
function s=target(m,logtarget,p)
prior=a3_CALC_LOG_PRIOR(m,p); s=struct('model',m,'logP',prior,'logL',logtarget-prior);
end
