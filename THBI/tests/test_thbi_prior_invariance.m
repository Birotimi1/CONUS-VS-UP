function report=test_thbi_prior_invariance(N,dr_mode)
% One-step invariance of the ACTUAL continuous layered prior/proposals.
% Independent exact prior starts avoid confounding correctness with slow
% mixing of one long random-walk chain. No native solvers are called.
if nargin<1, N=6000; end
assert(N>=1000 && N==fix(N));
par=bayes_inv_parms('extended60'); par.mod.start.method='prior';
if nargin>1, par.proposal.dr.mode=dr_mode; end
par.mod.kmax=[4 5 4];
par.conditions.pos_sed2basement=false; par.conditions.pos_moho=false;
par.conditions.max_moho_fraction=Inf; par.conditions.max_moho_jump=Inf;
par.proposal.dr.same_zone_weight=[.2 .5 .8]; % non-cancelling targeting
saved=rng; cleanup=onCleanup(@() rng(saved)); %#ok<NASGU>
rng(811729,'twister'); validate_configuration(par);
nmetric=sum(par.mod.kmax)+3+3+numel(par.inv.datatypes);
before=zeros(N,nmetric); after=before; stages=zeros(N,1);
for i=1:N
    x.model=b1_INITIATE_MODEL(par); x.logP=a3_CALC_LOG_PRIOR(x.model,par); x.logL=0;
    before(i,:)=metrics(x.model,par);
    [y.model,p1]=b2_PERTURB_MODEL(x.model,par); y.logP=-Inf; y.logL=0;
    if p1.valid, y.logP=a3_CALC_LOG_PRIOR(y.model,par); end
    [accepted,~]=b6_IFACCEPT(x,y,p1);
    result=x;
    if accepted
        result=y; stages(i)=1;
    elseif p1.valid && isfinite(y.logP)
        request=struct();
        if strcmp(par.proposal.dr.mode,'targeted_layered'), request.previous_zone=p1.zone; end
        [z.model,z.ptb]=b2_PERTURB_MODEL(y.model,par,request);
        z.logP=-Inf; z.logL=0;
        if z.ptb.valid, z.logP=a3_CALC_LOG_PRIOR(z.model,par); end
        [accepted,~]=delayed_rejection_layered(x,y,z,p1,par);
        if accepted, result=z; stages(i)=2; end
    end
    after(i,:)=metrics(result.model,par);
    if mod(i,1000)==0, fprintf('Independent prior transitions: %d/%d\n',i,N); end
end
delta=after-before; drift=mean(delta,1); se=std(delta,0,1)/sqrt(N);
score=abs(drift)./max(se,1e-12);
report=struct('dr_mode',par.proposal.dr.mode,'independent_starts',N,'maximum_paired_z',max(score), ...
    'paired_drift',drift,'paired_standard_error',se,'paired_z',score, ...
    'first_accepts',sum(stages==1),'second_accepts',sum(stages==2));
assert(all(score<4.5),'thbi:test','Continuous prior ensemble changed after one RJ/DR step.');
assert(report.first_accepts>0 && report.second_accepts>0);
fprintf('PASS: %d exact-prior starts; largest paired drift %.3f standard errors; stage 1/2 accepts %d/%d.\n', ...
    N,report.maximum_paired_z,report.first_accepts,report.second_accepts);
end

function values=metrics(m,p)
values=zeros(1,sum(p.mod.kmax)); first=0;
for z=1:3
    values(first+(1:p.mod.kmax(z)))=double(m.K(z)==(1:p.mod.kmax(z)));
    first=first+p.mod.kmax(z);
end
for z=1:3
    v=m.(p.mod.zone_fields{z}).VS;
    values(end+1)=(mean(v)-p.mod.vsmin(z))/(p.mod.vsmax(z)-p.mod.vsmin(z)); %#ok<AGROW>
end
values=[values,m.zsed,m.zmoh,m.halfspaceVS,m.datahparm.eta'];
end
