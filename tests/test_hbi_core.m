function test_hbi_core()
% Statistical/regression checks with an analytic forward model (no CPS).
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'HBI')); a0_STARTUP_BAYES();
par=bayes_inv_parms('pilot');
par.inv.nchains=2; par.inv.niter=120; par.inv.burnin=40;
par.inv.print_every=1000; par.inv.checkpoint_every=17;
par.mod.start_jitter=0;
root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
par.res.dir=fullfile(root,'likelihood'); mkdir(par.res.dir);
model=make_layer_model(par.mod.start,par);
for k=1:5
    d=par.inv.datatypes{k};
    data.(d)=struct('x',(1:3)','obs',[1;2;3],'sig',[.5;.8;1.1], ...
        'R',[1 .2 .1;.2 1 .3;.1 .3 1]);
    predata.(d)=[.8;2.3;2.7];
end
data=a2_PREPARE_DATA(data,par);
m=b4_CALC_MISFIT(data,predata,par);
[L,~]=b5_CALC_LIKELIHOOD(m,data,model.datahparm,par);
manual=0;
for k=1:5
    d=par.inv.datatypes{k}; C=data.(d).L0*data.(d).L0'; r=data.(d).obs-predata.(d);
    manual=manual-.5*(3*log(2*pi)+log(det(C))+r'*(C\r));
end
assert(abs(L-manual)<1e-11,'Covariance likelihood is wrong.');
e=model.datahparm; e.eta(5)=log(2);
[L2,~]=b5_CALC_LIKELIHOOD(m,data,e,par);
assert(abs((L2-L)-(-3*log(2)+.375*m.E2(5)))<1e-11,'Noise normalization is wrong.');
model.datahparm.eta=[-.2;.3;.7;-1;1.2];
physical_logprior=-.5*sum((model.datahparm.eta/2).^2)-sum(model.datahparm.eta);
assert(abs(calc_log_prior(model,par)-(physical_logprior+sum(model.datahparm.eta)))<1e-12, ...
    'Lognormal transformation Jacobian is wrong.');
bad=data; bad.rf.sig(1)=0;
expect_error(@() a2_PREPARE_DATA(bad,par),'hbi:data');
bad=data; bad.rf.R(1,2)=2;
expect_error(@() a2_PREPARE_DATA(bad,par),'hbi:covariance');

% Target synthetic model for cheap forward tests; independent of true Earth.
model=make_layer_model(par.mod.start,par);
pred=toy_forward(model,data,par,'');
for k=1:5, d=par.inv.datatypes{k}; data.(d).obs=pred.(d); end
par.inv.noise_move_probability=1; par.res.dir=fullfile(root,'noise_only');
c=run_one_chain(par,data,1,[],@toy_forward);
assert(c.nforward==1,'Noise-only updates recomputed the forward model.');
assert(all(all(c.theta==repmat(par.mod.start',par.inv.niter,1))));
assert(any(c.accepted) && any(~c.accepted),'Need both accepted and repeated states.');

% Exact continuation, including adaptation and RNG across interruption.
par.inv.noise_move_probability=.25; par.res.dir=fullfile(root,'continuous');
a=run_one_chain(par,data,1,[],@toy_forward);
par.res.dir=fullfile(root,'resumed');
run_one_chain(par,data,1,53,@toy_forward);
b=run_one_chain(par,data,1,[],@toy_forward);
assert(isequal(a.theta,b.theta) && isequal(a.eta,b.eta) && ...
    isequal(a.logL,b.logL) && isequal(a.accepted,b.accepted) && ...
    isequal(a.final_steps,b.final_steps),'Resumed chain differs from continuous chain.');
assert(a.nforward==1+sum(a.move<=19),'Forward-cache accounting is wrong.');
par.inv.niter=121;
expect_error(@() run_one_chain(par,data,1,[],@toy_forward),'hbi:resume');
par.inv.niter=120;
c=run_one_chain(par,data,2,[],@toy_forward);
assert(~isequal(b.theta,c.theta),'Independent chain seeds produced identical histories.');
assert(exist(fullfile(par.res.dir,'chain_01','checkpoint.mat'),'file')==2 && ...
    exist(fullfile(par.res.dir,'chain_02','checkpoint.mat'),'file')==2);

% Narrow errors force structural rejection; repeated states cannot be dropped.
par.res.dir=fullfile(root,'rejection'); par.inv.noise_move_probability=0;
for k=1:5, d=par.inv.datatypes{k}; data.(d).sig(:)=1e-12; end
data=a2_PREPARE_DATA(data,par);
c=run_one_chain(par,data,1,[],@toy_forward);
assert(~any(c.accepted) && size(c.theta,1)==par.inv.niter);
assert(all(all(c.theta==repmat(par.mod.start',par.inv.niter,1))));
fprintf('PASS: likelihood, covariance, Jacobian, data checks, cache, RNG resume, independent chains, repeated states.\n');
end

function predata=toy_forward(model,data,par,~)
for k=1:5
    d=par.inv.datatypes{k}; predata.(d)=sum(model.theta)*(1:data.(d).N)'/100;
end
end

function expect_error(f,id)
try, f(); catch e, assert(strcmp(e.identifier,id),'Unexpected error identifier.'); return; end
error('hbi:test','Expected %s.',id);
end
