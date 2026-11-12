function test_thbi_reference_chain(profile)
% Short injected-forward test of reference RF cache, kappa storage and restart.
% Analytic responses below are TEST FIXTURES, never scientific predictions.
if nargin<1, profile='reference100'; end
par=bayes_inv_parms(profile); par.inv.niter=240; par.inv.burnin=80;
nmoves=numel(par.proposal.names);
par.inv.print_every=1000; par.inv.checkpoint_every=37;
root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
for d=1:4
    name=par.inv.datatypes{d}; data.(name)=struct('x',[1;2;3],'obs',[1;2;3],'sig',ones(3,1));
end
data.rf=struct('x',1,'N',1,'obs',1,'dof',15,'nrecords',1, ...
    'waves',struct('t',[-1;0;1],'rf',[0;1;0],'rayp_s_km',.06));
z=par.rf.ccp.depth_km; parent=exp(-.5*(z/3).^2);
data.rf_sp=struct('x',z,'N',numel(z),'obs',exp(-.5*((z-80)/10).^2), ...
    'parent',parent,'dof',43,'nrecords',1);
data=a2_PREPARE_DATA(data,par);
par.res.dir=fullfile(root,'continuous'); a=run_one_chain(par,data,1,[],@toy_forward);
par.res.dir=fullfile(root,'resumed');
initial=run_one_chain(par,data,1,0,@toy_forward);
assert(initial.last_iter==0 && initial.nforward==1);
assert(exist(fullfile(par.res.dir,'chain_01','progress.txt'),'file')==2);
paused=run_one_chain(par,data,1,[],@toy_forward,struct('max_wall_seconds',eps));
assert(paused.last_iter==0);
run_one_chain(par,data,1,91,@toy_forward);
saved=load(fullfile(par.res.dir,'chain_01','checkpoint.mat'),'checkpoint');
assert(isscalar(saved.checkpoint.spec) && ...
    isequal(saved.checkpoint.spec.source,source_snapshot()), ...
    'thbi:test','Checkpoint must preserve the complete source cell in one scalar spec.');
% The native path attaches a build signature after constructing the spec.
spec=saved.checkpoint.spec; spec.native_signature=struct('test_fixture',true);
assert(isscalar(spec));
b=run_one_chain(par,data,1,[],@toy_forward);
for field={'fractions','velocities','K','boundaries','halfspace','vpvs','eta','logL','logP','E2', ...
        'bestfit_retained','accepted_stage','dr_move','dr_zone','dr_logalpha', ...
        'first_forward_called','dr_forward_called','move','accepted','predictions'}
    assert(isequaln(a.(field{1}),b.(field{1})),'thbi:test','Reference restart mismatch: %s',field{1});
end
assert(any(b.accepted_stage==2) && all(ismember(1:nmoves,unique(b.move))));
if isfield(b,'xi'), assert(isequal(a.xi,b.xi) && any(diff(b.xi(:,2))~=0)); end
for i=1:par.inv.niter
    model=model_from_chain(b,i,par); prediction=toy_forward(model,data,par,'');
    assert(a1_TEST_CONDITIONS(model,par) && model.vpvs==b.vpvs(i));
    if isfield(b,'xi'), assert(isequal(model.xi,b.xi(i,:))); end
    m=b4_CALC_MISFIT(data,prediction,par); L=b5_CALC_LIKELIHOOD(m,data,model.datahparm,par);
    assert(max(abs(m.E2'-b.E2(i,:)))<1e-10 && abs(L-b.logL(i))<1e-9);
end
c=run_one_chain(par,data,2,[],@toy_forward);
assert(~isequal(a.initial.model.matrix,c.initial.model.matrix));
[posterior,retained]=c2_BUILD_POSTERIOR({a,c},par);
assert(numel(posterior.vpvs)==2*(par.inv.niter-par.inv.burnin));
if isfield(b,'xi'), assert(isequal(size(posterior.xi),[320 3])); end
suite=c3_BUILD_MODEL_SUITE(retained,par); final=c4_FINAL_MODEL(posterior,suite,retained,par);
assert(isfield(final.representative_model,'vpvs') && final.representative_iteration>par.inv.burnin);
if exist('histcounts','file') || exist('histcounts','builtin')
    summary=c1_PROCESS_RESULTS({a,c},data,par);
    assert(size(summary.move_attempts,2)==nmoves && ...
        numel(summary.metric_names)==13+isfield(b,'xi'));
    assert(summary.bestfit_iteration==final.representative_iteration && summary.bestfit_chain==final.representative_chain);
end
par.proposal.probability=zeros(1,nmoves); par.proposal.probability(6)=1;
par.res.dir=fullfile(root,'noise'); noise=run_one_chain(par,data,1,[],@toy_forward);
assert(noise.nforward==1 && all(noise.vpvs==noise.initial.model.vpvs));
assert(all(all(noise.E2==repmat(noise.E2(1,:),par.inv.niter,1))));
fprintf('PASS: reference RF per-state likelihood/cache, all %d moves, parameter storage, exact restart, independent starts and retained reconstruction.\n',nmoves);
end

function pre=toy_forward(m,data,par,~)
v=sum(m.matrix(:,1).*m.matrix(:,3))/par.mod.maxz;
for d=1:4
    name=par.inv.datatypes{d}; pre.(name)=v*(1:data.(name).N)'/10;
end
pre.rf=exp(-((m.zmoh-45)/15)^2-((m.vpvs-1.76)/.2)^2);
pre.rf_sp=(m.mantmparm.VS(end)/4.5)*exp(-.5*((data.rf_sp.x-(m.zmoh+40))/10).^2);
pre.rf_sp_parent=data.rf_sp.parent;
if isfield(m,'xi')
    pre.rf=pre.rf*exp(-((m.xi(2)-1.06)/.1)^2);
    pre.cph=pre.cph*sqrt(m.xi(2));
end
end
