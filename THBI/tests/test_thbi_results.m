function test_thbi_results()
% Manual MATLAB checks for retained weights, interfaces, counts and summaries.
% No native solvers or MCMC. Constructs a short deterministic saved history.
par=bayes_inv_parms('pilot'); par.mod.start.jitter=false;
par.inv.niter=5; par.inv.burnin=2; par.inv.nchains=1;
par.res.depth_km=[0;.02;.25;36;70;85];
base=b1_INITIATE_MODEL(par);
[child,ptb]=b2_PERTURB_MODEL(base,par,struct('type','birth','zone',1, ...
    'index',1,'split_f',.05,'delta_vs',.02));
assert(ptb.valid,'Test split should be feasible.');
[N,K,Nall]=count_N_layers_in_mod(child,par);
assert(N==sum(base.K)+1 && Nall==N+1 && K(1)==base.K(1)+1);
chain=struct('fractions',{cell(5,3)},'velocities',{cell(5,3)},'K',zeros(5,3), ...
    'boundaries',zeros(5,2),'halfspace',zeros(5,1),'eta',zeros(5,5), ...
    'logL',zeros(5,1),'E2',ones(5,5),'last_iter',5);
for i=1:5
    model=base; if i==5, model=child; end
    for z=1:3
        zone=model.(par.mod.zone_fields{z});
        chain.fractions{i,z}=zone.f; chain.velocities{i,z}=zone.VS;
    end
    chain.K(i,:)=model.K; chain.boundaries(i,:)=[model.zsed model.zmoh];
    chain.halfspace(i)=model.halfspaceVS; chain.eta(i,:)=model.datahparm.eta';
end
chain.E2(1,:)=0; % burn-in has a better fit, but must not enter final selection
chain.E2(5,:)=2;
[posterior,retained]=c2_BUILD_POSTERIOR({chain},par);
assert(posterior.Nstored==3 && isequal(posterior.iteration,[3;4;5]));
assert(all(posterior.chain==1),'Single-chain collation lost identity.');
count_index=find(posterior.count_values{1}==base.K(1));
assert(abs(posterior.count_probability{1}(count_index)-2/3)<1e-12, ...
    'Repeated retained states were lost or reweighted.');
suite=c3_BUILD_MODEL_SUITE(retained,par);
assert(isequal(suite.VS(:,1),suite.VS(:,2)),'Identical retained states disagree.');
final_model=c4_FINAL_MODEL(posterior,suite,retained,par);
assert(final_model.representative_iteration==3 && final_model.representative_chain==1);
assert(isequal(final_model.representative_model.matrix,base.matrix));
assert(isequal(final_model.VSquantiles,sample_quantiles(suite.VS',[.05 .5 .95])'));

%% Summary/fit plots must use the SAME retained Earth, even if burn-in fits better
par.inv.prior_only=true; par.res.prediction_stride=0; par.res.dir=tempname;
state=struct('model',base,'predata',struct(),'misfit',struct('E2',zeros(5,1)), ...
    'logL',0,'logP',a3_CALC_LOG_PRIOR(base,par));
chain.initial=state; chain.bestfit=state;
state.misfit.E2=ones(5,1); chain.bestfit_retained=state;
chain.bestfit_retained_iteration=3;
chain.move=zeros(5,1); chain.zone=zeros(5,1);
chain.accepted=false(5,1); chain.proposal_valid=false(5,1); chain.forward_failed=false(5,1);
% Last stored Earth came from a rejected noise step followed by a birth.
chain.move(5)=6; chain.accepted(5)=true; chain.proposal_valid(5)=true;
chain.accepted_stage=[0;0;0;0;2]; chain.first_proposal_valid=chain.proposal_valid;
chain.dr_attempted=chain.accepted; chain.dr_proposal_valid=chain.accepted;
chain.dr_move=[0;0;0;0;4]; chain.dr_zone=[0;0;0;0;1];
summary=c1_PROCESS_RESULTS({chain},struct(),par);
assert(sum(summary.move_accepts(:))==0 && summary.dr_move_accepts(1,4)==1 && ...
    summary.dr_zone_move_accepts(1,4,1)==1 && summary.dr_attempts(1)==1, ...
    'thbi:test','Second-stage birth was counted as a first-stage noise acceptance.');
assert(summary.bestfit_iteration==final_model.representative_iteration && ...
    summary.bestfit_chain==final_model.representative_chain);
assert(sum(summary.bestfit.misfit.E2)==5 && isequal(summary.bestfit.model.matrix,base.matrix), ...
    'thbi:test','Burn-in best fit leaked into the retained-fit plot.');
legacy=rmfield(chain,{'bestfit_retained','bestfit_retained_iteration'});
old=c1_PROCESS_RESULTS({legacy},struct(),par);
assert(isnan(old.bestfit_iteration) && contains(old.bestfit_label,'legacy'), ...
    'thbi:test','Legacy cached fits must not be mislabeled as retained.');

%% Exact interface values and the half-space: no linear interpolation
matrix=[1 4 2 2;0 6 3.5 2.7];
values=sample_layer_model(matrix,[0;.999;1;1000]);
assert(isequal(values,matrix([1 1 2 2],2:4)));
chain.last_iter=4;
try
    collate_allmodels_perchain({chain},par);
    error('thbi:test','An incomplete chain was accepted.');
catch problem
    assert(strcmp(problem.identifier,'thbi:results'));
end
fprintf('PASS: retained weights/counts, step evaluation, retained fit/final-model agreement and legacy labeling.\n');
end
