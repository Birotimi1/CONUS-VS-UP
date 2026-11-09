function chain = run_one_chain(par, data, ichain, stop_iter, forward_fun)
% One zoned RJ-HBI chain, with exact RNG/checkpoint continuation.
% forward_fun is injectable for implementation checks. No adaptation or
% reheating changes the target/proposal during this initial RJ implementation.
if nargin<4 || isempty(stop_iter), stop_iter=par.inv.niter; end
if nargin<5, forward_fun=@b3_FORWARD_MODEL; end
validate_configuration(par);
nd=numel(par.inv.datatypes);
layered_dr=isfield(par.proposal,'dr') && par.proposal.dr.enabled && ...
    any(strcmp(par.proposal.dr.mode,{'layered_path','targeted_layered'}));
assert(isscalar(ichain) && isreal(ichain) && isfinite(ichain) && ...
    ichain==fix(ichain) && ichain>=1 && ichain<=par.inv.nchains, ...
    'thbi:chain','Chain index must be an integer in 1:nchains.');
assert(isscalar(stop_iter) && isreal(stop_iter) && isfinite(stop_iter) && ...
    stop_iter>=1 && stop_iter==fix(stop_iter), ...
    'thbi:chain','stop_iter must be a positive integer.');
stop_iter=min(stop_iter,par.inv.niter);

%% =================== INDEPENDENT CHAIN FILES AND LOCAL SOLVER SCRATCH ===================
outdir=fullfile(par.res.dir,sprintf('chain_%02d',ichain));
if ~exist(outdir,'dir'), mkdir(outdir); end
filename=fullfile(outdir,'checkpoint.mat');
spec=struct('schema','zoned_splitmerge_dr_v8_vti','par',par,'data',data, ...
    'forward',func2str(forward_fun),'source',source_snapshot());
if strcmp(func2str(forward_fun),'b3_FORWARD_MODEL') && ~par.inv.prior_only && ...
        isfield(par.mod,'xi_bounds')
    spec.native_signature=check_native_vti_build(par.forc);
end
workdir=tempname(par.res.scratch_root); mkdir(workdir);
cleanup=onCleanup(@() rmdir(workdir,'s')); %#ok<NASGU>
clock0=tic;

%% =================== RESTORE OR INITIALIZE ===================
if exist(filename,'file')
    loaded=load(filename,'checkpoint'); cp=loaded.checkpoint;
    if ~isequaln(cp.spec,spec)
        error('thbi:resume','Code/configuration/data changed. Start a new run instead of reusing this checkpoint.');
    end
    rng(cp.rng); state=cp.state; chain=cp.chain; first=cp.iter+1;
else
    rng(par.inv.seed+1009*ichain,'twister');
    state.model=b1_INITIATE_MODEL(par);
    state.logP=a3_CALC_LOG_PRIOR(state.model,par);
    if par.inv.prior_only
        state.predata=struct(); state.logL=0;
        state.misfit=struct('E2',zeros(nd,1),'rms',zeros(nd,1),'logL_indivdat',zeros(nd,1),'chi2',zeros(nd,1));
    else
        try
            state.predata=forward_fun(state.model,data,par,workdir);
            state.misfit=b4_CALC_MISFIT(data,state.predata,par);
            [state.logL,state.misfit]=b5_CALC_LIKELIHOOD(state.misfit,data,state.model.datahparm,par);
        catch problem
            context=struct('stage','initial model','iteration',0);
            preserve_forward_failure(fullfile(outdir,'initial_forward_failure'),state.model,data,par,problem,workdir,context);
            rethrow(problem);
        end
    end
    assert(isfinite(state.logP+state.logL),'thbi:start','Nonfinite starting posterior.');
    n=par.inv.niter;
    chain=struct();
    chain.fractions=cell(n,3); chain.velocities=cell(n,3);
    chain.K=NaN(n,3); chain.boundaries=NaN(n,2); chain.halfspace=NaN(n,1);
    if isfield(par.mod,'vpvs_bounds'), chain.vpvs=NaN(n,1); end
    if isfield(par.mod,'xi_bounds'), chain.xi=NaN(n,3); end
    chain.eta=NaN(n,nd); chain.logL=NaN(n,1); chain.logP=NaN(n,1);
    chain.E2=NaN(n,nd); chain.logL_indivdat=NaN(n,nd);
    chain.accepted=false(n,1); chain.proposal_valid=false(n,1);
    chain.forward_failed=false(n,1);
    chain.move=zeros(n,1); chain.zone=zeros(n,1);
    chain.logalpha=NaN(n,1); chain.log_q_ratio=NaN(n,1); chain.log_jacobian=NaN(n,1);
    chain.accepted_stage=zeros(n,1); chain.dr_attempted=false(n,1);
    chain.dr_proposal_valid=false(n,1);
    chain.first_proposal_valid=false(n,1); chain.first_logalpha=NaN(n,1);
    chain.dr_logalpha=NaN(n,1);
    chain.dr_mode=par.proposal.dr.mode;
    chain.dr_move=zeros(n,1); chain.dr_zone=zeros(n,1);
    chain.dr_K_before=NaN(n,3); chain.dr_K_after=NaN(n,3);
    chain.dr_log_q_ratio=NaN(n,1); chain.dr_log_jacobian=NaN(n,1);
    chain.dr_log_targeting_ratio=NaN(n,1);
    chain.dr_reverse_first_logalpha=NaN(n,1);
    chain.first_K_after=NaN(n,3);
    chain.first_forward_called=false(n,1); chain.dr_forward_called=false(n,1);
    % Save a fixed, iteration-based subsample of cached responses, including
    % repeated states. This costs no extra native forward calls.
    chain.prediction_iterations=zeros(0,1); chain.predictions=struct();
    if ~par.inv.prior_only && par.res.prediction_stride>0
        chain.prediction_iterations=(par.inv.burnin+1:par.res.prediction_stride:n)';
        for d=1:numel(par.inv.datatypes)
            name=par.inv.datatypes{d};
            chain.predictions.(name)=NaN(numel(chain.prediction_iterations),data.(name).N);
        end
    end
    chain.prediction_saved=false(numel(chain.prediction_iterations),1);
    chain.initial=state; chain.bestfit=state;
    chain.bestfit_iteration=0;
    chain.bestfit_retained=[]; chain.bestfit_retained_iteration=NaN;
    chain.nforward=double(~par.inv.prior_only); chain.failures=0;
    chain.consecutive_failures=0; chain.last_failure=''; chain.last_proposal_reason='';
    chain.elapsed_seconds=0; first=1;
end

%% =================== PROPOSE / FORWARD / ACCEPT / STORE ===================
for iter=first:stop_iter
    trial=struct('logP',-Inf,'logL',-Inf);
    [trial.model,ptb]=b2_PERTURB_MODEL(state.model,par);
    chain.last_proposal_reason=ptb.reason;
    if ptb.valid
        trial.logP=a3_CALC_LOG_PRIOR(trial.model,par);
        if isfinite(trial.logP)
            if par.inv.prior_only
                trial.predata=state.predata; trial.misfit=state.misfit; trial.logL=0;
            elseif ptb.structure
                chain.nforward=chain.nforward+1; chain.first_forward_called(iter)=true;
                try
                    trial.predata=forward_fun(trial.model,data,par,workdir);
                    trial.misfit=b4_CALC_MISFIT(data,trial.predata,par);
                    [trial.logL,trial.misfit]=b5_CALC_LIKELIHOOD(trial.misfit,data,trial.model.datahparm,par);
                    chain.consecutive_failures=0;
                catch problem
                    context=struct('stage','proposal','iteration',iter,'proposal',ptb);
                    if ~strcmp(problem.identifier,'thbi:forward_failure')
                        preserve_forward_failure(fullfile(outdir,'fatal_forward_failure'),trial.model,data,par,problem,workdir,context);
                        rethrow(problem);
                    end
                    preserve_forward_failure(fullfile(outdir,'first_forward_failure'),trial.model,data,par,problem,workdir,context);
                    chain.failures=chain.failures+1;
                    chain.forward_failed(iter)=true;
                    chain.consecutive_failures=chain.consecutive_failures+1;
                    chain.last_failure=problem.message;
                    if chain.consecutive_failures>=par.inv.max_consecutive_failures
                        preserve_forward_failure(fullfile(outdir,'fatal_forward_failure'),trial.model,data,par,problem,workdir,context);
                        rethrow(problem);
                    end
                end
            else
                % Hierarchical noise changes no structural forward prediction.
                trial.predata=state.predata; trial.misfit=state.misfit;
                [trial.logL,trial.misfit]=b5_CALC_LIKELIHOOD(trial.misfit,data,trial.model.datahparm,par);
            end
        end
    end
    [accepted,logalpha]=b6_IFACCEPT(state,trial,ptb);
    chain.first_proposal_valid(iter)=ptb.valid; chain.first_logalpha(iter)=logalpha;
    if accepted, chain.accepted_stage(iter)=1; end

    chain.first_K_after(iter,:)=ptb.K_after;
    %% =================== SECOND PROPOSAL AND REVERSE-PATH CORRECTION ===================
    % Layered mode mirrors Zach's two-try structure: propose z FROM rejected y
    % with the full b2 kernel, including split/merge. The smaller scalar mode
    % remains selectable for comparison. Native failures never trigger DR.
    if ~accepted && ptb.dr_eligible && ~chain.forward_failed(iter) && ...
            (~layered_dr || isfinite(trial.logP+trial.logL))
        chain.dr_attempted(iter)=true;
        second=struct('logP',-Inf,'logL',-Inf);
        if layered_dr
            base=trial;
            if strcmp(par.proposal.dr.mode,'targeted_layered')
                request=struct('previous_zone',ptb.zone);
                [second.model,second.ptb]=b2_PERTURB_MODEL(base.model,par,request);
            else
                [second.model,second.ptb]=b2_PERTURB_MODEL(base.model,par);
            end
        else
            base=state;
            request=struct('type',ptb.type,'index',ptb.index, ...
                'step_scale',par.proposal.dr.scale);
            if strcmp(ptb.type,'velocity') && ptb.zone==0
                request.halfspace=true; request.zone=1;
            elseif any(strcmp(ptb.type,{'velocity','interface'}))
                request.zone=ptb.zone;
            end
            [second.model,second.ptb]=b2_PERTURB_MODEL(base.model,par,request);
        end
        chain.dr_move(iter)=second.ptb.type_index; chain.dr_zone(iter)=second.ptb.zone;
        chain.dr_K_before(iter,:)=second.ptb.K_before;
        chain.dr_K_after(iter,:)=second.ptb.K_after;
        chain.dr_log_q_ratio(iter)=second.ptb.log_q_reverse-second.ptb.log_q_forward;
        chain.dr_log_jacobian(iter)=second.ptb.log_abs_jacobian;
        if second.ptb.valid
            second.logP=a3_CALC_LOG_PRIOR(second.model,par);
            if isfinite(second.logP)
                if par.inv.prior_only
                    second.predata=base.predata; second.misfit=base.misfit; second.logL=0;
                elseif second.ptb.structure
                    chain.nforward=chain.nforward+1; chain.dr_forward_called(iter)=true;
                    try
                        second.predata=forward_fun(second.model,data,par,workdir);
                        second.misfit=b4_CALC_MISFIT(data,second.predata,par);
                        [second.logL,second.misfit]=b5_CALC_LIKELIHOOD(second.misfit,data,second.model.datahparm,par);
                        chain.consecutive_failures=0;
                    catch problem
                        context=struct('stage','delayed-rejection proposal','iteration',iter,'proposal',second.ptb, ...
                            'first_proposal',ptb,'accepted_model',state.model,'rejected_model',trial.model);
                        if ~strcmp(problem.identifier,'thbi:forward_failure')
                            preserve_forward_failure(fullfile(outdir,'fatal_forward_failure'),second.model,data,par,problem,workdir,context);
                            rethrow(problem);
                        end
                        preserve_forward_failure(fullfile(outdir,'first_forward_failure'),second.model,data,par,problem,workdir,context);
                        chain.failures=chain.failures+1; chain.forward_failed(iter)=true;
                        chain.consecutive_failures=chain.consecutive_failures+1; chain.last_failure=problem.message;
                        if chain.consecutive_failures>=par.inv.max_consecutive_failures
                            preserve_forward_failure(fullfile(outdir,'fatal_forward_failure'),second.model,data,par,problem,workdir,context);
                            rethrow(problem);
                        end
                    end
                else
                    second.predata=base.predata; second.misfit=base.misfit;
                    [second.logL,second.misfit]=b5_CALC_LIKELIHOOD(second.misfit,data,second.model.datahparm,par);
                end
            end
        end
        chain.dr_proposal_valid(iter)=second.ptb.valid && isfinite(second.logP);
        if layered_dr
            [accepted,logalpha2,dr_terms]=delayed_rejection_layered(state,trial,second,ptb,par);
            if isfield(dr_terms,'log_targeting_ratio')
                chain.dr_log_targeting_ratio(iter)=dr_terms.log_targeting_ratio;
                chain.dr_reverse_first_logalpha(iter)=dr_terms.logalpha_zy;
            end
        else
            [accepted,logalpha2]=delayed_rejection_damp(state,trial,second,ptb);
        end
        chain.dr_logalpha(iter)=logalpha2;
        % move/zone record the first leg; dr_move/dr_zone record the second.
        % proposal_valid means at least one candidate passed support.
        ptb.valid=ptb.valid || second.ptb.valid;
        if accepted
            chain.accepted_stage(iter)=2; trial=second; logalpha=logalpha2;
        end
    end
    if accepted, state=trial; end
    if model_fit_score(state.misfit.E2',par)<model_fit_score(chain.bestfit.misfit.E2',par)
        chain.bestfit=state; chain.bestfit_iteration=iter;
    end
    % Keep the exact response of the best RETAINED Earth, even when this
    % iteration is not on the response-saving stride. No extra forward call.
    if iter>par.inv.burnin && (isempty(chain.bestfit_retained) || ...
            model_fit_score(state.misfit.E2',par)<model_fit_score(chain.bestfit_retained.misfit.E2',par))
        chain.bestfit_retained=state; chain.bestfit_retained_iteration=iter;
    end
    chain=b9_SAVE_RESULT(chain,state,iter,ptb,accepted,logalpha,par);

    %% =================== PROGRESS AND RESTART STATE ===================
    if mod(iter,par.inv.print_every)==0 || iter==stop_iter
        fprintf('Chain %d | %d/%d | K [%d %d %d] | logL %.2f | accept %.2f | forward %d\n', ...
            ichain,iter,par.inv.niter,state.model.K,state.logL,mean(chain.accepted(1:iter)),chain.nforward);
    end
    if mod(iter,par.inv.checkpoint_every)==0 || iter==stop_iter
        chain.last_iter=iter;
        chain.elapsed_seconds=chain.elapsed_seconds+toc(clock0); clock0=tic;
        checkpoint=struct('spec',spec,'rng',rng,'iter',iter,'state',state,'chain',chain);
        save_checkpoint(filename,checkpoint);
    end
end
if stop_iter==par.inv.niter
    save(fullfile(outdir,'chain.mat'),'chain','-v7');
end
end
