function chain = run_one_chain(par, data, ichain, stop_iter, forward_fun)
% Fixed-dimensional HBI chain. Restore RNG, adaptation, cached predictions,
% and repeated states exactly on resume. forward_fun is injectable for tests.
if nargin < 4 || isempty(stop_iter), stop_iter = par.inv.niter; end
if nargin < 5, forward_fun = @b3_FORWARD_MODEL; end
stop_iter = min(stop_iter,par.inv.niter);
outdir = fullfile(par.res.dir,sprintf('chain_%02d',ichain));
if ~exist(outdir,'dir'), mkdir(outdir); end
filename = fullfile(outdir,'checkpoint.mat');
spec = struct('schema',1,'par',par,'data',data,'forward',func2str(forward_fun));
workdir = tempname(par.res.scratch_root); mkdir(workdir);
cleanup = onCleanup(@() remove_scratch(workdir)); %#ok<NASGU>
clock0 = tic;
if exist(filename,'file')
    loaded = load(filename,'checkpoint'); cp = loaded.checkpoint;
    if ~isequaln(cp.spec,spec)
        error('hbi:resume','Checkpoint configuration/data changed. Use a new run directory.');
    end
    rng(cp.rng); state=cp.state; chain=cp.chain; tuning=cp.tuning;
    first = cp.iter+1;
else
    rng(par.inv.seed+1009*ichain,'twister');
    state.model = b1_INITIATE_MODEL(par);
    state.predata = forward_fun(state.model,data,par,workdir);
    state.misfit = b4_CALC_MISFIT(data,state.predata,par);
    [state.logL,state.misfit] = b5_CALC_LIKELIHOOD(state.misfit,data,state.model.datahparm,par);
    state.logP = a3_CALC_LOG_PRIOR(state.model,par);
    if ~isfinite(state.logL+state.logP), error('hbi:start','Nonfinite initial posterior.'); end
    nt = numel(state.model.theta); nd = numel(par.inv.datatypes); n = par.inv.niter;
    chain = struct('theta',NaN(n,nt),'eta',NaN(n,nd),'logL',NaN(n,1), ...
        'logP',NaN(n,1),'E2',NaN(n,nd),'logL_indivdat',NaN(n,nd), ...
        'accepted',false(n,1),'move',zeros(n,1),'initial',state, ...
        'bestfit',state,'bestposterior',state,'nforward',1,'failures',0, ...
        'consecutive_failures',0,'last_failure','','elapsed_seconds',0);
    tuning.base = [par.mod.step;par.mod.noise.step];
    tuning.steps = tuning.base; tuning.attempts = zeros(nt+nd,1);
    first = 1;
end
for iter = first:stop_iter
    [trial.model,move] = b2_PERTURB_MODEL(state.model,par,tuning.steps);
    trial.logP = a3_CALC_LOG_PRIOR(trial.model,par); trial.logL = -Inf;
    if isfinite(trial.logP)
        if move.structure
            chain.nforward = chain.nforward+1;
            try
                trial.predata = forward_fun(trial.model,data,par,workdir);
                trial.misfit = b4_CALC_MISFIT(data,trial.predata,par);
                chain.consecutive_failures = 0;
            catch problem
                % Programming/configuration errors must stop, not masquerade
                % as rejected models. Known solver failures are counted.
                if ~any(strcmp(problem.identifier,{'hbi:forward','hbi:prediction'})) && ...
                        ~strncmp(problem.identifier,'fiveobs:',8)
                    rethrow(problem);
                end
                chain.failures = chain.failures+1;
                chain.consecutive_failures = chain.consecutive_failures+1;
                chain.last_failure = problem.message;
                if chain.consecutive_failures >= par.inv.max_consecutive_failures
                    copyfile(workdir,fullfile(outdir,'failed_forward'));
                    rethrow(problem);
                end
                trial.misfit = [];
            end
        else
            trial.predata = state.predata;
            trial.misfit = state.misfit; % same structural residuals and C0
        end
        if ~isempty(trial.misfit)
            [trial.logL,trial.misfit] = b5_CALC_LIKELIHOOD(trial.misfit,data,trial.model.datahparm,par);
        end
    end
    accepted = b6_IFACCEPT(state.logL+state.logP,trial.logL+trial.logP);
    if accepted, state = trial; end
    if sum(state.misfit.E2) < sum(chain.bestfit.misfit.E2), chain.bestfit=state; end
    if state.logL+state.logP > chain.bestposterior.logL+chain.bestposterior.logP
        chain.bestposterior=state;
    end
    chain = b9_SAVE_RESULT(chain,state,iter,move,accepted);
    if par.inv.adapt && iter <= par.inv.burnin
        j = move.index; tuning.attempts(j) = tuning.attempts(j)+1;
        gain = (tuning.attempts(j)+10)^(-.6);
        tuning.steps(j) = tuning.steps(j)*exp(gain*(double(accepted)-.44));
        tuning.steps(j) = min(10*tuning.base(j),max(.05*tuning.base(j),tuning.steps(j)));
    end
    if mod(iter,par.inv.print_every)==0 || iter==stop_iter
        fprintf('Chain %d | %d/%d | logL %.2f | base chi2 %.1f | accept %.2f | forward %d\n', ...
            ichain,iter,par.inv.niter,state.logL,sum(state.misfit.E2), ...
            mean(chain.accepted(1:iter)),chain.nforward);
    end
    if mod(iter,par.inv.checkpoint_every)==0 || iter==stop_iter
        chain.last_iter = iter; chain.final_steps=tuning.steps;
        chain.elapsed_seconds = chain.elapsed_seconds+toc(clock0); clock0=tic;
        checkpoint = struct('spec',spec,'rng',rng,'iter',iter,'state',state, ...
            'chain',chain,'tuning',tuning);
        save_checkpoint(filename,checkpoint);
    end
end
if stop_iter == par.inv.niter
    save(fullfile(outdir,'chain.mat'),'chain','-v7');
end
end
