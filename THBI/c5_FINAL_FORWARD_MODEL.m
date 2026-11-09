function final_forward = c5_FINAL_FORWARD_MODEL(final_model, par, data, forward_fun)
% Explicitly recompute all five observables of the selected retained Earth.
% This stage invokes CPS/HVf when the user enables it in the DATA script.
if nargin<4, forward_fun=@b3_FORWARD_MODEL; end
assert(~par.inv.prior_only,'thbi:forward','Final data prediction requires likelihood data.');

%% =================== VALIDATE THE ACTUAL RETAINED MODEL ===================
model=final_model.representative_model;
assert(a1_TEST_CONDITIONS(model,par),'thbi:forward','Final model is outside prior support.');
workdir=fullfile(par.res.dir,'final_forward');
if ~exist(workdir,'dir'), mkdir(workdir); end

%% =================== FIVE RESPONSES, MISFIT AND HIERARCHICAL LIKELIHOOD ===================
final_forward.model=model;
final_forward.selection=final_model.selection;
final_forward.chain=final_model.representative_chain;
final_forward.iteration=final_model.representative_iteration;
if nargout(forward_fun)>1
    [final_forward.predata,final_forward.native]=forward_fun(model,data,par,workdir);
else
    final_forward.predata=forward_fun(model,data,par,workdir);
end
final_forward.misfit=b4_CALC_MISFIT(data,final_forward.predata,par);
[final_forward.logL,final_forward.misfit]=b5_CALC_LIKELIHOOD( ...
    final_forward.misfit,data,model.datahparm,par);
final_forward.logP=a3_CALC_LOG_PRIOR(model,par);
save(fullfile(workdir,'final_forward.mat'),'final_forward','-v7');
end
