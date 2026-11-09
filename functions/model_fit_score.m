function score = model_fit_score(E2,par)
% One display-model ranking, fixed before sampling; NOT a posterior density.
% Input rows are states and columns are data blocks. Legacy profiles use sum E2.
assert(size(E2,2)==numel(par.inv.datatypes),'thbi:results','Fit score shape mismatch.');
weights=ones(numel(par.inv.datatypes),1);
if isfield(par,'rf') && strcmp(par.rf.method,'hk_ccp')
    weights(5:6)=exp(-2*par.mod.noise.log_center(5:6));
end
score=E2*weights;
end
