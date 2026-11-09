function chain = b9_SAVE_RESULT(chain, state, iter, ptb, accepted, logalpha, par)
% Variable-size fractions/velocities are stored as cells, including repeats.
% Never align layer index j across models with different dimensions.

%% =================== PHYSICAL MODEL AND LAYER COUNTS ===================
fields={'sedmparm','crustmparm','mantmparm'};
for z=1:3
    chain.fractions{iter,z}=state.model.(fields{z}).f;
    chain.velocities{iter,z}=state.model.(fields{z}).VS;
end
chain.K(iter,:)=state.model.K;
chain.boundaries(iter,:)=[state.model.zsed state.model.zmoh];
chain.halfspace(iter,1)=state.model.halfspaceVS;
chain.eta(iter,:)=state.model.datahparm.eta';
if isfield(state.model,'vpvs'), chain.vpvs(iter,1)=state.model.vpvs; end

%% =================== TARGET AND PROPOSAL AUDIT TRAIL ===================
chain.logL(iter,1)=state.logL; chain.logP(iter,1)=state.logP;
chain.E2(iter,:)=state.misfit.E2';
chain.logL_indivdat(iter,:)=state.misfit.logL_indivdat';
chain.accepted(iter,1)=accepted;
chain.move(iter,1)=ptb.type_index; chain.zone(iter,1)=ptb.zone;
chain.logalpha(iter,1)=logalpha;
chain.log_q_ratio(iter,1)=ptb.log_q_reverse-ptb.log_q_forward;
chain.log_jacobian(iter,1)=ptb.log_abs_jacobian;
chain.proposal_valid(iter,1)=ptb.valid;

%% =================== SAMPLED RESPONSES ON THE FIXED DATA GRIDS ===================
if ~par.inv.prior_only && par.res.prediction_stride>0 && iter>par.inv.burnin && ...
        mod(iter-par.inv.burnin-1,par.res.prediction_stride)==0
    row=1+(iter-par.inv.burnin-1)/par.res.prediction_stride;
    assert(chain.prediction_iterations(row)==iter,'thbi:storage','Prediction iteration mismatch.');
    for d=1:numel(par.inv.datatypes)
        name=par.inv.datatypes{d};
        chain.predictions.(name)(row,:)=state.predata.(name)(:)';
    end
    chain.prediction_saved(row)=true;
end
end
