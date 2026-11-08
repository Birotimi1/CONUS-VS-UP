function chain = b9_SAVE_RESULT(chain, state, iter, move, accepted)
% Save the CURRENT state at EVERY step, including rejected/repeated states.

%% Structural parameters and hierarchical noise
chain.theta(iter,:) = state.model.theta';
chain.eta(iter,:) = state.model.datahparm.eta';

%% Posterior terms and per-observable diagnostics
chain.logL(iter,1) = state.logL;
chain.logP(iter,1) = state.logP;
chain.E2(iter,:) = state.misfit.E2';
chain.logL_indivdat(iter,:) = state.misfit.logL_indivdat';
%% Proposal bookkeeping; rejected steps retain the previous state's values
chain.accepted(iter,1) = accepted;
chain.move(iter,1) = move.index;
end
