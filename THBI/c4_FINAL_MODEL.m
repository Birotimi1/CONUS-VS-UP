function final_model = c4_FINAL_MODEL(posterior, suite_of_models, retained, par)
% Pointwise posterior summaries plus a separate actual retained Earth model.
% A pointwise median across different layerings need not be a valid sampled
% Earth. Forward evaluation uses the explicitly identified retained model.

%% =================== DEPTH-DEPENDENT MATERIAL DISTRIBUTIONS ===================
assert(posterior.Nstored==size(suite_of_models.VS,2),'thbi:results','Inconsistent sample counts.');
final_model.Z=suite_of_models.Z;
final_model.quantile_probabilities=[.05 .5 .95];
for variable={'VS','VP','rho'}
    name=variable{1}; values=suite_of_models.(name);
    q=sample_quantiles(values',[.05 .5 .95])';
    final_model.([name 'quantiles'])=q;
    final_model.([name 'av'])=q(:,2); % legacy THBI field name: median, not mean
end
final_model.boundary_quantiles=sample_quantiles([posterior.zsed posterior.zmoh],[.05 .5 .95]);
final_model.gamma_quantiles=sample_quantiles(posterior.gamma,[.05 .5 .95]);
final_model.count_values=posterior.count_values;
final_model.count_probability=posterior.count_probability;

%% =================== REPRESENTATIVE PHYSICAL MODEL FROM RETAINED STATES ===================
assert(all(isfinite(posterior.base_energy)),'thbi:results','Invalid retained residual energies.');
[~,row]=min(posterior.base_energy);
final_model.representative_model=model_from_chain(retained,row,par);
assert(a1_TEST_CONDITIONS(final_model.representative_model,par),'thbi:results','Selected model violates prior support.');
final_model.representative_row=row;
final_model.representative_chain=posterior.chain(row);
final_model.representative_iteration=posterior.iteration(row);
final_model.selection='Minimum fixed-reference-noise fit score among ALL retained states; first occurrence breaks ties.';
final_model.note='Depth medians are summaries, not a sampled Earth or cross-dimension MAP; forward modelling uses representative_model.';
end
