function model = model_from_chain(chain, iter, par)
% Reconstruct one saved variable-dimensional state, never an averaged layer list.
model.zsed=chain.boundaries(iter,1); model.zmoh=chain.boundaries(iter,2);
model.halfspaceVS=chain.halfspace(iter);
if isfield(chain,'vpvs'), model.vpvs=chain.vpvs(iter); end
if isfield(chain,'xi'), model.xi=chain.xi(iter,:); end
model.datahparm.eta=chain.eta(iter,:)';
for z=1:3
    model.(par.mod.zone_fields{z})=struct('f',chain.fractions{iter,z},'VS',chain.velocities{iter,z});
end
model=make_mod_from_parms(model,par);
end
