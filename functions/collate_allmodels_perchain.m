function retained = collate_allmodels_perchain(chains, par)
% Collate every post-burn-in state, including rejected repeats and one chain.
% Adapted from Zach/Brennan's helper to the layered column/cell history format.
% No best-fit/accepted-only filtering, chain ranking, or implicit thinning.
assert(iscell(chains) && ~isempty(chains),'thbi:results','Supply a nonempty chain cell array.');
keep=(par.inv.burnin+1):par.inv.niter;
assert(~isempty(keep),'thbi:results','No retained iterations.');
n=numel(keep); nc=numel(chains);
fields={'fractions','velocities','K','boundaries','halfspace','eta','logL','E2'};
if isfield(par.mod,'vpvs_bounds'), fields{end+1}='vpvs'; end
if isfield(par.mod,'xi_bounds'), fields{end+1}='xi'; end
retained=struct();
for k=1:numel(fields)
    name=fields{k}; pieces=cell(nc,1);
    for c=1:nc
        chain=chains{c};
        assert(chain.last_iter==par.inv.niter,'thbi:results','Chain %d is incomplete.',c);
        assert(size(chain.(name),1)>=par.inv.niter,'thbi:results','Missing history: %s.',name);
        pieces{c}=chain.(name)(keep,:);
    end
    retained.(name)=vertcat(pieces{:});
end
retained.chain=kron((1:nc)',ones(n,1));
retained.iteration=repmat(keep(:),nc,1);
retained.Nstored=n*nc;
retained.note='All retained states, with repeated states carrying their full iteration weight.';
end
