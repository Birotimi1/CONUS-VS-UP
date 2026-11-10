function [selected,residual] = select_cps_fundamental(table,binary,matrix,periods,strict)
% Retain CPS for VTI phase/group, rejecting its branch-boundary artefacts.
% The lowest phase velocity satisfying the free-surface boundary condition
% is the physical fundamental. Raw CPS mode0 can be a repeated-root artefact.
if nargin<5, strict=true; end
selected=NaN(numel(periods),size(table,2)); residual=NaN(numel(periods),1);
for j=1:numel(periods)
    ix=find(abs(binary.period-periods(j))<1e-6*max(1,periods(j)));
    assert(isscalar(ix),'fiveobs:cps','CPS binary period mismatch.');
    c=binary.velocity{ix};
    rows=table(abs(table(:,3)-periods(j))<6e-5*max(1,periods(j)),:);
    [~,order]=sort(rows(:,1)); rows=rows(order,:);
    assert(size(rows,1)==numel(c) && isequal(rows(:,1),(0:numel(c)-1)'), ...
        'fiveobs:cps','CPS modal table/binary mismatch.');
    [~,order]=sort(c);
    for q=order(:)'
        if c(q)>=binary.halfspace_shear(1)*(1-1e-9), continue; end
        r=cps_rayleigh_boundary_check(matrix,periods(j),c(q));
        if r<1e-3
            selected(j,:)=rows(q,:); selected(j,5)=c(q); residual(j)=r;
            break;
        end
    end
    if strict
        assert(isfinite(residual(j)),'fiveobs:cps', ...
            'No CPS mode satisfies the Rayleigh surface condition at %.8g s.',periods(j));
    end
end
end
