function dof = rf_effective_dof(daughter)
% Silver-Chan (1991) A12 spectral estimate, as in Eilon scdofcalc / dofTHBI.
% Adapted from eilonzach/seis_tools commit29260b271be6ed2d1d3e1a82250dbddc149417d2.
% Compute ONCE from the observed, cropped daughter, before any candidate.
x=daughter(:); assert(isreal(x) && all(isfinite(x)),'thbi:dof','Invalid RF daughter.');
nonzero=find(x~=0); assert(numel(nonzero)>=4,'thbi:dof','Insufficient nonzero RF samples.');
x=x(nonzero(1):nonzero(end)); power=abs(fft(x)).^2;
power=power(1:floor(numel(power)/2)); fourth=power.^2;
e2=sum(power)-.5*(power(1)+power(end));
e4=sum(fourth)-.5*(fourth(1)+fourth(end));
dof=2*(2*e2^2/e4-1);
assert(isfinite(dof) && dof>0,'thbi:dof','RF spectral DOF estimate is not positive.');
end
