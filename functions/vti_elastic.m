function e = vti_elastic(matrix)
% THBI velocity convention -> Love elastic constants for vertical TI.
% Rows: [h, Vp, Vs, rho, xi, phi, eta], or four isotropic columns.
% xi=N/L=(Vsh/Vsv)^2, phi=C/A=(Vpv/Vph)^2, eta=F/(A-2L).
% Vs^2=(2 Vsv^2+Vsh^2)/3; Vp^2=(Vpv^2+4 Vph^2)/5.
% These are the reference THBI parameter averages, not a general TI Voigt
% isotropic projection. Constants have units GPa for km/s and g/cm^3.
assert(isnumeric(matrix) && isreal(matrix) && ~isempty(matrix) && ...
    any(size(matrix,2)==[4 7]) && all(isfinite(matrix(:))), ...
    'fiveobs:model','Expected four isotropic or seven VTI columns.');
n=size(matrix,1); xi=ones(n,1); phi=xi; eta=xi;
if size(matrix,2)==7, xi=matrix(:,5); phi=matrix(:,6); eta=matrix(:,7); end
assert(all(xi>0 & phi>0 & matrix(:,2)>0 & matrix(:,3)>0 & matrix(:,4)>0), ...
    'fiveobs:model','Nonpositive elastic parameters.');
e.Vsv=matrix(:,3).*sqrt(3./(xi+2)); e.Vsh=e.Vsv.*sqrt(xi);
e.Vph=matrix(:,2).*sqrt(5./(phi+4)); e.Vpv=e.Vph.*sqrt(phi);
e.rho=matrix(:,4); e.A=e.rho.*e.Vph.^2; e.C=e.rho.*e.Vpv.^2;
e.L=e.rho.*e.Vsv.^2; e.N=e.rho.*e.Vsh.^2; e.F=eta.*(e.A-2*e.L);
assert(all(e.C>0 & e.L>0 & e.N>0 & e.A>e.N & ...
    (e.A-e.N).*e.C>e.F.^2),'fiveobs:model','TI stiffness is not positive definite.');
e.xi=xi; e.phi=phi; e.eta=eta;
end
