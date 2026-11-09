function [qp, qs, qsh] = vti_vertical_slowness(matrix, p)
% Positive, propagating qP/qSV/qSH slownesses (s/km), flat VTI Christoffel.
% One ray parameter and any number of materials. qP is the smaller root;
% the linked weak-anisotropy THBI prior keeps the P/S branches separated.
e=vti_elastic(matrix);
assert(isscalar(p) && isreal(p) && isfinite(p) && p>=0,'fiveobs:rayp','Invalid ray parameter.');
% Density-normalized stiffness avoids avoidable units/scale cancellation.
a=e.A./e.rho; c=e.C./e.rho; l=e.L./e.rho; f=e.F./e.rho;
B=(a*p^2-1).*c+(l*p^2-1).*l-(f+l).^2*p^2;
C=(a*p^2-1).*(l*p^2-1); A=l.*c;
D=B.^2-4*A.*C;
assert(all(D>0 & B<0 & C>0),'fiveobs:rayp','P/SV must both propagate for RF migration.');
% Avoid subtraction of nearly equal numbers for the smaller root.
qs2=(-B+sqrt(D))./(2*A); qp2=C./(A.*qs2);
qsh2=(1-(e.N./e.rho)*p^2)./l;
assert(all(qp2>0 & qs2>qp2 & qsh2>0),'fiveobs:rayp','Nonpropagating or crossed VTI branches.');
qp=sqrt(qp2); qs=sqrt(qs2); qsh=sqrt(qsh2);
end
