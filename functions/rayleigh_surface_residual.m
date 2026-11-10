function residual = rayleigh_surface_residual(matrix,period,velocity)
% Independent free-surface boundary test using expm and QR, without layer
% P/SV eigenvector division. This remains regular at a finite-layer repeated
% vertical eigenvalue where CPS TI dispersion can emit spurious roots.
% State [u,v,tau_xz,tau_zz] uses the same real transformed convention as HVf.
e=vti_elastic(matrix); omega=2*pi/period; k=omega/velocity;
scale=e.rho(end)*velocity^2*k;
S=diag([1 1 1/scale 1/scale]);
generator=@(j) S*[0 k 1/e.L(j) 0; -e.F(j)*k/e.C(j) 0 0 1/e.C(j); ...
    (e.A(j)-e.F(j)^2/e.C(j))*k^2-e.rho(j)*omega^2 0 0 e.F(j)*k/e.C(j); ...
    0 -e.rho(j)*omega^2 -k 0]/S;
M=generator(size(matrix,1)); [V,D]=eig(M); lambda=diag(D);
keep=real(lambda)<-1e-12*max(abs(lambda));
if sum(keep)~=2, residual=Inf; return; end % radiating half-space, not trapped
Y=V(:,keep); [Y,~]=qr(Y,0);
for j=size(matrix,1)-1:-1:1
    M=generator(j);
    % Limit exponential growth per QR step to prevent column collapse.
    n=max(1,ceil(matrix(j,1)*max(abs(real(eig(M))))/2));
    step=expm(-matrix(j,1)*M/n);
    for q=1:n, [Y,~]=qr(step*Y,0); end
end
s=svd(Y(3:4,:)); residual=min(s);
end
