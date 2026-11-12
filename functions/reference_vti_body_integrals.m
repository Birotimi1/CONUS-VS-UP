function body = reference_vti_body_integrals(matrix,frequency,nks,apsv,ash)
% Independent surface impedance integrated on HVf's documented contour mesh.
% Output columns match GBW.dat: horizontal PSV, horizontal SH, vertical PSV.
% Comparing the same quadrature isolates elastic physics from mesh error.
assert(isscalar(nks) && isfinite(nks) && nks>=1 && nks==fix(nks));
frequency=frequency(:); xi=matrix(end,5); vs=matrix(end,3);
vsv=vs*sqrt(3/(xi+2)); vsh=vsv*sqrt(xi);
t=(0:nks)'; t(1)=.1; % original rectangle endpoint convention, not k=0
ps=t/(vsv*nks); ss=t/(vsh*nks); % s/km; integrals below use SI slowness
Gp=reference_vti_compliance(matrix,frequency,ps,apsv,'PSV');
Gs=reference_vti_compliance(matrix,frequency,ss,ash,'SH');
body=zeros(numel(frequency),3);
for j=1:numel(frequency)
    omega=2*pi*frequency(j);
    body(j,1)=omega^2*sum(ps.*imag(reshape(Gp(1,1,:,j),[],1)))/(4*pi*vsv*nks*1e6);
    body(j,2)=omega^2*sum(ss.*imag(reshape(Gs(1,1,:,j),[],1)))/(4*pi*vsh*nks*1e6);
    body(j,3)=omega^2*sum(ps.*imag(reshape(Gp(2,2,:,j),[],1)))/(2*pi*vsv*nks*1e6);
end
assert(all(isfinite(body(:))) && all(body(:)<0),'thbi:reference','Reference spectral sign/finite check failed.');
end
