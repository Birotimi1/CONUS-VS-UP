function G = reference_vti_compliance(matrix,frequency,slowness,damping,wave,method)
% Independent forced surface response from physical displacement/traction.
% Input: [h Vp_average Vs_average rho xi phi eta], km/km/s/g/cm^3;
% frequency Hz, horizontal slowness s/km. Output compliance is SI, with
% traction on a plane of positive downward normal (opposite surface load).
% G(:,:,ip,ifreq) is PSV [ux,uz] or scalar SH. Validation use only.
% No HVf eigenvectors, Wang normalization or project tensor helper is used.
if nargin<6, method='bounded_modes'; end
assert(any(strcmp(wave,{'PSV','SH'})) && any(strcmp(method,{'bounded_modes','expm'})));
frequency=frequency(:); slowness=slowness(:);
assert(size(matrix,2)==7 && matrix(end,1)==0 && all(matrix(1:end-1,1)>0));
assert(all(isfinite(matrix(:))) && all(frequency>0) && all(isfinite(frequency)) && ...
    all(slowness>=0) && all(isfinite(slowness)) && isscalar(damping) && damping>0 && isfinite(damping));

%% =================== DIMENSIONLESS ELASTIC EQUATIONS ===================
% exp(i*k*x-i*w*t), z positive downward; state [ux,uz,txz,tzz].
% txz=L*(ux'+i*k*uz), tzz=F*i*k*ux+C*uz'. Traction is scaled
% by rho_ref*v_ref*omega; depth by v_ref/omega. Positive Im(w)
% selects the causal, decaying subspace of each physical ODE matrix.
vp=matrix(:,2)*1000; vs=matrix(:,3)*1000; rho=matrix(:,4)*1000;
xi=matrix(:,5); phi=matrix(:,6); eta=matrix(:,7);
A=rho.*vp.^2.*5./(phi+4); C=A.*phi;
L=rho.*vs.^2.*3./(xi+2); N=L.*xi; F=eta.*(A-2*L);
assert(all(L>0 & N>0 & C>0 & A>N & (A-N).*C>F.^2),'thbi:reference','Unstable VTI reference tensor.');
vref=vs(end); rhoref=rho(end); scale=rhoref*vref^2;
A=A/scale; C=C/scale; F=F/scale; L=L/scale; N=N/scale; rho=rho/rhoref;
w=1+1i*damping; nl=size(matrix,1); nd=2; if strcmp(wave,'SH'), nd=1; end
G=complex(zeros(nd,nd,numel(slowness),numel(frequency)));

%% =================== GENERIC EIGENSYSTEM, NOT HVf POLARIZATION FORMULAS ===================
for ip=1:numel(slowness)
    k=slowness(ip)*vref/1000;
    modes=cell(nl,1); operators=cell(nl,1);
    for j=1:nl
        if nd==2
            M=[0 -1i*k 1/L(j) 0; -1i*k*F(j)/C(j) 0 0 1/C(j); ...
                (A(j)-F(j)^2/C(j))*k^2-rho(j)*w^2 0 0 -1i*k*F(j)/C(j); ...
                0 -rho(j)*w^2 -1i*k 0];
        else
            M=[0 1/L(j); N(j)*k^2-rho(j)*w^2 0];
        end
        [V,D]=eig(M); lambda=diag(D);
        down=find(real(lambda)<0); up=find(real(lambda)>0);
        assert(numel(down)==nd && numel(up)==nd,'thbi:reference','Causal subspace is not separated.');
        modes{j}=struct('Ud',V(1:nd,down),'Td',V(nd+1:end,down), ...
            'Uu',V(1:nd,up),'Tu',V(nd+1:end,up),'ld',lambda(down),'lu',lambda(up));
        operators{j}=M;
    end
    for it=1:numel(frequency)
        omega=2*pi*frequency(it); m=modes{nl}; Z=m.Td/m.Ud;
        for j=nl-1:-1:1
            m=modes{j}; thickness=matrix(j,1)*1000*omega/vref;
            if strcmp(method,'expm')
                % Independent small-stack control; never grow a thick/deep
                % transfer matrix where exponential cancellation is expected.
                assert(max(abs(real([m.ld;m.lu])))*thickness<12, ...
                    'thbi:reference','Use bounded modes for a thick evanescent layer.');
                P=expm(operators{j}*thickness);
                Z=(P(nd+1:end,nd+1:end)-Z*P(1:nd,nd+1:end))\ ...
                    (Z*P(1:nd,1:nd)-P(nd+1:end,1:nd));
            else
                % Reflect at the bottom impedance; transport both directions
                % using only decaying exponentials, then form top impedance.
                R=(m.Tu-Z*m.Uu)\(Z*m.Ud-m.Td);
                R=diag(exp(-m.lu*thickness))*R*diag(exp(m.ld*thickness));
                Z=(m.Td+m.Tu*R)/(m.Ud+m.Uu*R);
            end
        end
        G(:,:,ip,it)=(Z\eye(nd))/(rhoref*vref*omega);
    end
end
assert(all(isfinite(G(:))),'thbi:reference','Nonfinite reference compliance.');
end
