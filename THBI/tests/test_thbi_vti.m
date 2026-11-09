function test_thbi_vti()
% Solver-free checks: independent Christoffel eigenvalues, RF timing and IO.
% No external executable, wavefield forward calculation or inversion.
par=bayes_inv_parms('reference100_vti'); truth=a2_BUILD_REFERENCE_TRUTH(par);
m=truth.zoned_parameters; e=vti_elastic(m.matrix);
assert(isequal(m.xi,[1 1.06 1]) && size(m.matrix,2)==7);
assert(max(abs((2*e.Vsv.^2+e.Vsh.^2)/3-m.VS.^2))<1e-12);
assert(max(abs((e.Vpv.^2+4*e.Vph.^2)/5-m.VP.^2))<1e-12);
assert(max(abs(e.phi-(2-e.xi)))<1e-12);

%% =================== INDEPENDENT CHRISTOFFEL DETERMINANT ===================
for p=[0 .04 .075]
    [qp,qs,qsh]=vti_vertical_slowness(m.matrix,p);
    for j=1:numel(qp)
        for q=[qp(j),qs(j)]
            G=[e.A(j)*p^2+e.L(j)*q^2,(e.F(j)+e.L(j))*p*q; ...
                (e.F(j)+e.L(j))*p*q,e.L(j)*p^2+e.C(j)*q^2]/e.rho(j);
            assert(min(abs(eig(G)-1))<1e-12);
        end
        assert(abs((e.N(j)*p^2+e.L(j)*qsh(j)^2)/e.rho(j)-1)<1e-12);
    end
end
iso=m.matrix; iso(:,5:7)=1;
[qp,qs]=vti_vertical_slowness(iso,.06);
assert(max(abs(qp-sqrt(iso(:,2).^-2-.06^2)))<1e-12);
assert(max(abs(qs-sqrt(iso(:,3).^-2-.06^2)))<1e-12);
bad=m; bad.xi(1)=1.01; assert(~a1_TEST_CONDITIONS(bad,par));
bad=m; bad.xi(2)=1.21; assert(~a1_TEST_CONDITIONS(bad,par));
% Legacy isotropic entry points must reject VTI BEFORE any native call or IO.
legacy={@() forward_five_observables(m.matrix,struct()), ...
    @() hv_prediction(m.matrix,1,'hvsr',struct(),''), ...
    @() dispersion_prediction(m.matrix,1,'cph',struct(),'')};
for j=1:numel(legacy)
    rejected=false;
    try, legacy{j}(); catch problem, rejected=strcmp(problem.identifier,'fiveobs:anisotropy'); end
    assert(rejected,'thbi:test','A legacy isotropic entry point accepted VTI.');
end

%% =================== NUMERICAL SUBDIVISION RETAINS ANISOTROPY INTERFACES ===================
fixture=[1 6.3 3.6 2.7 1 1 1;1 6.3 3.6 2.7 1.06 .94 1;0 8 4.5 3.3 1 1 1];
refined=refine_layer_model(fixture,.4);
assert(size(refined,1)==7 && sum(refined(1:3,1))==1 && all(refined(4:6,5)==1.06));
z=[-.5;0;.5;1;1.5;2;3];
lag=layered_conversion_time(fixture,z,.06,'Sp');
assert(max(abs(lag-layered_conversion_time(refined,z,.06,'Sp')))<1e-12);
[qp,qs]=vti_vertical_slowness(fixture,.06);
assert(abs(lag(end)+sum(qs-qp))<1e-12);

%% =================== CPS MODEL96 READ-BACK ===================
filename=[tempname '.mod']; cleanup=onCleanup(@() delete(filename)); %#ok<NASGU>
write_elastic_cps_model(filename,fixture); f=fopen(filename,'r');
for j=1:13, line=fgetl(f); if j==3, assert(strcmp(line,'TRANSVERSE ISOTROPIC')); end, end
ee=vti_elastic(fixture);
for j=1:3
    first=sscanf(fgetl(f),'%f'); second=sscanf(fgetl(f),'%f');
    assert(numel(first)==10 && numel(second)==3);
    assert(max(abs(first(1:4)-[fixture(j,1);ee.Vpv(j);ee.Vsv(j);ee.rho(j)]))<1e-9);
    assert(abs(second(3)^2*ee.rho(j)-ee.F(j))<1e-8);
end
fclose(f);

%% =================== VTI H-KAPPA ARRIVALS FROM INDEPENDENT POLYNOMIAL ROOTS ===================
H=m.zmoh; h=m.matrix(1:4,1); w=h/sum(h);
vs=1/sum(w./m.VS(1:4)); xi=sum(w.*e.xi(1:4)); phi=sum(w.*e.phi(1:4));
vp=vs*m.vpvs; a=vp^2*5/(phi+4); c=a*phi; l=vs^2*3/(xi+2); f=a-2*l;
t=(-100:.2:500)'; ps=par.rf.hk.rayp_s_km; rf=zeros(numel(t),numel(ps));
for j=1:numel(ps)
    p=ps(j); roots_q2=sort(roots([l*c,(a*p^2-1)*c+(l*p^2-1)*l-(f+l)^2*p^2, ...
        (a*p^2-1)*(l*p^2-1)]));
    q=sqrt(roots_q2); times=H*[q(2)-q(1),q(2)+q(1),2*q(2)]; amp=[1 .6 -.4];
    for k=1:3, rf(:,j)=rf(:,j)+amp(k)*exp(-.5*((t-times(k))/.35).^2); end
end
observed.waves=struct('t',t,'rf',rf,'rayp_s_km',ps);
hk=b3_FORWARD_MODEL_HK(m,observed,par); assert(hk.E_by_Emax>.99);
wrong=m; wrong.xi(2)=.88; wrong=make_mod_from_parms(wrong,par);
bad_hk=b3_FORWARD_MODEL_HK(wrong,observed,par);
assert(bad_hk.E_by_Emax<hk.E_by_Emax-.1);

%% =================== OBSERVATION-ONLY SPECTRAL DEGREES OF FREEDOM ===================
x=sin((0:1023)'*.2).*exp(-((0:1023)'-512).^2/10000);
dof=rf_effective_dof(x);
assert(isfinite(dof) && dof>0 && abs(dof-rf_effective_dof(7*x))<1e-10);
fprintf('PASS: VTI tensor, Christoffel roots, isotropic limit, HK arrivals, migration, model96 IO and RF DOF.\n');
end
