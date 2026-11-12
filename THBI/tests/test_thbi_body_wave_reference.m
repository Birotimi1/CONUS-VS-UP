function report = test_thbi_body_wave_reference(par,folder,run_native)
% Independent VTI body-wave validation; no inversion or surface-mode search.
% Closed-form half-spaces -> physical expm -> bounded impedance -> native GBW.
if nargin<3, run_native=true; end
if ~exist(folder,'dir'), mkdir(folder); end
previous=pwd; restore=onCleanup(@() cd(previous)); %#ok<NASGU>
saved_rng=rng; restore_rng=onCleanup(@() rng(saved_rng)); %#ok<NASGU>
report=struct('passed',false,'native_executed',run_native,'started',datestr(now,30));
report.reference_sources=struct('compliance',fileread(which('reference_vti_compliance')), ...
    'integrals',fileread(which('reference_vti_body_integrals')),'test',fileread([mfilename('fullpath') '.m']));
frequency=[.1;1;25]; p=[0;.04;.12;.23;.4]; damping=.001;

%% =================== ANALYTIC ISOTROPIC HALF-SPACE PSV AND SH ===================
half=[0 6.2 3.5 2.7 1 1 1];
G=reference_vti_compliance(half,frequency,p,damping,'PSV');
S=reference_vti_compliance(half,frequency,p,damping,'SH');
analytic=G; analytic_sh=S;
rho=2700; vp=6200; vs=3500; mu=rho*vs^2;
for j=1:numel(frequency)
    w=2*pi*frequency(j)*(1+1i*damping);
    for k=1:numel(p)
        kh=2*pi*frequency(j)*p(k)/1000;
        qp=sqrt(w^2/vp^2-kh^2); qs=sqrt(w^2/vs^2-kh^2);
        U=[kh qs;qp -kh];
        T=1i*[2*mu*kh*qp rho*w^2-2*mu*kh^2;rho*w^2-2*mu*kh^2 -2*mu*kh*qs];
        analytic(:,:,k,j)=U/T;
        analytic_sh(:,:,k,j)=1/(1i*mu*qs);
    end
end
report.isotropic_psv_error=relative_norm(G,analytic);
report.isotropic_sh_error=relative_norm(S,analytic_sh);
assert(max([report.isotropic_psv_error report.isotropic_sh_error])<1e-10);

%% =================== VTI SH CLOSED FORM AND NORMAL-INCIDENCE PSV ===================
vti=[0 6.2 3.5 2.7 .8 1.2 .93];
A=rho*vp^2*5/(vti(6)+4); C=A*vti(6);
L=rho*vs^2*3/(vti(5)+2); N=L*vti(5);
S=reference_vti_compliance(vti,frequency,p,damping,'SH'); analytic_sh=S;
G0=reference_vti_compliance(vti,frequency,0,damping,'PSV'); zero=G0;
for j=1:numel(frequency)
    omega=2*pi*frequency(j); w=omega*(1+1i*damping);
    q=sqrt((rho*w^2-N*(omega*p/1000).^2)/L);
    analytic_sh(1,1,:,j)=1./(1i*L*q);
    zero(:,:,1,j)=diag(1./(1i*w*sqrt(rho*[L C])));
end
report.vti_sh_error=relative_norm(S,analytic_sh);
report.vti_normal_incidence_error=relative_norm(G0,zero);
assert(max([report.vti_sh_error report.vti_normal_incidence_error])<1e-10);

%% =================== INDEPENDENT EXPM AND BOUNDED-MODE LAYER RECURSIONS ===================
two=[.1 1.2 .5 1.8 1.2 .8 .9;0 6.2 3.5 2.7 1.15 .9 1.05];
report.expm_error=zeros(1,2);
waves={'PSV','SH'};
for k=1:2
    bounded=reference_vti_compliance(two,[.2;1;3],p,damping,waves{k});
    direct=reference_vti_compliance(two,[.2;1;3],p,damping,waves{k},'expm');
    report.expm_error(k)=relative_norm(bounded,direct);
    for j=1:numel(p)*3
        response=bounded(:,:,j); dissipative=-(response-response')/(2i);
        assert(min(real(eig(dissipative)))>=-1e-10*norm(response),'thbi:test','Non-passive reference response.');
    end
end
assert(max(report.expm_error)<1e-8);

%% =================== CLOSED-FORM INTEGRATED VTI SH HALF-SPACE ===================
integrated=reference_vti_body_integrals(vti,1,12000,damping,damping);
omega=2*pi; w=omega*(1+1i*damping); kmax=omega*sqrt(rho/N);
expected=imag(-(sqrt(N*kmax^2-rho*w^2)-sqrt(-rho*w^2))/(N*sqrt(L)))/(4*pi);
report.sh_integral_quadrature_error=abs(integrated(2)-expected)/abs(expected);
assert(report.sh_integral_quadrature_error<.003,'thbi:test','SH quadrature/normalization fails its analytic integral.');
fprintf('PASS independent compliance: isotropic/VTI half-space, expm, passivity and analytic SH integral.\n');

%% =================== NATIVE BODY CONTRIBUTIONS WITHOUT SURFACE-MODE CANCELLATION ===================
if run_native
    report.native_signature=check_native_vti_build(par.forc);
    truth=a2_BUILD_REFERENCE_TRUTH(par); base=truth.zoned_parameters.matrix;
    top=[0;cumsum(base(1:end-1,1))]; crust=top>=truth.zoned_parameters.zsed & top<truth.zoned_parameters.zmoh;
    low=base; low(crust,5)=.8; low(crust,6)=1.2;
    high=base; high(crust,5)=1.2; high(crust,6)=.8;
    lvz=base; inside=top>=77.5 & top<165; lvz(inside,2:3)=.90*lvz(inside,2:3);
    models={half,vti,two,base,low,high,lvz};
    names={'isotropic_halfspace','vti_halfspace','general_vti_layers', ...
        'ecm1_ak135_truth','crust_xi_min','crust_xi_max','mantle_lvz'};
    % Fixed independent prior draws exercise valid structures beyond the
    % truth. Restore the caller's RNG so validation cannot alter an experiment.
    rng(12112026,'twister');
    for j=1:4
        draw=b1_INITIATE_MODEL(par); assert(a1_TEST_CONDITIONS(draw,par));
        models{end+1}=draw.matrix; names{end+1}=sprintf('prior_draw_%02d',j); %#ok<AGROW>
    end
    report.cases=cell(numel(models),1);
    for k=1:numel(models)
        f=frequency; nks=512;
        if k==4, f=[.1;.2;.5;1;5;18;25]; nks=par.forc.hvsr.nks; end
        matrix=models{k}; work=fullfile(folder,names{k}); mkdir(work); cd(work);
        fprintf('Independent body-wave case %d/%d: %s, %d nodes, %d frequencies\n',k,numel(models),names{k},nks,numel(f));
        % Exercise the same numerical subdivision as the production wrapper.
        % The independent impedance reference keeps the physical layers: its
        % bounded exponentials do not require that numerical subdivision.
        native_matrix=refine_layer_model(matrix,par.forc.max_layer_thickness_km);
        exchange=native_matrix; exchange(:,1:4)=1000*exchange(:,1:4);
        fid=fopen('model.txt','w'); assert(fid>=0); fprintf(fid,'%d\n',size(native_matrix,1));
        fprintf(fid,[repmat('%.15g ',1,7) '\n'],exchange'); fclose(fid);
        fid=fopen('frequencies.txt','w'); assert(fid>=0); fprintf(fid,'%.15g\n',f); fclose(fid);
        command=sprintf('%s -vti -f model.txt -ff frequencies.txt -nmr 0 -nml 0 -nks %d -apsv %.15g -ash %.15g -rep -hv', ...
            quote_shell_path(fullfile(par.forc.hvf_vti_bin,'HVf_VTI')),nks,par.forc.hvsr.apsv,par.forc.hvsr.ash);
        [status,message]=run_native_command([command ' > HV.dat 2> HVf.log'],par.forc,['body reference/' names{k}]);
        assert(status==0,'thbi:test','Body-only native call failed: %s',message);
        native=dlmread('GBW.dat','',1,0);
        assert(isequal(size(native),[numel(f) 4]) && max(abs(native(:,1)-f))<1e-8);
        reference=reference_vti_body_integrals(matrix,f,nks,par.forc.hvsr.apsv,par.forc.hvsr.ash);
        error=abs(native(:,2:4)-reference)./abs(reference);
        result=struct('name',names{k},'matrix',matrix,'native_matrix',native_matrix,'frequency_Hz',f,'nks',nks, ...
            'native',native(:,2:4),'reference',reference,'relative_error',error);
        report.cases{k}=result; save(fullfile(folder,'body_reference_report.mat'),'report','-v7');
        fprintf('  maximum component relative error %.6g\n',max(error(:)));
        assert(all(error(:)<1e-5),'thbi:test','Independent anisotropic body-wave disagreement in %s.',names{k});
    end
end
report.passed=true; report.completed=datestr(now,30);
save(fullfile(folder,'body_reference_report.mat'),'report','-v7');
fprintf('PASS independent VTI body-wave validation; native comparison executed=%d.\n',run_native);
end

function e=relative_norm(a,b)
e=norm(a(:)-b(:))/norm(b(:));
end
