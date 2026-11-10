function report = test_thbi_vti_native(par,folder,anis)
% Native validation gate, no MCMC. Keep every model/SAC/log.
% Compare VTI xi=phi=eta=1 with original isotropic operators, then compare
% anisotropic HVf fundamental phase/group/ellipticity with independent CPS.
assert(strcmp(par.mod.anisotropy,'vti'),'thbi:preflight','VTI profile required.');
signature=check_native_vti_build(par.forc);
if ~exist(folder,'dir'), mkdir(folder); end
[~,paths]=a0_STARTUP_BAYES();
iso_par=bayes_inv_parms('reference100'); iso_par.forc.cps_bin=par.forc.cps_bin;
iso_par.forc.hvf_bin=paths.hvf_bin; iso_par.forc.verbose=false;
par.forc.verbose=false;
T=[.5;1;5;10;30;60;100];
data.hvsr.x=[.1;.5;1;5;20];
for name={'ellip','cph','ugr'}, data.(name{1}).x=T; end
if nargin<3
    truth=a2_BUILD_REFERENCE_TRUTH(par); anis=truth.zoned_parameters;
end
assert(a1_TEST_CONDITIONS(anis,par),'thbi:preflight','Native fixture is outside prior support.');
iso_vti=anis; iso_vti.xi=ones(1,3); iso_vti=make_mod_from_parms(iso_vti,par);
iso=iso_vti; iso.matrix=iso.matrix(:,1:4);
% Keeping data.rf absent requests native Ps waves, not an HK scalar.
[pi,fi]=b3_FORWARD_MODEL(iso,data,iso_par,fullfile(folder,'isotropic'));
[pv,fv]=b3_FORWARD_MODEL(iso_vti,data,par,fullfile(folder,'vti_isotropic_limit'));
[pa,fa]=b3_FORWARD_MODEL(anis,data,par,fullfile(folder,'vti_nonzero'));
names=par.inv.datatypes;
report=struct('passed',false,'isotropic_limit',struct(),'anisotropic_change',struct(), ...
    'native_signature',signature);
for j=1:numel(names)
    name=names{j};
    report.isotropic_limit.(name)=norm(pv.(name)-pi.(name))/max(norm(pi.(name)),1e-12);
    report.anisotropic_change.(name)=norm(pa.(name)-pv.(name))/max(norm(pv.(name)),1e-12);
    if j<=4
        % A large HVSR resonance must not hide a mismatch at other frequencies.
        report.isotropic_max_pointwise.(name)=max(abs(pv.(name)-pi.(name))./abs(pi.(name)));
    end
end
save(fullfile(folder,'comparisons.mat'),'report','pi','pv','pa','fi','fv','fa','-v7');
% These are numerical agreement gates, not posterior recovery criteria.
limits=[.01 .01 .003 .015 .02 .02];
for j=1:numel(names)
    assert(report.isotropic_limit.(names{j})<limits(j),'thbi:preflight', ...
        'VTI isotropic-limit mismatch for %s: %.4g. Inspect %s.', ...
        names{j},report.isotropic_limit.(names{j}),folder);
    assert(report.anisotropic_change.(names{j})>1e-7,'thbi:preflight', ...
        '%s did not respond to anisotropy.',names{j});
    if j<=4
        assert(report.isotropic_max_pointwise.(names{j})<limits(j),'thbi:preflight', ...
            'VTI isotropic-limit pointwise mismatch for %s: %.4g.', ...
            names{j},report.isotropic_max_pointwise.(names{j}));
    end
end

%% =================== INDEPENDENT CPS/HVf NONZERO-ANISOTROPY PARITY ===================
% CPS tdpegn96 -ASC columns: mode,index,T,f,C,U,energy,gamma,ellipticity.
cps=dlmread(fullfile(fa.output_dir,'cph','TREGN.ASC'),'',1,0);
assert(size(cps,2)==9 && all(cps(:,1)==0),'thbi:preflight','Unexpected CPS TI ASCII schema.');
[~,order]=sort(cps(:,3)); cps=cps(order,:);
assert(numel(T)==size(cps,1) && max(abs(cps(:,3)-T))<1e-3);
previous=pwd; cleanup=onCleanup(@() cd(previous)); %#ok<NASGU>
cd(fullfile(fa.output_dir,'ellip'));
cmd=sprintf('%s -vti -f model.txt -ff frequencies.txt -nmr 1 -nml 0 -nks 0 -prec %.10g -ph -gr -hv', ...
    quote_shell_path(fullfile(par.forc.hvf_vti_bin,'HVf_VTI')),fa.config.hvsr.prec);
[status,message]=system([cmd ' > parity_HV.dat 2> parity_HVf.log']);
assert(status==0,'thbi:preflight','HVf parity run failed: %s',message);
phase=read_modal_slowness('Rph.dat',numel(T));
group=read_modal_slowness('Rgr.dat',numel(T));
% HVf stores frequency ascending, hence requested period descending.
phase=flipud(1./phase)/1000; group=flipud(1./group)/1000;
report.cps_hvf_phase=norm(phase-cps(:,5))/norm(cps(:,5));
report.cps_hvf_group=norm(group-cps(:,6))/norm(cps(:,6));
report.cps_hvf_ellipticity=norm(pa.ellip-abs(cps(:,9)))/norm(cps(:,9));
cd(previous);
save(fullfile(folder,'comparisons.mat'),'report','pi','pv','pa','fi','fv','fa','cps','phase','group','-v7');
assert(report.cps_hvf_phase<.003 && report.cps_hvf_group<.015 && ...
    report.cps_hvf_ellipticity<.02,'thbi:preflight', ...
    'Anisotropic CPS/HVf parity failed. Inspect comparisons.mat before any inversion.');
report.passed=true;
report.note='Native numerical comparisons passed; not evidence of posterior convergence or 250-km resolution.';
save(fullfile(folder,'report.mat'),'report','-v7');
fprintf('PASS native VTI gate: phase %.3g, group %.3g, ellipticity %.3g relative differences.\n', ...
    report.cps_hvf_phase,report.cps_hvf_group,report.cps_hvf_ellipticity);
end

function slow=read_modal_slowness(filename,n)
text=fileread(filename); token=regexp(strtrim(text),'\s+','split');
assert(numel(token)==2+2*n && str2double(token{1})==n && str2double(token{2})==1, ...
    'thbi:preflight','Incomplete single-mode output: %s',filename);
slow=str2double(token(3:2+n))'; flags=token(3+n:end);
assert(all(isfinite(slow) & slow>0) && all(strcmpi(flags,'T')), ...
    'thbi:preflight','Missing/nonphysical mode in %s',filename);
end
