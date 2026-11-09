function [data,TRUEmodel,forward] = a2_LOAD_REFERENCE_DATA(par,truth_only)
% Actual synthetic loading blocks for four curves + Ps H-kappa + Sp CCP.
% Observed RF products are constructed ONCE before MCMC, never per proposal.

%% =================== EXTRACTED KNOWN STRUCTURE ===================
TRUEmodel=a2_BUILD_REFERENCE_TRUTH(par);
data=struct(); forward=struct('status','not_computed');
if truth_only, return; end
model=TRUEmodel.zoned_parameters;
scratch=tempname(par.res.scratch_root); mkdir(scratch);
cleanup=onCleanup(@() rmdir(scratch,'s')); %#ok<NASGU>
old_rng=rng; restore_rng=onCleanup(@() rng(old_rng)); %#ok<NASGU>
rng(par.synth.seed,'twister');

%% =================== FOUR SURFACE CURVES AND NATIVE Ps/Sp ===================
coords={'f','T','T','T'};
for k=1:4
    name=par.inv.datatypes{k}; data.(name).x=par.synth.grids.(name).(coords{k});
end
try
    [~,forward]=b3_FORWARD_MODEL(model,data,par,scratch);
catch problem
    preserve_forward_failure(fullfile(par.res.dir,'truth_forward_failure'), ...
        model,data,par,problem,scratch,struct('stage','reference synthetic truth','iteration',0));
    rethrow(problem);
end
for k=1:4
    name=par.inv.datatypes{k}; v=data.(name); v.clean=forward.pred.(name).val(:);
    v.N=numel(v.x); v.sig=max(par.synth.sigma_floor(k),par.synth.relative_sigma(k)*abs(v.clean));
    v.R=eye(v.N); v.L0=diag(v.sig); v.logdetC0=2*sum(log(v.sig));
    v.obs=v.clean+par.synth.gamma(k)*v.L0*randn(v.N,1); data.(name)=v;
end

%% =================== OBSERVED Ps WAVES FOR H-KAPPA STACKING ===================
rayp=par.rf.hk.rayp_s_km(:); waves=struct();
for j=1:numel(rayp)
    cfg=par.forc.rf; cfg.rayp_s_km=rayp(j);
    if abs(rayp(j)-par.forc.rf.rayp_s_km)<1e-12
        trace=forward.pred.rf;
    else
        trace=forward_receiver_function_cps(forward.model,cfg,forward.config.cps_bin, ...
            fullfile(scratch,sprintf('Ps_event_%02d',j)));
        trace=receiver_function_bands(trace,cfg);
    end
    ix=find(trace.t>=par.synth.rf_window(1) & trace.t<=par.synth.rf_window(2));
    ix=ix(1:par.synth.rf_stride:end); waves.t=trace.t(ix);
    clean=trace.val(ix);
    % Stationary AR(1) Gaussian noise has the configured exponential covariance.
    a=exp(-trace.dt*par.synth.rf_stride/par.synth.rf_corr_seconds);
    noise=randn(numel(ix),1); noise(1)=noise(1)/sqrt(1-a^2);
    noise=sqrt(1-a^2)*filter(1,[1 -a],noise);
    waves.clean_rf(:,j)=clean; %#ok<AGROW>
    waves.rf(:,j)=clean+par.synth.gamma(5)*par.synth.sigma_floor(5)*noise; %#ok<AGROW>
end
waves.rayp_s_km=rayp;
data.rf=struct('x',1,'N',1,'obs',1,'clean',1,'kind','Ps H-kappa stack', ...
    'waves',waves,'dof',par.rf.hk.dof,'nrecords',1,'processing',par.rf.hk);
% obs=1 is the normalized maximum-energy reference, not a waveform sample.
% The true model's ratio need not be exactly one for the H-kappa approximation.
data.rf.true_stack=b3_FORWARD_MODEL_HK(model,data.rf,par);
data.rf.clean=data.rf.true_stack.E_by_Emax;

%% =================== OBSERVED Sp DEPTH STACK AND PARENT ===================
ccp=forward.pred.rf_sp; z=ccp.x; N=numel(z);
% Controlled stack-domain noise; parent is known in this synthetic experiment.
% Real event stacking, migration reference and observed uncertainties belong
% in a separate real-data loader, not hidden inside the candidate likelihood.
a=exp(-median(diff(z))/3); noise=randn(N,1); noise(1)=noise(1)/sqrt(1-a^2);
noise=sqrt(1-a^2)*filter(1,[1 -a],noise);
noise=noise.*ccp_depth_window(z,par.rf.ccp.window_km,par.rf.ccp.taper_km);
data.rf_sp=struct('x',z,'N',N,'obs',ccp.val+par.synth.gamma(6)*par.synth.sigma_floor(6)*noise, ...
    'clean',ccp.val,'parent',ccp.parent,'dof',par.rf.ccp.dof,'nrecords',1, ...
    'kind','Sp CCP depth stack','processing',par.rf.ccp);
% Reference 43-DOF override belongs to the Hopper observed dataset. For this
% synthetic, use dofTHBI's underlying Silver-Chan estimate on its own daughter.
if isnan(data.rf_sp.dof), data.rf_sp.dof=rf_effective_dof(data.rf_sp.obs); end
TRUEmodel.noise_added=any(par.synth.gamma~=0); TRUEmodel.noise_seed=par.synth.seed;
TRUEmodel.noise_type='Four Gaussian curve errors; correlated Ps waveform and Sp depth-stack noise, fixed Sp parent';
TRUEmodel.noise_gamma=par.synth.gamma;
data=a2_PREPARE_DATA(data,par);
copyfile(scratch,fullfile(par.res.dir,'truth_forward'));
forward.output_dir=fullfile(par.res.dir,'truth_forward');
end
