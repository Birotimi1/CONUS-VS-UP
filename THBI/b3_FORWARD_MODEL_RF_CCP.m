function ccp = b3_FORWARD_MODEL_RF_CCP(matrix, components, par)
% CPS R/Z -> P/SV -> simple source parent -> bandpass -> Sp depth stack.
% Layered isotropic adaptation of Zach/Brennan b3_FORWARD_MODEL_RF_ccp.

%% =================== PHYSICAL COMPONENT ROTATION AND SIGNED PARENT ===================
cfg=par.rf.ccp; native=par.forc.rf_sp;
r=components{1}; z=components{2}; t=r.t(:);
assert(numel(t)==numel(z.t) && max(abs(t-z.t))<1e-6, ...
    'thbi:ccp','CPS component coordinates disagree.');
% CPS Uz is positive UP. Eilon rotation requires Z positive DOWN.
[P,SV]=Rotate_XZ_to_PSV(r.val(:),-z.val(:),matrix(1,2),matrix(1,3),native.rayp_s_km);
[amplitude,peak]=max(abs(SV));
assert(isfinite(amplitude) && amplitude>0,'thbi:ccp','Missing Sp parent arrival.');
P=P/SV(peak); SV=SV/SV(peak); t=t-t(peak);
sigma=cfg.source_period_s/6;
parent=exp(-.5*(t/sigma).^2);
% CPS Gaussian exp(-omega^2/(4*alpha^2)) must match that source pulse.
assert(abs(native.component_gaussian_alpha-1/(sqrt(2)*sigma))<1e-10, ...
    'thbi:ccp','CPS component source and simple parent pulse widths disagree.');

%% =================== COMMON PERIOD FILTER, BEFORE MIGRATION ===================
wave=struct('t',t,'dt',r.dt,'val',P);
wave=receiver_function_bands(wave,native);
parent_wave=receiver_function_bands(struct('t',t,'dt',r.dt,'val',parent),native);
% Keep the padded record for the 100-s band; do not crop it to a 50-s trace.
depth=cfg.depth_km(:);
lag=layered_conversion_time(matrix,depth,native.rayp_s_km,'Sp');
assert(all(lag>=t(1) & lag<=t(end)),'thbi:ccp','Native record does not cover CCP depths.');
daughter=interp1(t,wave.val,lag,'linear');
parent_depth=interp1(t,parent_wave.val,lag,'linear');

%% =================== DAUGHTER WINDOW AND FIXED DEPTH COORDINATES ===================
daughter=daughter.*ccp_depth_window(depth,cfg.window_km,cfg.taper_km);
ccp=struct('x',depth,'val',daughter,'parent',parent_depth, ...
    'PSV',[daughter parent_depth],'kind','Sp CCP depth stack: P daughter / SV parent', ...
    'status','ok','model_input',matrix,'observed',false, ...
    'waveform',wave,'unfiltered_parent',SV,'native_time',t, ...
    'rayp_s_km',native.rayp_s_km);
end
