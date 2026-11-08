function report = test_thbi_receiver_functions(par, outdir)
% Manual native preflight. CPS only; no MCMC. Preserve every SAC and log.
% Known homogeneous/bilayer Earths test timing, signed amplitudes and stability.
% The production Sp calculation remains visible in b3_FORWARD_MODEL.
if ~exist(outdir,'dir'), mkdir(outdir); end
homogeneous=[35 6.2 3.6 2.7;0 6.2 3.6 2.7];
bilayer=[35 6.2 3.6 2.7;0 8.0 4.5 3.3];
report=struct(); records=struct();

%% =================== Ps: NO FALSE LATE PHASE, CORRECT PRIMARY CONVERSION ===================
cfg=par.forc.rf;
records.ps_homogeneous=forward_receiver_function_cps(homogeneous,cfg,par.forc.cps_bin,fullfile(outdir,'ps_homogeneous'));
records.ps_bilayer=forward_receiver_function_cps(bilayer,cfg,par.forc.cps_bin,fullfile(outdir,'ps_bilayer'));
q=records.ps_homogeneous; late=abs(q.t)>3 & abs(q.t)<90;
report.ps_homogeneous_late_ratio=max(abs(q.val(late)))/max(abs(q.val));
delay=35*(sqrt(1/3.6^2-cfg.rayp_s_km^2)-sqrt(1/6.2^2-cfg.rayp_s_km^2));
q=records.ps_bilayer; ix=find(abs(q.t-delay)<.9); [amp,j]=max(abs(q.val(ix)));
report.ps_expected_delay=delay; report.ps_picked_delay=q.t(ix(j));
report.ps_conversion_amplitude=amp;

%% =================== Sp: COMPONENTS, HOMOGENEOUS CONTROL AND NEGATIVE LAG ===================
cfg=par.forc.rf_sp;
records.sp_homogeneous=components(homogeneous,cfg,par.forc.cps_bin,fullfile(outdir,'sp_homogeneous'));
records.sp_bilayer=components(bilayer,cfg,par.forc.cps_bin,fullfile(outdir,'sp_bilayer'));
% The direct CPS -S ratio is used ONLY as a stable homogeneous control here.
% It is never used as the production layered Sp calculator.
control=cfg; control.component_gaussian_alpha=cfg.gaussian_alpha;
records.sp_direct_control=raw_component(homogeneous,control,par.forc.cps_bin, ...
    fullfile(outdir,'sp_direct_homogeneous'),'');
q=records.sp_homogeneous; direct=records.sp_direct_control;
report.sp_homogeneous_direct_error=norm(q.val-direct.val)/norm(direct.val);
late=abs(q.t)>4 & abs(q.t)<90;
report.sp_homogeneous_late_ratio=max(abs(q.val(late)))/max(abs(q.val));
delay=-35*(sqrt(1/3.6^2-cfg.rayp_s_km^2)-sqrt(1/6.2^2-cfg.rayp_s_km^2));
q=records.sp_bilayer; ix=find(abs(q.t-delay)<.9); [amp,j]=max(abs(q.val(ix)));
report.sp_expected_delay=delay; report.sp_picked_delay=q.t(ix(j));
report.sp_conversion_amplitude=amp;

%% =================== LONG-RECORD AND WATER-LEVEL SENSITIVITY ===================
% Double dt at fixed 4096 samples to double duration; the Sp high-frequency
% corner remains resolved. Compare the configured full passband, through
% 60 or 100 s. This also tests sensitivity to the coarser component sampling.
long=cfg; long.dt=2*cfg.dt;
records.sp_long=components(bilayer,long,par.forc.cps_bin,fullfile(outdir,'sp_long_record'));
records.sp_bilayer=receiver_function_bands(records.sp_bilayer,cfg);
records.sp_long=receiver_function_bands(records.sp_long,long);
q=records.sp_bilayer; selected=q.t>=-60 & q.t<=60;
reference=q.val(selected);
long_values=interp1(records.sp_long.t,records.sp_long.val,q.t(selected),'linear');
report.sp_long_record_relative_error=norm(long_values-reference)/norm(reference);
R=read_cps_sac(fullfile(outdir,'sp_bilayer','r','hrftn96.sac'));
Z=read_cps_sac(fullfile(outdir,'sp_bilayer','z','hrftn96.sac'));
low=cfg; low.water_level=cfg.water_level/2;
records.sp_lower_waterlevel=receiver_function_bands(waterlevel_receiver_function(R,Z,low),low);
report.sp_waterlevel_relative_error=norm(records.sp_lower_waterlevel.val(selected)-reference)/norm(reference);

%% =================== DEEP INTERFACE CONTROL FOR THE 250-KM PROFILE ===================
% An independent two-material control, not a modification of the ECM1 truth.
% Check that native calculations retain a conversion from the new model base.
if par.mod.maxz>=250
    deep=[par.mod.maxz 8.0 4.5 3.3;0 8.2 4.65 3.4];
    for name={'rf','rf_sp'}
        d=name{1}; cfg=par.forc.(d);
        if strcmp(d,'rf')
            q=forward_receiver_function_cps(deep,cfg,par.forc.cps_bin,fullfile(outdir,'deep_ps'));
            sign_delay=1;
        else
            q=components(deep,cfg,par.forc.cps_bin,fullfile(outdir,'deep_sp'));
            sign_delay=-1;
        end
        delay=sign_delay*par.mod.maxz*(sqrt(1/4.5^2-cfg.rayp_s_km^2)-sqrt(1/8^2-cfg.rayp_s_km^2));
        ix=find(abs(q.t-delay)<1); [amp,j]=max(abs(q.val(ix)));
        report.(d).deep_expected_delay=delay;
        report.(d).deep_picked_delay=q.t(ix(j));
        report.(d).deep_conversion_amplitude=amp;
        records.(['deep_' d])=receiver_function_bands(q,cfg);
    end
end
report.note='Numerical controls, not proof of geological recovery or uncertainty calibration.';
save(fullfile(outdir,'receiver_function_checks.mat'),'report','records','homogeneous','bilayer','-v7');

%% =================== GATES BEFORE ANY INVERSION ===================
assert(report.ps_homogeneous_late_ratio<.02 && report.sp_homogeneous_late_ratio<.02, ...
    'thbi:rf_test','Homogeneous Earth produced late conversion energy; inspect saved SAC files.');
assert(abs(report.ps_picked_delay-report.ps_expected_delay)<.7 && report.ps_conversion_amplitude>1e-4, ...
    'thbi:rf_test','Bilayer Ps timing/amplitude failed.');
assert(abs(report.sp_picked_delay-report.sp_expected_delay)<.7 && report.sp_conversion_amplitude>1e-4, ...
    'thbi:rf_test','Bilayer Sp negative-lag conversion failed.');
assert(report.sp_homogeneous_direct_error<.05, ...
    'thbi:rf_test','Sp component deconvolution disagrees with homogeneous CPS Z/R.');
assert(report.sp_long_record_relative_error<.05 && report.sp_waterlevel_relative_error<.10, ...
    'thbi:rf_test','Sp changes materially with record length/water level; inspect preflight before sampling.');
if par.mod.maxz>=250
    for name={'rf','rf_sp'}
        q=report.(name{1});
        assert(abs(q.deep_picked_delay-q.deep_expected_delay)<.7 && q.deep_conversion_amplitude>1e-4, ...
            'thbi:rf_test','Deep Ps/Sp conversion control failed; inspect the saved native records.');
    end
end
fprintf('PASS: native Ps/Sp homogeneous, bilayer timing, amplitude, duration and water-level controls.\n');
end

function rf=components(model,cfg,cps_bin,folder)
R=raw_component(model,cfg,cps_bin,fullfile(folder,'r'),'r');
Z=raw_component(model,cfg,cps_bin,fullfile(folder,'z'),'z');
rf=waterlevel_receiver_function(R,Z,cfg);
end

function trace=raw_component(model,cfg,cps_bin,folder,flag)
if ~exist(folder,'dir'), mkdir(folder); end
old=pwd; restore=onCleanup(@() cd(old)); %#ok<NASGU>
cd(folder); write_elastic_cps_model('rf_model.mod',model);
option=''; if ~isempty(flag), option=['-' flag]; end
cmd=sprintf('%s -S %s -2 -M rf_model.mod -RAYP %.10g -ALP %.10g -DT %.10g -NSAMP %d -D %.10g', ...
    quote_shell_path(fullfile(cps_bin,'hrftn96')),option,cfg.rayp_s_km, ...
    cfg.component_gaussian_alpha,cfg.dt,cfg.npts,cfg.delay);
if exist('hrftn96.sac','file'), delete('hrftn96.sac'); end
[status,message]=system([cmd ' > hrftn96.log 2>&1']);
assert(status==0 && exist('hrftn96.sac','file')==2,'thbi:rf_test', ...
    'CPS component failed: %s. See %s.',message,folder);
trace=read_cps_sac('hrftn96.sac');
assert(numel(trace.val)==cfg.npts && abs(trace.dt-cfg.dt)<1e-6 && ...
    abs(trace.t(1)+cfg.delay)<1e-5 && abs(trace.rayp_s_km-cfg.rayp_s_km)<1e-6 && ...
    abs(trace.gaussian_alpha-cfg.component_gaussian_alpha)<1e-5, ...
    'thbi:rf_test','CPS returned unexpected component headers.');
end
