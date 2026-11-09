function report=test_thbi_extended100()
% Solver-free checks of the ACTUAL 250-km truth, grids and 100-s RF filters.
% Native CPS/HVf controls are in DATA/check_thbi_extended100.m.
par=bayes_inv_parms('extended100'); validate_configuration(par);
old=bayes_inv_parms('extended60');
assert(par.mod.maxz==250 && par.inv.nchains==2 && ...
    par.inv.niter==10000 && par.inv.burnin==4000);
assert(strcmp(par.proposal.dr.mode,'layered_path') && par.proposal.dr.same_zone_weight==0);
assert(strcmp(old.proposal.dr.mode,'layered_path') && old.proposal.dr.same_zone_weight==0);

%% =================== ACTUAL KNOWN EARTH AND PRIOR SUPPORT ===================
[data,truth,forward]=a2_LOAD_DATA(par,true);
assert(isempty(fieldnames(data)) && strcmp(forward.status,'not_computed'));
assert(truth.total_depth==250 && isequal(truth.zone_counts,[3 4 3]));
assert(abs(sum(truth.model_matrix(:,1))-250)<1e-10);
assert(max(abs(cumsum(truth.model_matrix(1:4,1))-[.01;.05;.13;1]))<1e-12);
assert(truth.model_matrix(end-1,1)==180 && ...
    isequal(truth.model_matrix(end-1,2:4),truth.model_matrix(end,2:4)));
assert(a1_TEST_CONDITIONS(truth.zoned_parameters,par));
assert(isfinite(a3_CALC_LOG_PRIOR(truth.zoned_parameters,par)));
refined=refine_layer_model(truth.model_matrix,par.forc.max_layer_thickness_km);
assert(all(refined(:,1)<=2) && size(refined,1)<200);
% Identical finite mantle/half-space are merged by numerical refinement;
% the physical inversion domain still reaches 250 km.
start_par=par; start_par.mod.start.jitter=false;
start=b1_INITIATE_MODEL(start_par);
assert(abs(sum(start.matrix(:,1))-250)<1e-10 && ~isequal(start.K,truth.zone_counts));
deep_refined=refine_layer_model(start.matrix,par.forc.max_layer_thickness_km);
assert(size(deep_refined,1)<200 && abs(sum(deep_refined(:,1))-250)<1e-10);

%% =================== LONG PERIODS WITHOUT LOSING SHALLOW SAMPLES ===================
f=par.synth.grids.hvsr.f;
assert(abs(min(f)-.01)<1e-14 && abs(max(f)-25)<1e-12 && all(diff(f)>0));
assert(all(ismember(old.synth.grids.hvsr.f,f)));
for name={'ellip','cph','ugr'}
    d=name{1}; T=par.synth.grids.(d).T;
    assert(max(T)==100 && all(diff(T)>0) && all(ismember(old.synth.grids.(d).T,T)));
end

%% =================== TRUE 250-KM CONVERSION LAGS INSIDE DATA WINDOWS ===================
report=struct('model_base_km',250,'longest_period_s',100);
for name={'rf','rf_sp'}
    d=name{1}; cfg=par.forc.(d); stride=par.synth.([d '_stride']);
    window=par.synth.([d '_window']);
    assert(cfg.npts==4096 && cfg.npts*cfg.dt>=800 && ...
        cfg.dt*stride==.2 && cfg.period_band_s(2)==100);
    assert(isequal(cfg.diagnostic_period_bands_s(end,:),[60 100]));
    m=truth.model_matrix; p=cfg.rayp_s_km;
    delay=sum(m(:,1).*(sqrt(1./m(:,3).^2-p^2)-sqrt(1./m(:,2).^2-p^2)));
    if strcmp(d,'rf_sp'), delay=-delay; end
    assert(isreal(delay) && delay>window(1) && delay<window(2));
    report.([d '_base_conversion_lag_s'])=delay;
    % The base is not a contrast in this truth, so this is a travel-time
    % coverage check, not a claimed detectable converted arrival.

    %% Isolate an 80-s sinusoid and reject 10-s energy in the new 60-100 band
    t=((0:cfg.npts-1)'-cfg.npts/2)*cfg.dt;
    trace=struct('val',sin(2*pi*t/80)+sin(2*pi*t/10),'dt',cfg.dt,'t',t);
    filtered=receiver_function_bands(trace,cfg);
    center=abs(t)<150;
    basis=[sin(2*pi*t/80),cos(2*pi*t/80),sin(2*pi*t/10),cos(2*pi*t/10)];
    fit=basis(center,:)\filtered.band_val(center,end);
    expected=1/sqrt(1+(80/100)^8)/sqrt(1+(60/80)^8);
    assert(abs(hypot(fit(1),fit(2))-expected)<.04 && hypot(fit(3),fit(4))<.01);
    report.([d '_80s_band_gain'])=hypot(fit(1),fit(2));
end

%% =================== SIGNED Sp AT THE NEW NATIVE SAMPLING ===================
cfg=par.forc.rf_sp; N=cfg.npts;
parent=struct('val',zeros(N,1),'dt',cfg.dt,'t',(0:N-1)'*cfg.dt-cfg.delay);
daughter=parent; index=round(cfg.delay/cfg.dt)+1;
parent.val(index)=1; daughter.val(index-round(12/cfg.dt))=-.4;
rf=waterlevel_receiver_function(parent,daughter,cfg);
[amp,j]=min(rf.val);
assert(abs(rf.t(j)+12)<cfg.dt && abs(amp+.4*cfg.gaussian_alpha/sqrt(pi))<1e-5);
fprintf('PASS: 250-km true/prior Earth, all 100-s grids, RF duration, 80-s gain and signed Sp.\n');
fprintf('Base travel-time coverage: Ps %.3f s, Sp %.3f s (no true contrast at 250 km).\n', ...
    report.rf_base_conversion_lag_s,report.rf_sp_base_conversion_lag_s);
end
