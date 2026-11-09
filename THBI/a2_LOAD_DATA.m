function [data, TRUEmodel, forward] = a2_LOAD_DATA(par,truth_only)
% Known thin Earth -> all configured exact predictions -> seeded synthetic noise.
% Truth is saved for evaluation, never an argument of run_one_chain.
% Optional truth_only=true exposes the complete known structure for solver-free
% checks; no data are generated and forward.status is explicitly not_computed.
if nargin<2, truth_only=false; end
assert(islogical(truth_only) && isscalar(truth_only),'thbi:truth','truth_only must be a logical scalar.');
%% =================== KNOWN THIN TRUE MODEL ===================
% Prescribed synthetic truth from the existing HBI ECM1 thin test, not a
% measured/recovered real-site model. h in km and Vs in km/s.
h = [0.010; 0.040; 0.080; 0.870; 8; 14; 17; 8; 22];
vs = [0.20; 0.55; 1.50; 2.80; 3.55; 3.75; 3.98; 4.30; 4.45; 4.54];
assert(par.mod.maxz>=70,'thbi:truth','The thin truth requires a model base at or below 70 km.');
zone_counts=[3 4 2];
if par.mod.maxz>70
    % Continue the old half-space unchanged. Its new finite segment adds no
    % impedance jump at the deeper model base and preserves the thin cap.
    h=[h;par.mod.maxz-70]; vs=[vs;vs(end)]; zone_counts(3)=3;
end
% Preserve the original HBI Brocher conversion and floor order exactly.
vp = 0.9409 + 2.0947*vs - 0.8206*vs.^2 + 0.2683*vs.^3 - 0.0251*vs.^4;
rho = 1.6612*vp - 0.4721*vp.^2 + 0.0671*vp.^3 - 0.0043*vp.^4 + 0.000106*vp.^5;
vp = max(vp, 1.5*vs);
rho = max(rho, 1.5);
model = [[h; 0], vp, vs, rho];
validate_layer_model(model);
assert(abs(sum(h)-par.mod.maxz) < 1e-10 && abs(sum(h(1:7))-40) < 1e-10);
TRUEmodel = struct('case_name', 'ecm1_thin', 'description', ...
    sprintf('ECM1 thin truth: 50 m soft cap, 40 km Moho, model base %g km',par.mod.maxz), ...
    'model_matrix', model, 'theta', [h; vs], 'sed_thickness', .05, ...
    'moho_depth', 40, 'moho_interface', 7, 'total_depth', par.mod.maxz, ...
    'noise_added', false, 'noise_type', 'none', ...
    'source_note', 'Assigned h and Vs from HBI ECM1 thin; Vp/density derived using the same empirical equations.');
TRUEmodel.sediment_base_depth = .13;
TRUEmodel.zone_counts = zone_counts; % finite layers; half-space is separate

%% =================== CHECK THAT THE TRUTH IS IN THE STRUCTURAL PRIOR ===================
assert(abs(TRUEmodel.total_depth-par.mod.maxz)<1e-10, ...
    'thbi:truth','The thin truth and the fixed half-space depth must agree.');
true_parameters.zsed=TRUEmodel.sediment_base_depth;
true_parameters.zmoh=TRUEmodel.moho_depth;
true_parameters.halfspaceVS=vs(end);
true_parameters.datahparm.eta=par.mod.noise.log_center;
zone_edges=[0,true_parameters.zsed,true_parameters.zmoh,par.mod.maxz];
first=1;
for z=1:3
    rows=first:first+TRUEmodel.zone_counts(z)-1;
    assert(abs(sum(h(rows))-diff(zone_edges(z:z+1)))<1e-10, ...
        'thbi:truth','True layer counts do not match the zone boundary depths.');
    fractions=cumsum(h(rows(1:end-1)))/diff(zone_edges(z:z+1));
    true_parameters.(par.mod.zone_fields{z})=struct('f',fractions(:),'VS',vs(rows));
    first=rows(end)+1;
end
assert(first==numel(h)+1,'thbi:truth','Zone counts do not cover all finite layers.');
[supported,why]=a1_TEST_CONDITIONS(true_parameters,par);
assert(supported,'thbi:truth','True structure is outside the configured prior: %s',why);
true_parameters=make_mod_from_parms(true_parameters,par);
assert(max(abs(true_parameters.matrix(:)-model(:)))<1e-10, ...
    'thbi:truth','Zoned representation changed the true Earth.');
TRUEmodel.zoned_parameters=true_parameters;
if truth_only
    data=struct(); forward=struct('status','not_computed'); return;
end

%% =================== NATIVE FORWARD CALCULATION ===================
scratch = tempname(par.res.scratch_root); mkdir(scratch);
cleanup = onCleanup(@() rmdir(scratch,'s')); %#ok<NASGU>
% Use b3 for both generation and inversion: edits to its solver blocks
% therefore apply to both sides of the synthetic test.
forward_data.hvsr.x = par.synth.grids.hvsr.f;
forward_data.ellip.x = par.synth.grids.ellip.T;
forward_data.cph.x = par.synth.grids.cph.T;
forward_data.ugr.x = par.synth.grids.ugr.T;
true_structure.matrix = TRUEmodel.model_matrix;
try
    [~,forward] = b3_FORWARD_MODEL(true_structure,forward_data,par,scratch);
catch problem
    context=struct('stage','synthetic truth','iteration',0);
    preserve_forward_failure(fullfile(par.res.dir,'truth_forward_failure'), ...
        true_structure,forward_data,par,problem,scratch,context);
    rethrow(problem);
end
assert(forward.all_forward_succeeded,'Synthetic forward model failed.');
%% =================== OBSERVABLE GRIDS AND BASE ERRORS ===================
coords = {'f','T','T','T'};
for k = 1:numel(par.inv.datatypes)
    d = par.inv.datatypes{k};
    if k<5
        x=par.synth.grids.(d).(coords{k}); clean=forward.pred.(d).val;
    else
        t=forward.pred.(d).t; window=par.synth.([d '_window']);
        stride=par.synth.([d '_stride']);
        assert(window(1)>=t(1)-2e-5 && window(2)<=t(end)+2e-5, ...
            'thbi:truth','%s window exceeds the native record.',d);
        ix=find(t>=window(1)-2e-5 & t<=window(2)+2e-5);
        ix=ix(1:stride:end); x=t(ix); clean=forward.pred.(d).val(ix);
        data.(d).processing=par.forc.(d);
        data.(d).kind=forward.pred.(d).kind;
    end
    data.(d).x=x(:); data.(d).clean=clean(:); data.(d).obs=clean(:);
    data.(d).sig=max(par.synth.sigma_floor(k),par.synth.relative_sigma(k)*abs(clean(:)));
    data.(d).R=eye(numel(x));
    if k>=5 && par.synth.([d '_corr_seconds'])>0
        data.(d).R=exp(-abs(bsxfun(@minus,x(:),x(:)'))/par.synth.([d '_corr_seconds']));
    end
end
%% =================== FIXED COVARIANCE AND WHITENING ===================
% Fix coordinates, uncertainties, and covariance BEFORE starting any chain.
% No proposal-dependent masks, truncation, normalization, or dropped data.

%% Validate each dataset and build its fixed whitening matrix
for k = 1:numel(par.inv.datatypes)
    d = par.inv.datatypes{k};
    if ~isfield(data,d), error('hbi:data','Missing %s data.',d); end
    v = data.(d); v.x = v.x(:); v.obs = v.obs(:); v.sig = v.sig(:);
    N = numel(v.obs);
    if ~isreal(v.x) || ~isreal(v.obs) || ~isreal(v.sig) || ...
            N < 2 || numel(v.x) ~= N || numel(v.sig) ~= N || ...
            any(~isfinite([v.x;v.obs;v.sig])) || any(v.sig <= 0) || ...
            any(diff(v.x) <= 0)
        error('hbi:data','%s requires matching finite arrays, positive sig, increasing x.',d);
    end
    if ~any(strcmp(d,{'rf','rf_sp'})) && any(v.x <= 0)
        error('hbi:data','Positive frequency/period required.');
    end
    if ~isfield(v,'R'), v.R = eye(N); end
    if ~isreal(v.R) || ~isequal(size(v.R),[N N]) || any(~isfinite(v.R(:))) || ...
            norm(v.R-v.R','fro') > 1e-10 || any(abs(diag(v.R)-1)>1e-10)
        error('hbi:covariance','%s R must be a finite symmetric correlation matrix.',d);
    end
    % C0 = diag(sig)*R*diag(sig); L0*L0' = C0.
    [L,flag] = chol(v.R,'lower');
    if flag, error('hbi:covariance','%s R must be positive definite.',d); end
    v.N = N;
    v.L0 = bsxfun(@times,v.sig,L); % C0 = diag(sig)*R*diag(sig)
    v.logdetC0 = 2*sum(log(diag(v.L0)));
    assert(all(isfinite(v.L0(:))) && all(diag(v.L0)>0) && isfinite(v.logdetC0), ...
        'hbi:covariance','%s covariance scaling is not numerically representable.',d);
    data.(d) = v;
end
%% =================== ADD SEEDED SYNTHETIC NOISE ===================
old_rng=rng; restore_rng=onCleanup(@() rng(old_rng)); %#ok<NASGU>
rng(par.synth.seed,'twister');
for k = 1:numel(par.inv.datatypes)
    d=par.inv.datatypes{k};
    data.(d).obs=data.(d).clean+par.synth.gamma(k)*data.(d).L0*randn(data.(d).N,1);
end
TRUEmodel.noise_added=any(par.synth.gamma~=0);
TRUEmodel.noise_type='Gaussian with fixed within-dataset correlation; independent across datasets';
TRUEmodel.noise_gamma=par.synth.gamma;
TRUEmodel.noise_seed=par.synth.seed;
%% =================== SAVE TRUE-MODEL SOLVER INPUTS AND OUTPUTS ===================
copyfile(scratch,fullfile(par.res.dir,'truth_forward'));
forward.output_dir=fullfile(par.res.dir,'truth_forward');
end
