function [data, TRUEmodel, forward] = a2_LOAD_DATA(par)
% Known thin Earth -> five exact predictions -> seeded synthetic noise.
% Truth is saved for evaluation, never an argument of run_one_chain.
%% =================== KNOWN THIN TRUE MODEL ===================
% Prescribed synthetic truth from the existing HBI ECM1 thin test, not a
% measured/recovered real-site model. h in km and Vs in km/s.
h = [0.010; 0.040; 0.080; 0.870; 8; 14; 17; 8; 22];
vs = [0.20; 0.55; 1.50; 2.80; 3.55; 3.75; 3.98; 4.30; 4.45; 4.54];
% Preserve the original HBI Brocher conversion and floor order exactly.
vp = 0.9409 + 2.0947*vs - 0.8206*vs.^2 + 0.2683*vs.^3 - 0.0251*vs.^4;
rho = 1.6612*vp - 0.4721*vp.^2 + 0.0671*vp.^3 - 0.0043*vp.^4 + 0.000106*vp.^5;
vp = max(vp, 1.5*vs);
rho = max(rho, 1.5);
model = [[h; 0], vp, vs, rho];
validate_layer_model(model);
assert(abs(sum(h)-70) < 1e-10 && abs(sum(h(1:7))-40) < 1e-10);
TRUEmodel = struct('case_name', 'ecm1_thin', 'description', ...
    'ECM1 thin synthetic truth: 50 m soft cap, 40 km Moho, 70 km to half-space', ...
    'model_matrix', model, 'theta', [h; vs], 'sed_thickness', .05, ...
    'moho_depth', 40, 'moho_interface', 7, 'total_depth', 70, ...
    'noise_added', false, 'noise_type', 'none', ...
    'source_note', 'Assigned h and Vs from HBI ECM1 thin; Vp/density derived using the same empirical equations.');

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
[~,forward] = b3_FORWARD_MODEL(true_structure,forward_data,par,scratch);
assert(forward.all_forward_succeeded,'Synthetic forward model failed.');
%% =================== OBSERVABLE GRIDS AND BASE ERRORS ===================
coords = {'f','T','T','T'};
for k = 1:5
    d = par.inv.datatypes{k};
    if k<5
        x=par.synth.grids.(d).(coords{k}); clean=forward.pred.(d).val;
    else
        t=forward.pred.rf.t;
        ix=find(t>=par.synth.rf_window(1)-1e-5 & t<=par.synth.rf_window(2)+1e-5);
        ix=ix(1:par.synth.rf_stride:end); x=t(ix); clean=forward.pred.rf.val(ix);
    end
    data.(d).x=x(:); data.(d).clean=clean(:); data.(d).obs=clean(:);
    data.(d).sig=max(par.synth.sigma_floor(k),par.synth.relative_sigma(k)*abs(clean(:)));
    data.(d).R=eye(numel(x));
    if k==5 && par.synth.rf_corr_seconds>0
        data.(d).R=exp(-abs(bsxfun(@minus,x(:),x(:)'))/par.synth.rf_corr_seconds);
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
    if N < 2 || numel(v.x) ~= N || numel(v.sig) ~= N || ...
            any(~isfinite([v.x;v.obs;v.sig])) || any(v.sig <= 0) || ...
            any(diff(v.x) <= 0)
        error('hbi:data','%s requires matching finite arrays, positive sig, increasing x.',d);
    end
    if ~strcmp(d,'rf') && any(v.x <= 0)
        error('hbi:data','Positive frequency/period required.');
    end
    if ~isfield(v,'R'), v.R = eye(N); end
    if ~isequal(size(v.R),[N N]) || any(~isfinite(v.R(:))) || ...
            norm(v.R-v.R','fro') > 1e-10 || any(abs(diag(v.R)-1)>1e-10)
        error('hbi:covariance','%s R must be a finite symmetric correlation matrix.',d);
    end
    % C0 = diag(sig)*R*diag(sig); L0*L0' = C0.
    [L,flag] = chol(v.R,'lower');
    if flag, error('hbi:covariance','%s R must be positive definite.',d); end
    v.N = N;
    v.L0 = bsxfun(@times,v.sig,L); % C0 = diag(sig)*R*diag(sig)
    v.logdetC0 = 2*sum(log(diag(v.L0)));
    data.(d) = v;
end
%% =================== ADD SEEDED SYNTHETIC NOISE ===================
old_rng=rng; restore_rng=onCleanup(@() rng(old_rng)); %#ok<NASGU>
rng(par.synth.seed,'twister');
for k = 1:5
    d=par.inv.datatypes{k};
    data.(d).obs=data.(d).clean+par.synth.gamma(k)*data.(d).L0*randn(data.(d).N,1);
end
%% =================== SAVE TRUE-MODEL SOLVER INPUTS AND OUTPUTS ===================
copyfile(scratch,fullfile(par.res.dir,'truth_forward'));
forward.output_dir=fullfile(par.res.dir,'truth_forward');
end
