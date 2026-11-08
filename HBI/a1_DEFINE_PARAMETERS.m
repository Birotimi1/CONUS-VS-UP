function par = a1_DEFINE_PARAMETERS(profile)
% Fixed-layer HBI parameters. Units: km, km/s, g/cm^3, seconds, Hz.
% Organization follows THBI_MCMC; layer/noise priors follow our HBI.
% See docs/HBI_DESIGN.md for equations, attribution and deliberate changes.
if nargin == 0, profile = 'pilot'; end
[repo,paths] = a0_STARTUP_BAYES();
%% Chain settings
par.inv.datatypes = {'hvsr','ellip','cph','ugr','rf'};
par.inv.profile = profile;
par.inv.nchains = 2;
par.inv.niter = 240;
par.inv.burnin = 120;
par.inv.seed = 31000;
par.inv.parallel = false; % Set true for MATLAB process workers (PCT).
par.inv.checkpoint_every = 40;
par.inv.print_every = 40;
par.inv.noise_move_probability = .25;
par.inv.adapt = true; % Robbins-Monro proposal tuning, only during burn-in.
par.inv.max_consecutive_failures = 10;
if strcmp(profile,'production')
    par.inv.nchains = 4;
    par.inv.niter = 30000;
    par.inv.burnin = 10000;
    par.inv.checkpoint_every = 250;
    par.inv.print_every = 250;
elseif ~strcmp(profile,'pilot')
    error('hbi:profile','Choose pilot or production.');
end

%% Layered Earth: nine finite layers over one half-space
par.mod.nlayer = 10;
par.mod.idx.h = 1:9; par.mod.idx.vs = 10:19;
par.mod.names = [arrayfun(@(k)sprintf('h%d',k),1:9,'UniformOutput',false), ...
    arrayfun(@(k)sprintf('Vs%d',k),1:10,'UniformOutput',false)];
par.mod.lower = [.005;.02;.08;.30;1.5;2;4;6;8; ...
    .12;.25;.5;.8;2;2.6;3;3.4;3.8;4.1];
par.mod.upper = [.08;.30;1.5;4;10;15;18;22;28; ...
    .9;1.6;2.2;3.2;3.8;4;4.3;4.6;4.9;5];
par.mod.step = [.003;.010;.030;.15;.40;.50;.60;.80;.80; ...
    .02;.03;.05;.06;.06;.05;.05;.05;.05;.05];
par.mod.monotonic = true;
% Explicit reference start, different from the prescribed truth. Every
% chain is independently jittered. No truth is passed to the sampler.
par.mod.start = [.018;.070;.150;.80;7;12;16;9;20; ...
    .28;.70;1.70;2.65;3.40;3.65;3.90;4.22;4.40;4.65];
par.mod.start_jitter = 1;

%% Hierarchical noise: eta = ln(gamma), C_d = gamma_d^2 C0_d
par.mod.noise.log_center = zeros(5,1);
par.mod.noise.log_width = 2*ones(5,1); % natural-log units, not decades
par.mod.noise.lower = log(.1)*ones(5,1);
par.mod.noise.upper = log(50)*ones(5,1);
par.mod.noise.step = .12*ones(5,1);

%% Exact forward solvers shared by synthetic generation and inversion
par.forc.project_root = repo;
par.forc.shared_project = paths.shared_project;
par.forc.cps_bin = paths.cps_bin;
par.forc.hvf_bin = paths.hvf_bin;
par.forc.max_layer_thickness_km = 2;
par.forc.verbose = false;
par.forc.rf = struct('rayp_s_km',.06,'gaussian_alpha',2.5, ...
    'dt',.025,'npts',4096,'delay',10,'plot_window',[-2 25]);

%% Synthetic data and known base uncertainties
par.synth.seed = 9041;
par.synth.gamma = ones(5,1);
par.synth.grids.hvsr.f = logspace(log10(.2),log10(25),80)';
par.synth.grids.ellip.T = linspace(3,25,50)';
par.synth.grids.cph.T = linspace(9,35,35)';
par.synth.grids.ugr.T = linspace(3,25,40)';
par.synth.rf_window = [-2 25];
par.synth.rf_stride = 8; % select every eighth native point: 0.2 s spacing
par.synth.rf_corr_seconds = .4; % exponential covariance of ADDED RF noise
par.synth.relative_sigma = [.05;.03;0;0;0];
par.synth.sigma_floor = [.05;.015;.03;.04;.015];

%% Persistent output on filer, disposable solver files on local disk
par.res.root = paths.inversions;
par.res.scratch_root = paths.scratch;
par.res.show_figures = true;
par.res.depth_km = unique([linspace(0,.2,101),linspace(.2,2,91),linspace(2,100,197)])';
end
