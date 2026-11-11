function par = bayes_inv_parms(profile)
% Zoned, layered transdimensional HBI. Units: km, km/s, g/cm^3, s, Hz.
% Edit this file before running DATA/run_ecm1_thin_thbi.m in MATLAB.
if nargin == 0, profile = 'pilot'; end
%% =================== THBI LINKED RADIAL ANISOTROPY ===================
if strcmp(profile,'reference100_vti')
    par=bayes_inv_parms('reference100'); par.inv.profile=profile;
    par.mod.anisotropy='vti';
    % Mirror Zach/Brennan active example: isotropic sediments/mantle and one
    % free crustal xi. Solvers also support mantle xi if these bounds change.
    par.mod.xi_bounds=[1 1;.8 1.2;1 1]; par.mod.start.xi=[1 1 1];
    par.mod.psi_to_xi=-1; par.mod.elastic_eta=1;
    par.synth.xi=[1 1.06 1]; % assigned synthetic perturbation, NOT ECM1/AK135 data
    par.proposal.names{end+1}='anisotropy';
    par.proposal.probability=[.94*par.proposal.probability .06];
    par.proposal.xi_step=.01;
    par.forc.hvf_vti_bin=fullfile(par.forc.project_root,'native','HVf_VTI','bin');
    % Complete modal sum: 4096 is a checked capacity, not a truncation choice.
    % Native search stops at the physical cutoff and fails if this cap is hit.
    par.forc.hvsr=struct('nmr',4096,'nml',4096,'nks',4000,'prec',1e-8, ...
        'apsv',1e-3,'ash',1e-3,'dphase',pi/512,'require_complete',true);
    validate_configuration(par); return;
end
%% =================== LAYERED THBI RF PRODUCTS AND PUBLISHED TRUTH ===================
if strcmp(profile,'reference100')
    par=bayes_inv_parms('extended100'); par.inv.profile=profile;
    par.mod.material='thbi'; par.mod.vpvs_bounds=[1.6 2.1];
    par.mod.anisotropy='isotropic'; % current CPS/HVf implementation
    par.mod.start.vpvs=1.76; par.mod.start.method='prior';
    par.mod.start.max_attempts=10000;
    par.mod.boundary_lower=[.001 15]; par.mod.boundary_upper=[5 70];
    par.mod.vsmin=[.5 2.5 3.7]; par.mod.vsmax=[3.3 4.5 5.1];
    par.conditions.moho_gradient_window_km=4;
    par.conditions.moho_gradient_bounds=[-.008 .1];
    par.proposal.names=[par.proposal.names,{'vpvs','boundary_mean','boundary_jump','moho_hk'}];
    par.proposal.probability=[.23 .13 .08 .13 .13 .14 .06 .03 .03 .04];
    par.proposal.vpvs_step=.015;
    par.proposal.coupled_vs_step=.04; par.proposal.moho_hk_step=[.3 .008];
    par.rf.method='hk_ccp';
    par.rf.hk=struct('H_km',linspace(15,70,200),'kappa',linspace(1.6,2.1,201), ...
        'phase_weights',[.5 .3 .2],'dof',15,'error_scale',1,'min_error',.002, ...
        'rayp_s_km',[.045 .06 .075]);
    par.rf.ccp=struct('depth_km',(-15:.5:300)','dof',NaN, ...
        'dof_method','silver_chan_observed_daughter', ...
        'window_km',[20 250],'taper_km',10,'source_period_s',2.5, ...
        'weight_depth_val',[-15 .3;30 .3;70 1;6371 1],'weight_smooth_km',80);
    par.forc.rf_sp.component_gaussian_alpha=sqrt(18)/par.rf.ccp.source_period_s;
    par.forc.rf.delay=250; par.synth.rf_window=[-100 500];
    par.forc.rf.plot_window=[-5 35];
    par.synth.truth='ecm1_19178_ak135';
    % First four errors remain HBI relative multipliers; RFs use absolute sigma
    % and reference effective degrees of freedom, with proper log-uniform priors.
    par.mod.noise.log_uniform=[false(4,1);true(2,1)];
    par.mod.noise.log_center(5:6)=log([.3;.2]);
    par.mod.noise.lower(5:6)=log([.01;.01]);
    par.mod.noise.upper(5:6)=log([2;2]);
    validate_configuration(par);
    return;
end
[repo,paths] = a0_STARTUP_BAYES();

%% =================== CHAINS AND FIXED MOVE PROBABILITIES ===================
par.inv.datatypes = {'hvsr','ellip','cph','ugr','rf'};
par.inv.profile = profile;
par.inv.nchains = 2;
par.inv.niter = 1000;
par.inv.burnin = 500;
par.inv.seed = 42000;
par.inv.parallel = false;
par.inv.checkpoint_every = 50;
par.inv.print_every = 50;
par.inv.prior_only = false; % true omits likelihood/forward calls (for checks)
par.inv.max_consecutive_failures = 10;
if strcmp(profile,'production')
    par.inv.nchains = 4;
    par.inv.niter = 30000;
    par.inv.burnin = 10000;
    par.inv.checkpoint_every = 250;
    par.inv.print_every = 250;
elseif ~any(strcmp(profile,{'pilot','extended60','extended100'}))
    error('thbi:profile','Choose pilot, production, extended60 or extended100.');
end
% Never renormalize these probabilities at a count bound: invalid moves stay
% at the current state. This keeps forward/reverse selection probabilities explicit.
par.proposal.names = {'velocity','interface','boundary','birth','death','noise'};
par.proposal.probability = [.30 .15 .10 .15 .15 .15];
par.proposal.zone_probability = [1 1 1]/3;
par.proposal.halfspace_probability = .15; % within the velocity-move family
par.proposal.vs_step = [.04 .06 .06];
par.proposal.interface_step = [.04 .04 .04]; % zone-fraction coordinates
par.proposal.boundary_step = [.03 1.0];      % basement and Moho, km
par.proposal.birth_vs_sigma = [.08 .08 .08]; % untruncated Gaussian contrast
par.proposal.halfspace_step = .05;
par.proposal.noise_step = .12;
par.proposal.birth_broad_probability = 0; % legacy pilot: one Gaussian
par.proposal.birth_vs_sigma_broad = [.4 .3 .2];
par.proposal.dr = struct('enabled',false,'mode','fixed_scalar','scale',.25,'same_zone_weight',.5);
% Proposal scales are fixed for this first RJ implementation, including burn-in.

%% =================== ZONES AND VARIABLE PHYSICAL LAYER COUNTS ===================
par.mod.zone_names = {'sediment','crust','mantle'};
par.mod.zone_fields = {'sedmparm','crustmparm','mantmparm'};
par.mod.kmin = [1 1 1];
par.mod.kmax = [8 10 8];
par.mod.count_prior = 'jeffreys'; % base mass proportional to 1/K; or 'uniform'
par.mod.min_fraction = [.002 .001 .001];
% Fractions are relative to each zone's thickness, not absolute depths.
% With K layers, every fractional thickness is >= min_fraction.
par.mod.vsmin = [.12 2.4 3.6];
par.mod.vsmax = [3.2 4.3 5.0];
par.mod.boundary_lower = [.05 20];   % sediment base and Moho
par.mod.boundary_upper = [2.0 60];
par.mod.maxz = 70;                   % fixed half-space top in this thin test
par.mod.halfspace_bounds = [4.1 5.1];

%% =================== THBI-STYLE PHYSICAL CONDITIONS ===================
par.conditions.monotonic = [true false false]; % monotone sediment; crust/mantle LVZs allowed
par.conditions.pos_sed2basement = true;
par.conditions.pos_moho = true;
% Match the active tests in Zach/Brennan a1 and make_mod_from_parms:
% fractional jump uses the mean of the velocities on the two sides.
par.conditions.max_moho_fraction = .25;
par.conditions.max_moho_jump = .8; % km/s; Inf disables the absolute limit
% These cross-zone constraints condition the joint prior. Consequently the
% marginal count prior need not remain exactly 1/K; see README equations.

%% =================== STARTING MODEL, INDEPENDENT OF SYNTHETIC TRUTH ===================
par.mod.start.method = 'reference'; % 'prior' gives dispersed prior draws
par.mod.start.zsed = .25;
par.mod.start.zmoh = 36;
par.mod.start.f = {[.10;.40],[.25;.65],.35};
par.mod.start.vs = {[.30;.80;1.70],[3.20;3.65;3.95],[4.25;4.45]};
par.mod.start.halfspace = 4.65;
par.mod.start.jitter = true;
par.mod.start.max_attempts = 500;

%% =================== HIERARCHICAL DATA ERRORS ===================
par.mod.noise.log_center = zeros(5,1);
par.mod.noise.log_width = 2*ones(5,1);
par.mod.noise.lower = log(.1)*ones(5,1);
par.mod.noise.upper = log(50)*ones(5,1);

%% =================== CPS AND HVf: SAME MODEL FOR ALL FIVE DATASETS ===================
par.forc.project_root = repo;
par.forc.shared_project = paths.shared_project;
par.forc.cps_bin = paths.cps_bin;
par.forc.hvf_bin = paths.hvf_bin;
par.forc.max_layer_thickness_km = 2;
par.forc.verbose = false;
par.forc.rf = struct('rayp_s_km',.06,'gaussian_alpha',2.5, ...
    'dt',.025,'npts',4096,'delay',10,'plot_window',[-2 25]);

%% =================== KNOWN THIN SYNTHETIC DATA ===================
par.synth.seed = 9041;
par.synth.gamma = ones(5,1);
par.synth.grids.hvsr.f = logspace(log10(.2),log10(25),80)';
par.synth.grids.ellip.T = linspace(3,25,50)';
par.synth.grids.cph.T = linspace(9,35,35)';
par.synth.grids.ugr.T = linspace(3,25,40)';
par.synth.rf_window = [-2 25];
par.synth.rf_stride = 8;
par.synth.rf_corr_seconds = .4;
par.synth.relative_sigma = [.05;.03;0;0;0];
par.synth.sigma_floor = [.05;.015;.03;.04;.015];

%% =================== OUTPUT AND COMMON DEPTH GRID ===================
par.res.root = paths.inversions;
par.res.scratch_root = paths.scratch;
par.res.show_figures = true;
par.res.nprior = 1000; % independent conditioned-prior draws; 0 skips this stage
par.res.prediction_stride = 10; % save cached responses every N retained states; 0 disables
par.res.depth_km = unique([linspace(0,.2,101),linspace(.2,2,91),linspace(2,85,250)])';
par.res.interface_bins_km = unique([linspace(0,.2,41),linspace(.2,2,37),linspace(2,70,137)]);

%% =================== NEXT TEST: 60-S DISPERSION, Ps/Sp, CORRECTED DR ===================
if any(strcmp(profile,{'extended60','extended100'}))
    par.inv.datatypes={'hvsr','ellip','cph','ugr','rf','rf_sp'};
    par.inv.nchains=2; par.inv.niter=10000; par.inv.burnin=4000;
    par.inv.checkpoint_every=100; par.inv.print_every=100;
    par.inv.max_consecutive_failures=1; % expose the first native failure in this debug run
    par.proposal.dr.enabled=true;
    par.proposal.dr.mode='layered_path'; % paper-style random second move FROM rejected Earth
    par.proposal.dr.same_zone_weight=0; % no preference for the first move's zone
    % Both stages use the ordinary move/zone mixture, including birth/death.
    % Keep the full reverse-path DR acceptance, including rejection factors.
    % dr.scale applies only to the optional fixed_scalar comparison kernel.
    par.proposal.birth_broad_probability=.30;
    par.proposal.interface_step=[.04 .004 .01];
    par.proposal.boundary_step=[.005 .5];
    % Keep shallow samples; append long periods rather than replacing them.
    par.synth.grids.ellip.T=unique([par.synth.grids.ellip.T;(26:2:60)']);
    par.synth.grids.cph.T=unique([par.synth.grids.cph.T;(36:2:60)']);
    par.synth.grids.ugr.T=unique([par.synth.grids.ugr.T;(26:2:60)']);
    % RF frequency passbands through 60 s, not merely 60-s lag windows.
    % A 409.6-s native record resolves >6 cycles at the longest period.
    par.forc.rf=struct('rayp_s_km',.06,'gaussian_alpha',2.5, ...
        'dt',.1,'npts',4096,'delay',140,'plot_window',[-100 120], ...
        'period_band_s',[1.25 60],'diagnostic_period_bands_s',[1.25 10;10 30;30 60]);
    par.forc.rf_sp=struct('rayp_s_km',.10,'gaussian_alpha',1.0, ...
        'component_gaussian_alpha',10,'water_level',.01,'edge_taper_fraction',.05, ...
        'dt',.1,'npts',4096,'delay',160,'plot_window',[-120 100], ...
        'period_band_s',[3 60],'diagnostic_period_bands_s',[3 10;10 30;30 60]);
    par.synth.rf_window=[-100 120]; par.synth.rf_stride=2;
    par.synth.rf_sp_window=[-120 100]; par.synth.rf_sp_stride=2;
    par.synth.rf_sp_corr_seconds=1;
    nd=numel(par.inv.datatypes);
    par.mod.noise.log_center=zeros(nd,1); par.mod.noise.log_width=2*ones(nd,1);
    par.mod.noise.lower=log(.1)*ones(nd,1); par.mod.noise.upper=log(50)*ones(nd,1);
    par.synth.gamma=ones(nd,1);
    par.synth.relative_sigma=[.05;.03;0;0;0;0];
    par.synth.sigma_floor=[.05;.015;.03;.04;.015;.015];
    % The old 70-km half-space cannot represent deeper variable structure.
    % Continue the same thin truth's 4.54-km/s half-space through 200 km;
    % a2 adds that constant-velocity segment, without inventing a deep jump.
    par.mod.maxz=200;
    par.res.depth_km=unique([par.res.depth_km(:);(86:2:220)']);
    par.res.interface_bins_km=unique([par.res.interface_bins_km(:);(71:1:200)'])';
end

%% =================== 250-KM THIN EARTH WITH PERIODS THROUGH 100 S ===================
if strcmp(profile,'extended100')
    par.mod.maxz=250;
    % HVSR has its own analysis band: 0.1-25 Hz, independent of the 100-s
    % dispersion/RF coverage. Preserve the existing 0.2-25-Hz coordinates.
    par.synth.grids.hvsr.f=unique([logspace(log10(.1),log10(.2),12)';par.synth.grids.hvsr.f]);
    for name={'ellip','cph','ugr'}
        dt=name{1}; par.synth.grids.(dt).T=unique([par.synth.grids.(dt).T;(62:2:100)']);
    end
    for name={'rf','rf_sp'}
        dt=name{1}; cfg=par.forc.(dt);
        cfg.period_band_s(2)=100;
        cfg.diagnostic_period_bands_s=[cfg.diagnostic_period_bands_s;60 100];
        % CPS hrftn96 -2 uses an internal double-length FFT. Keep the
        % supported 4096 output samples; dt=.2 gives an 819.2-s record.
        % Retain every .2-s sample: same likelihood cadence as extended60.
        cfg.dt=.2; cfg.delay=300;
        par.forc.(dt)=cfg;
        par.synth.([dt '_stride'])=1;
    end
    % Conversion windows already include 250-km Ps/Sp arrivals. They are
    % lag windows, not passbands; the extra native padding supports 100 s.
    par.res.depth_km=unique([par.res.depth_km(:);(222:2:270)']);
    par.res.interface_bins_km=unique([par.res.interface_bins_km(:);(201:250)'])';
end
end
