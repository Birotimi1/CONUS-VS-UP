% CHECK_THBI_EXTENDED60 -- run this script yourself in MATLAB BEFORE inversion.
% Known thin Earth -> six noisy data blocks (HVSR, ellipticity, phase, group, Ps and Sp) -> zoned transdimensional HBI.
% Algebraic/toy-chain checks, native RF controls, thin true Earth and split/merge.
% No production inversion. All native outputs and band figures are preserved.

%% =================== a0: PROJECT AND SOLVER PATHS ===================
launcher=mfilename('fullpath');
repo=fileparts(fileparts(launcher));
addpath(fullfile(repo,'THBI'));
[repo,paths]=a0_STARTUP_BAYES();

%% =================== EXPLICIT TEST SETTINGS ===================
par=bayes_inv_parms('extended60');
par.inv.profile='extended60_preflight';
par.inv.nchains=2;
par.inv.niter=10000;
par.inv.burnin=4000;
par.inv.parallel=false;          % two serial chains; no Parallel Toolbox needed
par.inv.seed=42000;              % sampler uses seed + 1009*chain_number
par.inv.checkpoint_every=100;
par.inv.print_every=100;
par.inv.prior_only=false;        % real CPS/HVf likelihoods for all six data blocks
par.res.prediction_stride=10;    % 600 cached response draws per chain; no state thinning
par.res.nprior=1000;
par.res.show_figures=true;
run_implementation_checks=true; % short algebraic/toy checks BEFORE native work
validate_configuration(par);
assert(strcmp(par.proposal.dr.mode,'layered_path') && par.proposal.dr.same_zone_weight==0, ...
    'thbi:test','This launcher follows the paper: random second proposal from the rejected Earth.');
assert(all(par.mod.kmax>par.mod.kmin) && all(par.proposal.probability(4:5)>0), ...
    'thbi:test','This test requires variable layer counts and birth/death moves.');

%% =================== NEW OUTPUT DIRECTORY AND LAUNCHER RECORD ===================
if ~exist(par.res.root,'dir'), mkdir(par.res.root); end
result_dir=fullfile(par.res.root,['ECM1_thin_THBI_' par.inv.profile '_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(result_dir,'dir'),'thbi:test','Result directory already exists; wait a second and rerun.');
mkdir(result_dir); par.res.dir=result_dir;
copyfile([launcher '.m'],fullfile(result_dir,'launcher.m'));
runtime=struct('version',version,'computer',computer,'started',datestr(now,30));
save(fullfile(result_dir,'test_settings.mat'),'par','runtime','-v7');
fprintf('\nOutput: %s\n',result_dir);
fprintf('Forward-only preflight: no production MCMC.\n');

%% =================== IMPLEMENTATION CHECKS BEFORE EXPENSIVE SAMPLING ===================
% The failure check deliberately raises AND catches errors. Use dbstop if
% error for debugging; dbstop if caught error also stops at these expected errors.
if run_implementation_checks
    test_thbi_rj();
    test_thbi_extended60();
    test_thbi_layered_dr();
    test_thbi_dr_stationarity();
    test_thbi_chain();
    test_thbi_chain('extended60');
    test_thbi_chain('extended60','fixed_scalar');
    test_thbi_failures();
    test_thbi_failures('extended60');
    test_thbi_results();
end

%% =================== INDEPENDENT NATIVE Ps/Sp CONTROLS ===================
rf_checks=test_thbi_receiver_functions(par,fullfile(result_dir,'native_rf_checks'));

%% =================== a2: KNOWN THIN EARTH AND SIX SYNTHETIC DATA BLOCKS ===================
[data,TRUEmodel,forward]=a2_LOAD_DATA(par);
par.forc=forward.config;
plot_RF_BANDS(forward,par,'true');
save(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','forward','runtime','-v7');
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; v=data.(name);
    assert(numel(v.clean)==v.N && all(isfinite(v.clean)) && all(isfinite(v.obs)), ...
        'thbi:test','Invalid true/noisy synthetic data for %s.',name);
    fprintf('Synthetic %-5s: %d finite samples.\n',name,v.N);
end

%% =================== NATIVE SPLIT/MERGE CHECK ON THE ACTUAL TEST GRIDS ===================
% Exercise different physical layer counts before committing to 20,000 steps.
% Reuse one solver directory to check that CPS/HVf do not reuse stale models.
check_par=par; check_par.mod.start.jitter=false;
check_model=b1_INITIATE_MODEL(check_par);
[born,birth]=b2_PERTURB_MODEL(check_model,par,struct( ...
    'type','birth','zone',2,'index',1,'split_f',.12,'delta_vs',.08));
assert(birth.valid,'thbi:test','Native preflight birth is outside prior support.');
[back,death]=b2_PERTURB_MODEL(born,par,struct('type','death','zone',2,'index',1));
assert(death.valid && isequal(back.K,check_model.K) && ...
    max(abs(back.matrix(:)-check_model.matrix(:)))<1e-10, ...
    'thbi:test','Native preflight merge did not restore the physical Earth.');
preflight_dir=fullfile(result_dir,'native_split_merge'); mkdir(preflight_dir);
check_models={check_model,born,back}; check_predictions=cell(3,1);
for j=1:3
    try
        [check_predictions{j},check_forward]=b3_FORWARD_MODEL(check_models{j},data,par,preflight_dir);
        for d=1:numel(par.inv.datatypes)
            name=par.inv.datatypes{d}; y=check_predictions{j}.(name);
            assert(isreal(y) && numel(y)==data.(name).N && all(isfinite(y)), ...
                'thbi:test','Invalid native preflight response for %s.',name);
            assert(isequal(check_forward.pred.(name).model_input,check_forward.model), ...
                'thbi:test','Native solvers used different Earths for %s.',name);
        end
    catch problem
        context=struct('stage','native split/merge preflight','iteration',0,'model_index',j);
        preserve_forward_failure(fullfile(result_dir,'preflight_failure'), ...
            check_models{j},data,par,problem,preflight_dir,context);
        rethrow(problem);
    end
end
native_check=struct();
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d};
    native_check.(name).split_change=max(abs(check_predictions{1}.(name)-check_predictions{2}.(name)));
    native_check.(name).roundtrip_error=max(abs(check_predictions{1}.(name)-check_predictions{3}.(name)));
end
save(fullfile(preflight_dir,'preflight.mat'),'check_models','check_predictions','native_check','-v7');
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d};
    assert(native_check.(name).split_change>1e-8 && native_check.(name).roundtrip_error<1e-6, ...
        'thbi:test','Native split/merge response check failed for %s; inspect preflight.mat.',name);
end
fprintf('PASS: all six native responses change after a split and return after the inverse merge.\n');

%% =================== SAVE THE PREFLIGHT OUTCOME ===================
preflight=struct('status','All programmed preflight checks passed', ...
    'receiver_functions',rf_checks,'split_merge',native_check, ...
    'note','Inspect the truth and RF band figures; these checks do not establish posterior recovery.');
save(fullfile(result_dir,'preflight_report.mat'),'preflight','-v7');
% Plot all curves and the actual true velocity structure, before inversion.
f=figure('Color','w','Position',[100 80 1400 900]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
nexttile; plot_layer_model(TRUEmodel.model_matrix,'k',1.7,'True thin Earth');
set(gca,'YDir','reverse'); ylim([0 par.mod.maxz]); grid on;
xlabel('Vs (km/s)'); ylabel('Depth (km)'); title('Known thin true Earth');
labels={'HVSR','Rayleigh ellipticity','Phase velocity','Group velocity','Ps R/Z','Sp Z/R'};
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; v=data.(name); nexttile;
    plot(v.x,v.obs,'.','Color',[.65 .65 .65]); hold on;
    plot(v.x,v.clean,'k-','LineWidth',1.4); grid on; title(labels{d});
    if d==1
        set(gca,'XScale','log'); xlabel('Frequency (Hz)');
    elseif d<5
        xlabel('Period (s)');
    else
        xlabel('Conversion lag (s)');
    end
    if d==3 || d==4, ylabel('km/s'); else, ylabel('Amplitude ratio'); end
end
sgtitle('Thin true Earth: extended periods, signed Ps and Sp');
exportgraphics(f,fullfile(result_dir,'true_observables.png'),'Resolution',180);
exportgraphics(f,fullfile(result_dir,'true_observables.pdf'),'ContentType','vector');
fprintf('\nPreflight complete. Inspect outputs in:\n%s\n',result_dir);
fprintf('Next run script: DATA/run_ecm1_thin_thbi_60s_ps_sp.m\n');
