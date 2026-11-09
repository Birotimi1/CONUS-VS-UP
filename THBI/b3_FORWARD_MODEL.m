function [predata, fwd] = b3_FORWARD_MODEL(model, data, par, workdir)
% Forward calculate the five observables from one variable-layer Earth.
% The preparation, HVf, CPS dispersion, RF and output blocks below execute
% directly in this main function, as in THBI's b3 forward-model stages.
% Units: h km, Vp/Vs km/s, density g/cm^3, frequency Hz, period/time s.
% Layer/IO utilities are in the repository functions/ folder; calculations are here.

%% =================== REQUESTED DATA GRIDS ===================
grids.hvsr.f = data.hvsr.x;
grids.ellip.T = data.ellip.x;
grids.cph.T = data.cph.x;
grids.ugr.T = data.ugr.x;
cfg = par.forc;
output_dir = workdir;
physical_model=model;
stack_mode=isfield(par,'rf') && strcmp(par.rf.method,'hk_ccp');
model = model.matrix;
isvti=size(model,2)==7;
rf_program='hrftn96'; hv_program='HVf'; hv_flags='';
if isvti, rf_program='trftn96'; hv_program='HVf_VTI'; hv_flags=' -vti'; end
rf_sac=[rf_program '.sac'];

%% =================== SOLVER CONFIGURATION AND PATHS ===================
% Evaluate one layered Earth for all five observables (no likelihood here).
% grids.hvsr.f (Hz); grids.ellip.T/cph.T/ugr.T (s). RF grid is cfg.rf.
% Model columns: thickness km, Vp km/s, Vs km/s, density g/cm^3.
validate_layer_model(model);
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
% Absolute solver paths are necessary because each executable changes cwd.
[ok, attrs] = fileattrib(output_dir);
if ~ok, error('fiveobs:io', 'Cannot access output directory.'); end
output_dir = attrs.Name;
previous_directory = pwd;
restore_directory = onCleanup(@() cd(previous_directory)); %#ok<NASGU>
cfg = merge_defaults(cfg, struct('cps_bin', '', 'hvf_bin', '', ...
    'project_root', fileparts(fileparts(mfilename('fullpath'))), ...
    'shared_project', '', 'max_layer_thickness_km', 2, 'verbose', true));
if ~isfield(cfg, 'hvsr'), cfg.hvsr = struct(); end
cfg.hvsr = merge_defaults(cfg.hvsr, struct('nmr', 10, 'nml', 10, 'nks', 1000, ...
    'prec', 1e-4, 'apsv', 1e-3, 'ash', 1e-3));
if ~isfield(cfg, 'rf'), cfg.rf = struct(); end
cfg.rf = merge_defaults(cfg.rf, struct('rayp_s_km', 0.06, 'gaussian_alpha', 2.5, ...
    'dt', 0.025, 'npts', 4096, 'delay', 10, 'plot_window', [-2 25]));
cfg.cps_bin = find_bin(cfg.cps_bin, 'HBI_CPS_BIN', ...
    {fullfile(cfg.project_root, 'AGU_models', 'bin_v3.30'), ...
     fullfile(cfg.shared_project, 'AGU_models', 'bin_v3.30')}, 'surf96');
if isvti
    if ~isfield(cfg,'hvf_isotropic_bin'), cfg.hvf_isotropic_bin=cfg.hvf_bin; end
    cfg.hvf_bin=cfg.hvf_vti_bin;
    assert(exist(fullfile(cfg.hvf_bin,hv_program),'file')==2,'fiveobs:dependency', ...
        'Build native/HVf_VTI/build.sh before the VTI preflight.');
    required={'tprep96','tdisp96','tregn96','tdpegn96',rf_program};
else
    cfg.hvf_bin = find_bin(cfg.hvf_bin, 'HBI_HVF_BIN', ...
        {fullfile(cfg.project_root,'AGU_models','bin_HV_DFA'), ...
         fullfile(cfg.shared_project,'AGU_models','bin_HV_DFA')}, 'HVf');
    required={rf_program};
end
for j=1:numel(required)
    assert(exist(fullfile(cfg.cps_bin,required{j}),'file')==2, ...
        'fiveobs:dependency','Missing CPS program %s.',required{j});
end
old_path = getenv('PATH');
restore_path = onCleanup(@() setenv('PATH', old_path)); %#ok<NASGU>
setenv('PATH', [cfg.cps_bin ':' cfg.hvf_bin ':' old_path]);
%% Preserve physical interfaces while subdividing thick layers
source_model = model;
model = refine_layer_model(model, cfg.max_layer_thickness_km);
forward = struct('source_model', source_model, 'model', model, ...
    'config', cfg, 'output_dir', output_dir, 'pred', struct());
forward.geometry = 'Flat elastic; identical layers and tensor for all five observables';
names = {'hvsr', 'ellip', 'cph', 'ugr'};
coords = {'f', 'T', 'T', 'T'};
% Validate all requested grids before invoking any executable.
for k = 1:4
    dt = names{k}; coord = coords{k};
    if ~isfield(grids, dt) || ~isfield(grids.(dt), coord)
        error('fiveobs:grid', 'Missing grids.%s.%s.', dt, coord);
    end
    x = grids.(dt).(coord);
    if ~isnumeric(x) || ~isreal(x) || ~isvector(x) || isempty(x) || any(~isfinite(x) | x <= 0)
        error('fiveobs:grid', '%s coordinates must be finite positive vectors.', dt);
    end
    if any(strcmp(dt, {'cph', 'ugr'})) && any(round(x * 1e4) <= 0)
        error('fiveobs:grid', 'Periods must remain positive at CPS exchange precision.');
    end
end
%% Calculate HVSR, ellipticity, phase velocity, and group velocity
started = tic;
solver_failures={};
for k = 1:4
    dt = names{k}; x = grids.(dt).(coords{k})(:);
    if cfg.verbose
        fprintf('Calculating %s from the shared %d-layer model...\n', dt, size(model, 1));
    end
    workdir = fullfile(output_dir, dt);
    if ~exist(workdir, 'dir'), mkdir(workdir); end
    try
        if any(strcmp(dt, {'hvsr', 'ellip'}))
            cd(workdir);
            fid = fopen('model.txt', 'w');
            if fid < 0, error('fiveobs:io', 'Cannot write HV-DFA model.'); end
            fprintf(fid, '%d\n', size(model, 1));
            exchange=model; exchange(:,1:4)=1000*exchange(:,1:4);
            format=[repmat('%.12g ',1,size(exchange,2)) '\n'];
            fprintf(fid,format,exchange');
            fclose(fid);
            %% HVSR: surface-wave modes plus the diffuse body-wave contribution
            if strcmp(dt, 'hvsr')
                frequencies = x;
                nmr = cfg.hvsr.nmr; nml = cfg.hvsr.nml; nks = cfg.hvsr.nks;
                apsv = cfg.hvsr.apsv; ash = cfg.hvsr.ash;
            else
                % ELLIPTICITY: fundamental Rayleigh H/V only
                frequencies = 1 ./ x;
                nmr = 1; nml = 0; nks = 0; apsv = 0; ash = 0;
            end
            [requested, ~, map] = unique(frequencies);
            fid = fopen('frequencies.txt', 'w');
            fprintf(fid, '%.15g\n', requested); fclose(fid);
            cmd = sprintf('%s%s -f model.txt -ff frequencies.txt -nmr %d -nml %d -nks %d -prec %.10g -apsv %.10g -ash %.10g -rep -hv', ...
                quote_shell_path(fullfile(cfg.hvf_bin,hv_program)), hv_flags, nmr, nml, nks, cfg.hvsr.prec, apsv, ash);
            [status, message] = system([cmd ' > HV.dat 2> HVf.log']);
            if status ~= 0, error('fiveobs:hvf', 'HVf failed (%d): %s. See %s', status, message, workdir); end
            fid = fopen('HV.dat', 'r');
            if fid < 0, error('fiveobs:hvf', 'HVf did not produce output.'); end
            values = fscanf(fid, '%f'); fclose(fid);
            if numel(values) ~= 2 * numel(requested)
                error('fiveobs:hvf', 'HVf returned an incomplete frequency grid; see %s', fullfile(workdir, 'HV.dat'));
            end
            values = reshape(values, 2, [])';
            if any(abs(values(:, 1) - requested) > 1e-6 * max(1, abs(requested)))
                error('fiveobs:hvf', 'HVf returned a different frequency grid.');
            end
            prediction = struct('val', values(map, 2), 'model_input', model, 'command', cmd, 'status', 'computed');
        elseif isvti
            cd(workdir);
            % CPS TI surface-wave sequence, all actual calculation blocks here.
            % Remove only known generated caches before each candidate Earth.
            cache_names={'tdisp96.dat','tdisp96.ray','tregn96.egn','TREGN.dsp','TREGN.ASC'};
            for ic=1:numel(cache_names)
                if exist(cache_names{ic},'file'), delete(cache_names{ic}); end
            end
            write_elastic_cps_model('start.mod',model);
            [requested,~,map]=unique(x);
            fid=fopen('periods.txt','w');
            if fid<0, error('fiveobs:io','Cannot write CPS periods.'); end
            fprintf(fid,'%.15g\n',requested); fclose(fid);
            commands={sprintf('%s -M start.mod -PARR periods.txt -R -NMOD 1', ...
                quote_shell_path(fullfile(cfg.cps_bin,'tprep96'))), ...
                quote_shell_path(fullfile(cfg.cps_bin,'tdisp96')), ...
                quote_shell_path(fullfile(cfg.cps_bin,'tregn96')), ...
                [quote_shell_path(fullfile(cfg.cps_bin,'tdpegn96')) ' -R -PER -S -ASC']};
            for ic=1:numel(commands)
                [status,message]=system([commands{ic} ' > cps_ti_' num2str(ic) '.log 2>&1']);
                if status~=0, error('fiveobs:cps','CPS TI stage %d failed: %s',ic,message); end
            end
            if ~exist('TREGN.dsp','file'), error('fiveobs:cps','CPS TI output is missing.'); end
            fid=fopen('TREGN.dsp','r');
            columns=textscan(fid,'%s %s %s %s %f %f %f %f'); fclose(fid);
            velocity_type='C'; if strcmp(dt,'ugr'), velocity_type='U'; end
            selected=strcmp(columns{2},'R') & strcmp(columns{3},velocity_type) & columns{5}==0;
            periods=columns{6}(selected); velocities=columns{7}(selected);
            predicted=NaN(size(requested));
            % tdpegn96 uses g14.5: tolerance follows its five significant digits.
            for ip=1:numel(requested)
                rows=find(abs(periods-requested(ip))<=6e-5*max(1,requested(ip)));
                if numel(rows)==1, predicted(ip)=velocities(rows); end
            end
            prediction=struct('val',predicted(map),'model_input',model, ...
                'command',{commands},'status','computed');
        else
            cd(workdir);
            % surf96 reuses tmpsrfi.* (including the old model) even after start.mod is
            % replaced. Equivalent to `surf96 39`: clear only CPS's generated cache.
            % Essential when MCMC evaluates successive models in the same directory.
            cache = dir('tmpsrfi.*');
            for icache = 1:numel(cache)
                if ~cache(icache).isdir, delete(cache(icache).name); end
            end
            write_elastic_cps_model('start.mod', model);
            % PHASE ('C') and GROUP ('U'): fundamental Rayleigh mode
            if strcmp(dt, 'cph'), velocity_type = 'C'; else, velocity_type = 'U'; end
            % CPS exchange files print periods to four decimals. Match by returned period,
            % never by position or by truncating the requested observations.
            [requested, ~, map] = unique(round(x * 1e4) / 1e4);
            fid = fopen('disp_obs.dsp', 'w');
            for iperiod = 1:numel(requested)
                fprintf(fid, 'SURF96 R %s X 0 %.4f 3.0 0.1\n', velocity_type, requested(iperiod));
            end
            fclose(fid);
            fid = fopen('sobs.d', 'w');
            fprintf(fid, '0.005 0.005 0.0 0.005 0.0\n0 0 0 0 0 1 1 0 1 0\nstart.mod\ndisp_obs.dsp\n');
            fclose(fid);
            if exist('predicted.dsp', 'file'), delete('predicted.dsp'); end
            cmd = [quote_shell_path(fullfile(cfg.cps_bin, 'surf96')) ' 1 27 predicted.dsp'];
            [status, message] = system([cmd ' > surf96.log 2>&1']);
            if status ~= 0 || ~exist('predicted.dsp', 'file')
                error('fiveobs:cps', 'surf96 failed (%d): %s. See %s', status, message, fullfile(workdir, 'surf96.log'));
            end
            fid = fopen('predicted.dsp', 'r');
            columns = textscan(fid, '%s %s %s %s %f %f %f %f'); fclose(fid);
            if isempty(columns{6}) || any(~strcmp(columns{2}, 'R')) || ...
                    any(~strcmp(columns{3}, velocity_type)) || any(columns{5} ~= 0)
                error('fiveobs:cps', 'Unexpected CPS wave/mode/velocity type or empty output.');
            end
            predicted = NaN(size(requested));
            for iperiod = 1:numel(requested)
                rows = find(abs(columns{6} - requested(iperiod)) < 6e-5);
                if numel(rows) == 1, predicted(iperiod) = columns{7}(rows); end
            end
            prediction = struct('val', predicted(map), 'model_input', model, 'command', cmd, 'status', 'computed');
        end
        if numel(prediction.val) ~= numel(x) || any(~isfinite(prediction.val) | prediction.val <= 0)
            error('fiveobs:incomplete', '%s has missing or nonphysical predictions; see %s.', dt, workdir);
        end
        prediction.status = 'ok';
    catch problem
        % Only known native-output failures can become rejected proposals.
        % Programming errors, filesystem errors and missing dependencies stop here.
        if ~any(strcmp(problem.identifier,{'fiveobs:hvf','fiveobs:cps','fiveobs:incomplete'}))
            rethrow(problem);
        end
        solver_failures{end+1}=problem; %#ok<AGROW>
        prediction = struct('val', NaN(size(x)), 'status', 'failed', ...
            'message', problem.message, 'model_input', model);
        warning('fiveobs:solver', '%s: %s', dt, problem.message);
    end
    prediction.x = x;
    forward.pred.(dt) = prediction;
end
%% Calculate the incident-P receiver function
if stack_mode && isfield(data,'rf') && isfield(data.rf,'waves')
    % Observed Ps waves are fixed; recompute their H-kappa stack using this
    % candidate's harmonic crustal Vs, and evaluate its H and sampled kappa.
    hk=b3_FORWARD_MODEL_HK(physical_model,data.rf,par);
    forward.pred.rf=struct('x',1,'val',hk.E_by_Emax,'hk',hk,'status','ok', ...
        'model_input',model,'kind','Ps H-kappa stack energy ratio','observed',false);
else
if cfg.verbose, fprintf('Calculating synthetic P-to-S receiver function...\n'); end
try
    rf_cfg = cfg.rf;
    rf_workdir = fullfile(output_dir,'rf');
    % Synthetic P receiver function from the shared layered model.
    % hrftn96 computes the radial/vertical frequency-domain response directly,
    % including multiples, then applies its Gaussian filter. No observed RF,
    % arbitrary amplitude normalization, or fabricated RF misfit is used.
    if ~isscalar(rf_cfg.rayp_s_km) || ~isfinite(rf_cfg.rayp_s_km) || rf_cfg.rayp_s_km <= 0
        error('fiveobs:rayp', 'Require 0 < rayp < 1/Vp in every layer (rayp in s/km).');
    end
    if ~isscalar(rf_cfg.npts) || rf_cfg.npts < 64 || rf_cfg.npts > 4096 ...
            || log2(rf_cfg.npts) ~= floor(log2(rf_cfg.npts))
        error('fiveobs:rf_config', 'npts must be a power of two, 64..4096; -2 uses twice this length.');
    end
    if ~isscalar(rf_cfg.dt) || ~isfinite(rf_cfg.dt) || rf_cfg.dt <= 0 || ...
            ~isscalar(rf_cfg.delay) || ~isfinite(rf_cfg.delay) || rf_cfg.delay <= 0 || ...
            rf_cfg.delay >= rf_cfg.npts * rf_cfg.dt || ...
            ~isscalar(rf_cfg.gaussian_alpha) || ~isfinite(rf_cfg.gaussian_alpha) || rf_cfg.gaussian_alpha <= 0
        error('fiveobs:rf_config', 'Invalid RF sample interval, delay, or Gaussian alpha.');
    end
    if ~exist(rf_workdir, 'dir'), mkdir(rf_workdir); end
    cd(rf_workdir);
    vti_vertical_slowness(model,rf_cfg.rayp_s_km);
    write_elastic_cps_model('rf_model.mod', model);
    cmd = sprintf('%s -P -2 -M rf_model.mod -RAYP %.10g -ALP %.10g -DT %.10g -NSAMP %d -D %.10g', ...
        quote_shell_path(fullfile(cfg.cps_bin,rf_program)), rf_cfg.rayp_s_km, ...
        rf_cfg.gaussian_alpha, rf_cfg.dt, rf_cfg.npts, rf_cfg.delay);
    if exist(rf_sac, 'file'), delete(rf_sac); end
    [status, message] = system([cmd ' > hrftn96.log 2>&1']);
    if status ~= 0 || ~exist(rf_sac, 'file')
        error('fiveobs:rf', 'hrftn96 failed (%d): %s. See %s', status, message, fullfile(rf_workdir, 'hrftn96.log'));
    end
    rf = read_cps_sac(rf_sac);
    if numel(rf.val) ~= rf_cfg.npts || abs(rf.dt - rf_cfg.dt) > 1e-6 || ...
            abs(rf.t(1) + rf_cfg.delay) > 1e-5 || ...
            abs(rf.rayp_s_km - rf_cfg.rayp_s_km) > 1e-6 || ...
            abs(rf.gaussian_alpha - rf_cfg.gaussian_alpha) > 1e-5
        error('fiveobs:rf', 'RF SAC header does not match the requested settings.');
    end
    rf.kind = 'Synthetic P-to-S RF: Gaussian-filtered radial/vertical transfer response';
    rf.model_input = model;
    rf.command = cmd;
    rf.observed = false;
    rf.status = 'ok';
    % Kinematic primary Ps delays to each interface, not full waveform picks.
    [qp,qs]=vti_vertical_slowness(model,rf_cfg.rayp_s_km);
    qp=qp(1:end-1); qs=qs(1:end-1);
    rf.interface_depth_km = cumsum(model(1:end-1, 1));
    rf.primary_ps_delay_s = cumsum(model(1:end-1, 1) .* (qs - qp));
    rf.interface_contrast = any(abs(diff(model(:, 2:end), 1, 1)) > 1e-10, 2);
    rf=receiver_function_bands(rf,rf_cfg);
    forward.pred.rf = rf;
catch problem
    if ~any(strcmp(problem.identifier,{'fiveobs:rf','fiveobs:sac'}))
        rethrow(problem);
    end
    solver_failures{end+1}=problem;
    forward.pred.rf = struct('t', [], 'val', [], 'status', 'failed', ...
        'message', problem.message, 'model_input', model, 'observed', false);
    warning('fiveobs:solver', 'Receiver function: %s', problem.message);
end
end

%% =================== Sp: CPS SV COMPONENTS, THEN REGULARIZED Z/R ===================
% CPS documents that its direct -S ratio can be unstable. Obtain the two
% physical waveforms with -S -r/-z, then deconvolve on the real frequency axis.
% Source/sign convention: CPS hrftn96.f (Uz positive upward). Sp precursors
% keep their negative lag; do not reverse time or multiply by an arbitrary sign.
if any(strcmp(par.inv.datatypes,'rf_sp'))
    rf_cfg=cfg.rf_sp; rf_workdir=fullfile(output_dir,'rf_sp');
    if ~exist(rf_workdir,'dir'), mkdir(rf_workdir); end
    try
        assert(isfinite(rf_cfg.rayp_s_km) && rf_cfg.rayp_s_km>0 && ...
            rf_cfg.rayp_s_km>0, ...
            'fiveobs:rf_config','Sp requires propagating P and S in every layer; no silent model truncation.');
        assert(rf_cfg.npts>=64 && rf_cfg.npts<=4096 && mod(log2(rf_cfg.npts),1)==0 && ...
            isfinite(rf_cfg.dt) && rf_cfg.dt>0 && rf_cfg.delay>0 && ...
            rf_cfg.delay<rf_cfg.npts*rf_cfg.dt && rf_cfg.gaussian_alpha>0 && ...
            rf_cfg.component_gaussian_alpha>=rf_cfg.gaussian_alpha, ...
            'fiveobs:rf_config','Invalid Sp sample/filter settings.');
        vti_vertical_slowness(model,rf_cfg.rayp_s_km);
        cd(rf_workdir); write_elastic_cps_model('rf_model.mod',model);
        components=cell(2,1); commands=cell(2,1); flags={'r','z'};
        for component=1:2
            flag=flags{component};
            cmd=sprintf('%s -S -%s -2 -M rf_model.mod -RAYP %.10g -ALP %.10g -DT %.10g -NSAMP %d -D %.10g', ...
                quote_shell_path(fullfile(cfg.cps_bin,rf_program)),flag,rf_cfg.rayp_s_km, ...
                rf_cfg.component_gaussian_alpha,rf_cfg.dt,rf_cfg.npts,rf_cfg.delay);
            if exist(rf_sac,'file'), delete(rf_sac); end
            [status,message]=system([cmd ' > hrftn96_' flag '.log 2>&1']);
            if status~=0 || ~exist(rf_sac,'file')
                error('fiveobs:rf','Sp component %s failed (%d): %s',flag,status,message);
            end
            trace=read_cps_sac(rf_sac);
            assert(numel(trace.val)==rf_cfg.npts && abs(trace.dt-rf_cfg.dt)<1e-6 && ...
                abs(trace.t(1)+rf_cfg.delay)<1e-5 && abs(trace.rayp_s_km-rf_cfg.rayp_s_km)<1e-6 && ...
                abs(trace.gaussian_alpha-rf_cfg.component_gaussian_alpha)<1e-5, ...
                'fiveobs:rf','Sp component SAC headers disagree with the request.');
            components{component}=trace; commands{component}=cmd;
            movefile(rf_sac,[upper(flag) '.sac'],'f');
        end
        if stack_mode
            rf=b3_FORWARD_MODEL_RF_CCP(model,components,par);
        else
            rf=waterlevel_receiver_function(components{1},components{2},rf_cfg);
            rf=receiver_function_bands(rf,rf_cfg);
        end
        rf.model_input=model; rf.command=commands; rf.observed=false; rf.status='ok';
        [qp,qs]=vti_vertical_slowness(model,rf_cfg.rayp_s_km);
        qp=qp(1:end-1); qs=qs(1:end-1);
        rf.interface_depth_km=cumsum(model(1:end-1,1));
        rf.primary_sp_delay_s=-cumsum(model(1:end-1,1).*(qs-qp));
        rf.interface_contrast=any(abs(diff(model(:,2:end),1,1))>1e-10,2);
        forward.pred.rf_sp=rf;
    catch problem
        if ~any(strcmp(problem.identifier,{'fiveobs:rf','fiveobs:sac'})), rethrow(problem); end
        solver_failures{end+1}=problem;
        forward.pred.rf_sp=struct('t',[],'val',[],'status','failed', ...
            'message',problem.message,'model_input',model,'observed',false);
        warning('fiveobs:solver','Sp receiver function: %s',problem.message);
    end
end
forward.all_forward_succeeded = strcmp(forward.pred.rf.status, 'ok');
if isfield(forward.pred,'rf_sp')
    forward.all_forward_succeeded=forward.all_forward_succeeded && strcmp(forward.pred.rf_sp.status,'ok');
end
for k = 1:4
    forward.all_forward_succeeded = forward.all_forward_succeeded && strcmp(forward.pred.(names{k}).status, 'ok');
end
forward.elapsed_seconds = toc(started);

fwd = forward;
if ~fwd.all_forward_succeeded
    problem=MException('thbi:forward_failure','Incomplete native forward calculation; inspect %s.',output_dir);
    for k=1:numel(solver_failures)
        problem=addCause(problem,solver_failures{k});
    end
    throw(problem);
end
%% =================== PUT PREDICTIONS ON THE DATA GRIDS ===================
for k = 1:numel(par.inv.datatypes)
    d = par.inv.datatypes{k}; p = fwd.pred.(d);
    assert(isequal(p.model_input,fwd.model),'Forward models used different Earths.');
    if stack_mode && strcmp(d,'rf_sp')
        if isfield(data,d)
            assert(isequal(p.x,data.(d).x),'thbi:ccp','CCP coordinates changed.');
        end
        val=p.val; predata.rf_sp_parent=p.parent;
    elseif stack_mode && strcmp(d,'rf') && isfield(p,'hk')
        val=p.val; predata.hk=p.hk;
    elseif any(strcmp(d,{'rf','rf_sp'})) && ~isfield(data,d)
        val = p.val;
    elseif any(strcmp(d,{'rf','rf_sp'}))
        ix = round((data.(d).x-p.t(1))/p.dt)+1;
        if any(ix<1 | ix>numel(p.t)) || any(abs(p.t(ix)-data.(d).x)>2e-5)
            error('hbi:forward','RF coordinates do not match native CPS samples.');
        end
        val = p.val(ix);
    else
        if numel(p.x)~=numel(data.(d).x) || any(abs(p.x-data.(d).x)>1e-8)
            error('hbi:forward','Prediction coordinates changed for %s.',d);
        end
        val = p.val;
    end
    predata.(d) = val(:);
end
end
