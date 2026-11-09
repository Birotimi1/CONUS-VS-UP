% PLOT_ECM1_THIN_THBI_COLORED_CHAINS -- run manually in MATLAB.
% Adapted layout/colors from our Figure_codes/plot_model_unified_colored_chains.m
% and plot_model_unified_all_2.m. Read saved THBI results; no MCMC or solvers.
% Every physical model has its OWN K, interfaces and Vs. No layer-slot averaging.
% Accepted proposals form the ensemble; ALL retained states enter summaries.

%% =================== SELECT RUN AND DISPLAY SETTINGS ===================
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'THBI')); a0_STARTUP_BAYES();
result_dir='';                 % absolute run folder; empty reads LATEST_RUN.txt
show_figures=true;
max_models_per_chain=Inf;      % Inf shows EVERY accepted post-burn-in model
include_noise_accepts=true;    % noise-only acceptances repeat the same Vs model
plot_saved_responses=true;    % cached predictions only; NEVER forward calculate
max_responses_per_chain=150;  % deterministic display subset, summaries use all saved draws
save_matlab_figures=false;    % PNG/PDF always saved; .fig optional
ensemble_strength=.25;        % RGB blend against white, not unsupported RGBA line colors
chain_clrs=[0 .4 .8; .9 .1 .1; 0 .7 .3; .8 .5 0; .5 0 .8; 0 .7 .7];
assert(isscalar(max_models_per_chain) && max_models_per_chain>=1 && ...
    (isinf(max_models_per_chain) || max_models_per_chain==fix(max_models_per_chain)), ...
    'thbi:plot','Set a positive integer display limit or Inf.');
assert(isscalar(max_responses_per_chain) && isfinite(max_responses_per_chain) && ...
    max_responses_per_chain>=1 && max_responses_per_chain==fix(max_responses_per_chain), ...
    'thbi:plot','Response display limit must be a positive integer.');
assert(isscalar(ensemble_strength) && ensemble_strength>0 && ensemble_strength<=1, ...
    'thbi:plot','ensemble_strength must be in (0,1].');

%% =================== LOAD COMPLETED CHAINS, TRUTH AND THEIR ACTUAL SETTINGS ===================
if isempty(result_dir)
    pointer=fullfile(repo,'DATA','transdimensional','LATEST_RUN.txt');
    assert(exist(pointer,'file')==2,'thbi:plot','Set result_dir to the completed synthetic run.');
    result_dir=strtrim(fileread(pointer));
end
assert(exist(fullfile(result_dir,'results.mat'),'file')==2, ...
    'thbi:plot','Missing results.mat in %s. Finish/resume the run first.',result_dir);
saved=load(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary');
par=saved.par; data=saved.data; TRUEmodel=saved.TRUEmodel; chains=saved.chains;
assert(iscell(chains) && numel(chains)==par.inv.nchains && ~isempty(chains), ...
    'thbi:plot','Saved chain count disagrees with the run settings.');
nc=numel(chains); keep=(par.inv.burnin+1:par.inv.niter)'; n=numel(keep);
assert(n>0,'thbi:plot','No post-burn-in states.');
z=unique([par.res.depth_km(:);cumsum(TRUEmodel.model_matrix(1:end-1,1))]);
zmax=max(z); assert(zmax>par.mod.maxz,'thbi:plot','Depth grid must extend into the half-space.');
if nc>size(chain_clrs,1), chain_clrs=[chain_clrs;lines(nc-size(chain_clrs,1))]; end
if show_figures, visibility='on'; else, visibility='off'; end
savedir=fullfile(result_dir,'Figures_THBI_Chains');
if ~exist(savedir,'dir'), mkdir(savedir); end
fprintf('Reading %d completed chains; burn-in %d; %d retained states/chain.\n%s\n', ...
    nc,par.inv.burnin,n,result_dir);

%% =================== COMMON-DEPTH PROFILES AND EXACT ACCEPTED STEP CURVES ===================
% Common depths are ONLY for mean/quantiles. Ensemble lines use exact physical
% step coordinates, so a very thin layer is not lost between grid samples.
% An accepted noise move is still accepted, but does not change the Earth.
profiles=zeros(n,numel(z),nc);
accepted_iterations=cell(nc,1); shown_iterations=cell(nc,1);
ensemble_vs=cell(nc,1); ensemble_z=cell(nc,1);
best_models=cell(nc,1); last_models=cell(nc,1);
best_vs=cell(nc,1); best_z=cell(nc,1); last_vs=cell(nc,1); last_z=cell(nc,1);
start_vs=cell(nc,1); start_z=cell(nc,1);
chain_mean=zeros(numel(z),nc); chain_quantiles=zeros(5,numel(z),nc);
chain_info=struct([]);
interface_edges=par.res.interface_bins_km(:)';
assert(interface_edges(1)==0 && interface_edges(end)>=par.mod.maxz, ...
    'thbi:plot','Interface bins must cover every inferred physical interface.');
interface_counts=zeros(nc,numel(interface_edges)-1);
for c=1:nc
    chain=chains{c};
    assert(chain.last_iter==par.inv.niter && size(chain.K,1)==par.inv.niter && ...
        numel(chain.accepted)==par.inv.niter && numel(chain.move)==par.inv.niter, ...
        'thbi:plot','Chain %d is incomplete or has incompatible storage.',c);
    assert(all(chain.accepted==0 | chain.accepted==1),'thbi:plot','Invalid acceptance flags.');
    accepted=logical(chain.accepted(keep));
    chain_info(c).accepted_proposals=sum(accepted);
    chain_info(c).accepted_structural=sum(accepted & chain.move(keep)~=6);
    chain_info(c).accepted_noise=sum(accepted & chain.move(keep)==6);
    if ~include_noise_accepts, accepted=accepted & chain.move(keep)~=6; end
    accepted_iterations{c}=keep(accepted);
    selected=accepted_iterations{c};
    if numel(selected)>max_models_per_chain
        selected=selected(unique(round(linspace(1,numel(selected),max_models_per_chain))));
    end
    shown_iterations{c}=selected;
    draw_row=zeros(par.inv.niter,1); draw_row(selected)=1:numel(selected);
    xv=cell(numel(selected),1); zv=cell(numel(selected),1);
    E2=chain.E2(keep,:);
    assert(isreal(E2) && all(isfinite(E2(:))) && all(E2(:)>=0), ...
        'thbi:plot','Invalid retained residual energies in chain %d.',c);
    [~,ibest]=min(sum(E2,2));
    chain_info(c).best_iteration=keep(ibest);
    chain_info(c).best_base_energy=sum(E2(ibest,:));
    chain_info(c).last_iteration=par.inv.niter;
    chain_info(c).retained=n; chain_info(c).shown=numel(selected);
    chain_info(c).forward_failures=chain.failures;
    for j=1:n
        i=keep(j); model=model_from_chain(chain,i,par);
        [ok,why]=a1_TEST_CONDITIONS(model,par);
        assert(ok && isequal(model.K,chain.K(i,:)), ...
            'thbi:plot','Invalid model in chain %d iteration %d: %s',c,i,why);
        material=sample_layer_model(model.matrix,z);
        profiles(j,:,c)=material(:,2)';
        interfaces=cumsum(model.matrix(1:end-1,1));
        % Exclude the fixed model-base interface from inferred-interface counts.
        interface_counts(c,:)=interface_counts(c,:)+histcounts(interfaces(1:end-1),interface_edges);
        if draw_row(i)>0
            edges=[0;interfaces;zmax];
            xv{draw_row(i)}=[reshape([model.matrix(:,3) model.matrix(:,3)]',[],1);NaN];
            zv{draw_row(i)}=[reshape([edges(1:end-1) edges(2:end)]',[],1);NaN];
        end
    end
    % NaN separates models in one line object per chain; no inter-model joins.
    if isempty(selected)
        ensemble_vs{c}=NaN; ensemble_z{c}=NaN;
    else
        ensemble_vs{c}=vertcat(xv{:}); ensemble_z{c}=vertcat(zv{:});
    end
    best_models{c}=model_from_chain(chain,keep(ibest),par);
    last_models{c}=model_from_chain(chain,par.inv.niter,par);
    chain_info(c).best_K=best_models{c}.K; chain_info(c).last_K=last_models{c}.K;
    model_list={best_models{c}.matrix,last_models{c}.matrix,chain.initial.model.matrix};
    step_vs=cell(3,1); step_z=cell(3,1);
    for j=1:3
        matrix=model_list{j}; edges=[0;cumsum(matrix(1:end-1,1));zmax];
        step_vs{j}=reshape([matrix(:,3) matrix(:,3)]',[],1);
        step_z{j}=reshape([edges(1:end-1) edges(2:end)]',[],1);
    end
    best_vs{c}=step_vs{1}; best_z{c}=step_z{1};
    last_vs{c}=step_vs{2}; last_z{c}=step_z{2};
    start_vs{c}=step_vs{3}; start_z{c}=step_z{3};
    chain_mean(:,c)=mean(profiles(:,:,c),1)';
    chain_quantiles(:,:,c)=sample_quantiles(profiles(:,:,c),[.025 .16 .5 .84 .975]);
    fprintf('Ch %d: %d accepted (%d structural, %d noise); showing %d; best retained iteration %d, K [%d %d %d].\n', ...
        c,chain_info(c).accepted_proposals,chain_info(c).accepted_structural,chain_info(c).accepted_noise, ...
        numel(selected),keep(ibest),best_models{c}.K);
    if isempty(selected)
        warning('thbi:plot','Chain %d has no accepted models matching the display selection; summaries still use its retained states.',c);
    end
end

%% =================== RESIDENCE-WEIGHTED SUMMARIES, TRUTH AND INTERFACE OCCURRENCE ===================
% The chain already samples the target: no exp(logL) weighting, no likelihood
% percentile cut, no accepted-only uncertainty calculation and no slot means.
pooled=reshape(permute(profiles,[1 3 2]),n*nc,numel(z));
pooled_mean=mean(pooled,1);
pooled_quantiles=sample_quantiles(pooled,[.025 .16 .5 .84 .975]);
interface_rate=bsxfun(@rdivide,interface_counts/n,diff(interface_edges));
interface_z=reshape([interface_edges(1:end-1);interface_edges(2:end)],[],1);
true_matrix=TRUEmodel.model_matrix;
true_interfaces=cumsum(true_matrix(1:end-1,1));
edges=[0;true_interfaces;zmax];
true_vs=reshape([true_matrix(:,3) true_matrix(:,3)]',[],1);
true_z=reshape([edges(1:end-1) edges(2:end)]',[],1);
[~,best_chain]=min([chain_info.best_base_energy]);
count_values=cell(1,3); count_probability=cell(1,3);
for zone=1:3
    count_values{zone}=par.mod.kmin(zone):par.mod.kmax(zone);
    count_probability{zone}=zeros(numel(count_values{zone}),nc);
    for c=1:nc
        for j=1:numel(count_values{zone})
            count_probability{zone}(j,c)=mean(chains{c}.K(keep,zone)==count_values{zone}(j));
        end
    end
end
vs_limit=1.02*max([par.mod.vsmax par.mod.halfspace_bounds(2) true_matrix(:,3)']);
% Extend the 1-km view below the target interface: a cutoff at exactly 1 km
% hides a recovered discontinuity that lies slightly deeper than the truth.
view_limits=[0 zmax;0 .2;0 .5;0 1.2;0 5;20 zmax];
view_names={'full','200m','500m','1km','5km','deep'};
view_titles={'Full depth','Upper 200 m','Upper 500 m','Around 1 km (0-1.2 km)','Upper 5 km','Crust and mantle'};
nfig=0; figures=gobjects(size(view_limits,1)+3+double(plot_saved_responses),1);
figure_names=cell(size(figures));
selection_note='Accepted post-burn-in proposals (including noise-only repeats)';
if ~include_noise_accepts, selection_note='Accepted post-burn-in structural proposals'; end

%% =================== CHAIN-COLORED ACCEPTED MODELS WITH INTERFACE PANELS ===================
% Layout follows the supplied scripts: broad Vs panel + narrow depth-aligned
% interface panel. Every thin curve is a physical step model, not a spline.
for v=1:size(view_limits,1)
    limits=view_limits(v,:);
    f=figure('Visible',visibility,'Color','w','Position',[70 80 1150 850]);
    nfig=nfig+1; figures(nfig)=f; figure_names{nfig}=['accepted_vs_' view_names{v}];
    axV=axes('Parent',f,'Position',[.08 .12 .58 .78]); hold(axV,'on'); box(axV,'on');
    for c=1:nc
        pale=1-ensemble_strength+ensemble_strength*chain_clrs(c,:);
        plot(axV,ensemble_vs{c},ensemble_z{c},'-','Color',pale,'LineWidth',.45, ...
            'DisplayName',sprintf('Ch %d accepted (%d shown)',c,numel(shown_iterations{c})));
    end
    for c=1:nc
        plot(axV,best_vs{c},best_z{c},'-','Color',chain_clrs(c,:),'LineWidth',2.2, ...
            'DisplayName',sprintf('Ch %d best retained: K [%d %d %d]',c,best_models{c}.K));
    end
    plot(axV,true_vs,true_z,'k-','LineWidth',2.4,'DisplayName','True Earth');
    xupper=vs_limit;
    if limits(2)<=5
        observed=pooled(:,z<=limits(2)); maximum_vs=max(observed(:));
        for c=1:nc
            visible_vs=ensemble_vs{c}(ensemble_z{c}<=limits(2));
            if ~isempty(visible_vs), maximum_vs=max(maximum_vs,max(visible_vs)); end
            maximum_vs=max(maximum_vs,max(best_vs{c}(best_z{c}<=limits(2))));
        end
        maximum_vs=max(maximum_vs,max(true_vs(true_z<=limits(2))));
        xupper=min(vs_limit,1.1*maximum_vs);
    end
    xlim(axV,[0 xupper]); ylim(axV,limits);
    set(axV,'YDir','reverse','FontSize',13,'LineWidth',1.3,'Layer','top'); grid(axV,'on');
    xlabel(axV,'Vs (km/s)'); ylabel(axV,'Depth (km)'); title(axV,view_titles{v});
    legend(axV,'Location','best','FontSize',9);
    axD=axes('Parent',f,'Position',[.73 .12 .23 .78]); hold(axD,'on'); box(axD,'on');
    for c=1:nc
        rate=reshape([interface_rate(c,:);interface_rate(c,:)],[],1);
        plot(axD,rate,interface_z,'-','Color',chain_clrs(c,:),'LineWidth',1.3, ...
            'DisplayName',sprintf('Chain %d',c));
    end
    for depth=true_interfaces(1:end-1)'
        if depth>=limits(1) && depth<=limits(2)
            yline(axD,depth,'k:','HandleVisibility','off');
        end
    end
    bins=interface_edges(1:end-1)<limits(2) & interface_edges(2:end)>limits(1);
    selected_rate=interface_rate(:,bins); xlim(axD,[0 max(1,1.05*max(selected_rate(:)))]);
    set(axD,'YDir','reverse','FontSize',12,'LineWidth',1.3,'YTickLabel',[],'Layer','top');
    ylim(axD,limits); grid(axD,'on'); xlabel(axD,{'Expected interfaces/km','All retained states'});
    title(axD,'Physical interfaces'); linkaxes([axV axD],'y');
    sgtitle(f,selection_note,'FontSize',13);
end

%% =================== PER-CHAIN FINAL EARTHS AND ALL-RETAINED UNCERTAINTY ===================
% "Best retained" minimizes sum(E2) within that chain after burn-in. It is
% an actual sampled model, not a cross-dimension MAP or an averaged layering.
f=figure('Visible',visibility,'Color','w','Position',[80 60 1450 max(650,360*nc)]);
nfig=nfig+1; figures(nfig)=f; figure_names{nfig}='final_vs_by_chain';
tiledlayout(f,nc,3,'TileSpacing','compact','Padding','compact');
final_limits=[zmax .2 5];
for c=1:nc
    for v=1:3
        ax=nexttile; hold(ax,'on'); box(ax,'on'); q=chain_quantiles(:,:,c);
        fill(ax,[q(1,:) fliplr(q(5,:))],[z' fliplr(z')], ...
            1-.15+.15*chain_clrs(c,:),'EdgeColor','none','DisplayName','Retained 95% interval');
        plot(ax,ensemble_vs{c},ensemble_z{c},'Color',1-.12+.12*chain_clrs(c,:), ...
            'LineWidth',.35,'HandleVisibility','off');
        plot(ax,start_vs{c},start_z{c},'--','Color',[.55 .55 .55], ...
            'LineWidth',1.2,'DisplayName','Starting Earth');
        plot(ax,q(3,:),z,':','Color',chain_clrs(c,:),'LineWidth',2,'DisplayName','Retained median Vs(z)');
        plot(ax,last_vs{c},last_z{c},'-.','Color',chain_clrs(c,:), ...
            'LineWidth',1.7,'DisplayName',sprintf('Last state: K [%d %d %d]',last_models{c}.K));
        plot(ax,best_vs{c},best_z{c},'-','Color',chain_clrs(c,:), ...
            'LineWidth',2.3,'DisplayName',sprintf('Best retained: K [%d %d %d]',best_models{c}.K));
        plot(ax,true_vs,true_z,'k-','LineWidth',2.3,'DisplayName','True Earth');
        set(ax,'YDir','reverse','FontSize',11,'LineWidth',1.2,'Layer','top'); grid(ax,'on');
        ylim(ax,[0 final_limits(v)]); xlim(ax,[0 vs_limit]);
        xlabel(ax,'Vs (km/s)'); ylabel(ax,'Depth (km)');
        title(ax,sprintf('Chain %d | upper %.3g km | best iteration %d',c,final_limits(v),chain_info(c).best_iteration));
        if v==1, legend(ax,'Location','best','FontSize',8); end
    end
end
sgtitle(f,'Final physical Earths and depth summaries by chain; inspect chain agreement','FontSize',14);

%% =================== POOLED MEAN, MEDIAN, INTERVALS AND BEST RETAINED MODEL ===================
f=figure('Visible',visibility,'Color','w','Position',[90 60 1100 900]);
nfig=nfig+1; figures(nfig)=f; figure_names{nfig}='retained_vs_summary';
tiledlayout(f,2,2,'TileSpacing','compact','Padding','compact');
for limit=[zmax .2 .5 5]
    ax=nexttile; hold(ax,'on'); box(ax,'on'); q=pooled_quantiles;
    fill(ax,[q(1,:) fliplr(q(5,:))],[z' fliplr(z')],[.92 .92 .92], ...
        'EdgeColor','none','DisplayName','Retained 95% interval');
    fill(ax,[q(2,:) fliplr(q(4,:))],[z' fliplr(z')],[.78 .78 .78], ...
        'EdgeColor','none','DisplayName','Retained 68% interval');
    plot(ax,pooled_mean,z,'-','Color',[.1 .3 .9],'LineWidth',1.8,'DisplayName','Pointwise mean');
    plot(ax,q(3,:),z,'--','Color',[.55 0 .55],'LineWidth',1.8,'DisplayName','Pointwise median');
    plot(ax,best_vs{best_chain},best_z{best_chain},'-','Color',[.9 .3 0], ...
        'LineWidth',1.8,'DisplayName',sprintf('Best retained Earth (Ch %d)',best_chain));
    plot(ax,true_vs,true_z,'k-','LineWidth',2.3,'DisplayName','True Earth');
    set(ax,'YDir','reverse','FontSize',12,'LineWidth',1.2,'Layer','top'); grid(ax,'on');
    xlim(ax,[0 vs_limit]); ylim(ax,[0 limit]); xlabel(ax,'Vs (km/s)'); ylabel(ax,'Depth (km)');
    title(ax,sprintf('Upper %.3g km',limit));
    if limit==zmax, legend(ax,'Location','best','FontSize',8); end
end
sgtitle(f,sprintf('All %d retained states: repeats preserved; pooled summaries require chain agreement',n*nc),'FontSize',13);

%% =================== LAYER-COUNT MASSES AND CHAINS, WITH THE SAME COLORS ===================
f=figure('Visible',visibility,'Color','w','Position',[90 80 1250 720]);
nfig=nfig+1; figures(nfig)=f; figure_names{nfig}='layer_counts_by_chain';
tiledlayout(f,ceil((numel(par.inv.datatypes)+1)/3),3,'TileSpacing','compact','Padding','compact');
for zone=1:3
    ax=nexttile; hold(ax,'on'); box(ax,'on');
    bars=bar(ax,count_values{zone},count_probability{zone},'grouped');
    for c=1:nc, bars(c).FaceColor=chain_clrs(c,:); bars(c).DisplayName=sprintf('Chain %d',c); end
    xline(ax,TRUEmodel.zone_counts(zone),'k--','True K','HandleVisibility','off');
    xticks(ax,count_values{zone}); ylim(ax,[0 1]); grid(ax,'on');
    xlabel(ax,'Finite physical layers'); ylabel(ax,'Retained probability'); title(ax,par.mod.zone_names{zone});
    if zone==1, legend(ax,'Location','best'); end
end
for zone=1:3
    ax=nexttile; hold(ax,'on'); box(ax,'on');
    for c=1:nc
        stairs(ax,1:par.inv.niter,chains{c}.K(:,zone),'Color',chain_clrs(c,:),'LineWidth',.9);
    end
    xline(ax,par.inv.burnin,'k:','Burn-in'); yline(ax,TRUEmodel.zone_counts(zone),'k--');
    ylim(ax,[par.mod.kmin(zone)-.5 par.mod.kmax(zone)+.5]); grid(ax,'on');
    xlabel(ax,'Iteration'); ylabel(ax,'K'); title(ax,par.mod.zone_names{zone});
end
sgtitle(f,'Variable physical layer counts: half-space and numerical subdivisions excluded','FontSize',14);

%% =================== OPTIONAL FIVE OBSERVABLES FROM SAVED RESPONSE DRAWS ===================
% Responses were cached every prediction_stride iterations, not at every
% accepted proposal. Match by prediction_iterations, not by ensemble row.
response_rows=cell(nc,1); response_iterations=cell(nc,1);
if plot_saved_responses
    f=figure('Visible',visibility,'Color','w','Position',[100 80 1400 850]);
    nfig=nfig+1; figures(nfig)=f; figure_names{nfig}='saved_five_observables_by_chain';
    tiledlayout(f,ceil((numel(par.inv.datatypes)+1)/3),3,'TileSpacing','compact','Padding','compact');
    labels={'HVSR','Rayleigh ellipticity','Phase velocity','Group velocity','Ps receiver function','Sp receiver function'};
    units={'H/V','|H/V|','km/s','km/s','R/Z amplitude','Z/R amplitude'};
    for d=1:numel(par.inv.datatypes)
        ax=nexttile; hold(ax,'on'); box(ax,'on'); name=par.inv.datatypes{d}; obs=data.(name);
        for c=1:nc
            chain=chains{c};
            assert(numel(chain.prediction_iterations)==numel(chain.prediction_saved), ...
                'thbi:plot','Invalid prediction indexing for chain %d.',c);
            rows=find(chain.prediction_saved & chain.prediction_iterations>par.inv.burnin & ...
                chain.prediction_iterations<=par.inv.niter);
            if isempty(rows), continue; end
            assert(isfield(chain.predictions,name),'thbi:plot','Missing cached %s predictions.',name);
            values=chain.predictions.(name)(rows,:);
            assert(size(values,2)==obs.N && isreal(values) && all(isfinite(values(:))), ...
                'thbi:plot','Invalid cached %s predictions in chain %d.',name,c);
            pick=unique(round(linspace(1,numel(rows),min(max_responses_per_chain,numel(rows)))));
            response_rows{c}=rows(pick); response_iterations{c}=chain.prediction_iterations(rows(pick));
            h=plot(ax,obs.x,values(pick,:)','Color',1-.15+.15*chain_clrs(c,:),'LineWidth',.4);
            set(h,'HandleVisibility','off');
            q=sample_quantiles(values,.5);
            plot(ax,obs.x,q,'Color',chain_clrs(c,:),'LineWidth',1.7, ...
                'DisplayName',sprintf('Ch %d saved-response median',c));
        end
        plot(ax,obs.x,obs.obs,'.','Color',[.45 .45 .45],'MarkerSize',6,'DisplayName','Noisy synthetic');
        plot(ax,obs.x,obs.clean,'k-','LineWidth',1.8,'DisplayName','True response');
        if d==1
            set(ax,'XScale','log'); xlabel(ax,'Frequency (Hz)');
        elseif d==5
            xlabel(ax,'Lag relative to P (s)');
        elseif d==6
            xlabel(ax,'Lag relative to S (s)');
        else
            xlabel(ax,'Period (s)');
        end
        xlim(ax,[obs.x(1) obs.x(end)]); ylabel(ax,units{d}); title(ax,labels{d}); grid(ax,'on');
        if d==1, legend(ax,'Location','best','FontSize',8); end
    end
    ax=nexttile; axis(ax,'off');
    text(ax,0,1,{'Saved retained response draws only'; ...
        'Not every accepted Vs model has a cached response.'; ...
        'Thick colored curves: per-chain response medians'; ...
        'Black: known synthetic response'; ...
        'No new forward calculations.'},'Units','normalized','VerticalAlignment','top','FontSize',11);
    sgtitle(f,'Joint observables by chain from saved model responses','FontSize',14);
end

%% =================== SAVE FIGURES AND EXACT MODEL/ITERATION SELECTIONS ===================
assert(nfig==numel(figures),'thbi:plot','Figure collection is incomplete.');
for j=1:numel(figures)
    exportgraphics(figures(j),fullfile(savedir,[figure_names{j} '.png']),'Resolution',200);
    exportgraphics(figures(j),fullfile(savedir,[figure_names{j} '.pdf']),'ContentType','image','Resolution',300);
    if save_matlab_figures, savefig(figures(j),fullfile(savedir,[figure_names{j} '.fig'])); end
    if ~show_figures, close(figures(j)); end
end
plot_summary=struct();
plot_summary.result_dir=result_dir; plot_summary.burnin=par.inv.burnin;
plot_summary.retained_iterations=keep; plot_summary.accepted_iterations=accepted_iterations;
plot_summary.shown_iterations=shown_iterations; plot_summary.chain_info=chain_info;
plot_summary.selection_note=selection_note; plot_summary.chain_colors=chain_clrs(1:nc,:);
plot_summary.depth_km=z; plot_summary.quantile_probabilities=[.025 .16 .5 .84 .975];
plot_summary.chain_mean=chain_mean; plot_summary.chain_quantiles=chain_quantiles;
plot_summary.pooled_mean=pooled_mean; plot_summary.pooled_quantiles=pooled_quantiles;
plot_summary.best_models=best_models; plot_summary.last_models=last_models;
plot_summary.interface_edges=interface_edges; plot_summary.interface_rate_per_km=interface_rate;
plot_summary.count_values=count_values; plot_summary.count_probability=count_probability;
plot_summary.response_rows=response_rows; plot_summary.response_iterations=response_iterations;
plot_summary.note='Summaries use every retained state with equal weight, including repeats; accepted-only lines are a display selection.';
plot_summary.representative_selection='Minimum sum(E2) among retained states in each chain; not MAP.';
plot_summary.figure_names=figure_names;
save(fullfile(savedir,'plot_summary.mat'),'plot_summary','-v7');
fid=fopen(fullfile(savedir,'plot_notes.txt'),'w'); assert(fid>=0);
fprintf(fid,'Source run: %s\nBurn-in %d; %d retained states per chain.\n',result_dir,par.inv.burnin,n);
fprintf(fid,'%s\n%s\n%s\n',selection_note,plot_summary.note,plot_summary.representative_selection);
for c=1:nc
    fprintf(fid,'Ch %d: accepted %d (structural %d, noise %d); shown %d; best iteration %d K [%d %d %d]; last K [%d %d %d]; native failures %d.\n', ...
        c,chain_info(c).accepted_proposals,chain_info(c).accepted_structural,chain_info(c).accepted_noise, ...
        chain_info(c).shown,chain_info(c).best_iteration,chain_info(c).best_K,chain_info(c).last_K,chain_info(c).forward_failures);
end
fprintf(fid,'Interface panel: expected number of inferred physical interfaces per km; excludes fixed half-space top.\n');
fprintf(fid,'Intervals: 2.5-97.5 and 16-84 percentiles; pointwise means/medians need not be a sampled Earth.\n');
fprintf(fid,'Saved response overlays are indexed by their own iterations, independent of accepted-model display selection.\n');
if isfield(saved,'summary') && isfield(saved.summary,'split_rhat')
    fprintf(fid,'Saved classical split R-hat (inspect chain agreement):\n');
    for j=1:numel(saved.summary.metric_names)
        fprintf(fid,'  %s: %.6g\n',saved.summary.metric_names{j},saved.summary.split_rhat(j));
    end
end
fprintf(fid,'Inspect mixing and native failures before interpreting pooled intervals as posterior uncertainty.\n');
fclose(fid);
fprintf('\nSaved %d figure sets and plot_summary.mat:\n%s\n',numel(figures),savedir);
