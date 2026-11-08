function c2_PLOT_RESULTS(par, data, TRUEmodel, chains, summary, posterior, suite_of_models, final_model)
% MATLAB plot blocks for transdimensional structure, fits, counts and mixing.
% Helper plot/report bodies are in functions/, following THBI organization.
nd=numel(par.inv.datatypes);
out=par.res.dir; best=summary.bestfit; initial=chains{1}.initial;
fit_label='Best visited fit (legacy)';
if isfield(summary,'bestfit_label'), fit_label=summary.bestfit_label; end
if nargin>=8 && isfield(summary,'bestfit_iteration') && isfinite(summary.bestfit_iteration)
    assert(summary.bestfit_chain==final_model.representative_chain && ...
        summary.bestfit_iteration==final_model.representative_iteration, ...
        'thbi:results','Fit plot and final-forward selection disagree.');
end
visibility='off'; if par.res.show_figures, visibility='on'; end
orange=[.87 .29 .08]; blue=[.12 .38 .75]; gray=[.55 .55 .55];
figures=gobjects(0);

%% =================== SHARED EARTH AND FIVE OBSERVABLE FITS ===================
f=figure('Visible',visibility,'Color','w','Position',[50 50 1400 860]); figures(end+1)=f;
tiledlayout(ceil((nd+1)/3),3,'TileSpacing','compact','Padding','compact');
nexttile; hold on; z=summary.depth_km; q=summary.vs_quantiles;
if ~isempty(summary.prior)
    pq=summary.prior.vs_quantiles;
    fill([pq(1,:) fliplr(pq(3,:))],[z' fliplr(z')],[.94 .94 .94], ...
        'EdgeColor','none','DisplayName','Conditioned prior 5-95%');
end
fill([q(1,:) fliplr(q(3,:))],[z' fliplr(z')],[.88 .92 .98], ...
    'EdgeColor','none','DisplayName','Retained 5-95% (inspect mixing)');
plot_layer_model(TRUEmodel.model_matrix,'k',1.8,'True');
plot_layer_model(initial.model.matrix,gray,1.1,'Start, chain 1');
plot_layer_model(best.model.matrix,orange,1.5,fit_label);
set(gca,'YDir','reverse'); ylim([0 max(summary.depth_km)]); xlim([0 5.2]); grid on;
xlabel('Vs (km/s)'); ylabel('Depth (km)'); title('Variable-layer Earth');
legend('Location','southwest','FontSize',8);
labels={'HVSR','Rayleigh ellipticity','Rayleigh phase velocity','Rayleigh group velocity','Ps receiver function','Sp receiver function'};
units={'H/V','|H/V|','km/s','km/s','R/Z amplitude','Z/R amplitude'};
for d=1:nd
    name=par.inv.datatypes{d}; v=data.(name);
    nexttile; hold on;
    if isfield(summary.response_quantiles,name)
        q=summary.response_quantiles.(name);
        fill([v.x' fliplr(v.x')],[q(1,:) fliplr(q(3,:))],[.88 .92 .98], ...
            'EdgeColor','none','DisplayName','Model responses 5-95%');
    end
    plot(v.x,v.obs,'.','Color',gray,'MarkerSize',6,'DisplayName','Noisy synthetic');
    plot(v.x,v.clean,'k-','LineWidth',1.5,'DisplayName','True response');
    plot(v.x,initial.predata.(name),'--','Color',blue,'LineWidth',1,'DisplayName','Start, chain 1');
    plot(v.x,best.predata.(name),'-','Color',orange,'LineWidth',1.4,'DisplayName',fit_label);
    if d==1
        set(gca,'XScale','log'); xlabel('Frequency (Hz)');
    elseif d==5
        xlabel('Lag relative to P (s)');
    elseif d==6
        xlabel('Lag relative to S (s)');
    else
        xlabel('Period (s)');
    end
    ylabel(units{d}); grid on; xlim([v.x(1) v.x(end)]);
    title(sprintf('%s | base chi2/N %.2f',labels{d},summary.base_reduced_chi2(d)));
    if d==1, legend('Location','best','FontSize',8); end
end
sgtitle(sprintf('Zoned transdimensional HBI | %s | convergence requires inspection',par.inv.profile));
exportgraphics(f,fullfile(out,'joint_fit.png'),'Resolution',170);
exportgraphics(f,fullfile(out,'joint_fit.pdf'),'ContentType','vector');

%% =================== POSTERIOR PHYSICAL LAYER COUNTS ===================
f=figure('Visible',visibility,'Color','w','Position',[80 80 1200 400]); figures(end+1)=f;
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for z=1:3
    nexttile;
    if isempty(summary.prior)
        bar(summary.count_values{z},summary.count_probability{z},'FaceColor',blue);
    else
        bar(summary.count_values{z},[summary.prior.count_probability{z}(:),summary.count_probability{z}(:)]);
        legend({'Conditioned prior','Retained posterior'},'Location','best','FontSize',8);
    end
    xline(TRUEmodel.zone_counts(z),'k--','True count','HandleVisibility','off');
    xlabel('Finite physical layers'); ylabel('Probability');
    title(par.mod.zone_names{z}); grid on;
    xticks(summary.count_values{z}); ylim([0 1]);
end
sgtitle('Layer-count distributions: half-space and numerical subdivisions excluded');
exportgraphics(f,fullfile(out,'layer_counts.png'),'Resolution',170);
exportgraphics(f,fullfile(out,'layer_counts.pdf'),'ContentType','vector');

%% =================== COUNT, BOUNDARY, NOISE AND LIKELIHOOD HISTORIES ===================
f=figure('Visible',visibility,'Color','w','Position',[80 80 1400 1000]); figures(end+1)=f;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for z=1:3
    nexttile; hold on;
    for c=1:numel(chains), stairs(chains{c}.K(:,z)); end
    yline(TRUEmodel.zone_counts(z),'k--'); xline(par.inv.burnin,'k:');
    xlabel('Iteration'); ylabel('K'); title([par.mod.zone_names{z} ' layer count']); grid on;
end
nexttile; hold on;
for c=1:numel(chains), plot(chains{c}.boundaries(:,1)); end
yline(TRUEmodel.sediment_base_depth,'k--'); xline(par.inv.burnin,'k:');
xlabel('Iteration'); ylabel('Depth (km)'); title('Sediment base'); grid on;
nexttile; hold on;
for c=1:numel(chains), plot(chains{c}.boundaries(:,2)); end
yline(TRUEmodel.moho_depth,'k--'); xline(par.inv.burnin,'k:');
xlabel('Iteration'); ylabel('Depth (km)'); title('Moho'); grid on;
nexttile; hold on;
for c=1:numel(chains), plot(chains{c}.logL); end
xline(par.inv.burnin,'k:'); xlabel('Iteration'); ylabel('log L'); title('Joint likelihood'); grid on;
nexttile; hold on; colors=lines(nd); hh=gobjects(nd,1);
for d=1:nd
    for c=1:numel(chains)
        style='-'; if c>1, style='--'; end
        p=plot(exp(chains{c}.eta(:,d)),'Color',colors(d,:),'LineStyle',style);
        if c==1, hh(d)=p; else, p.HandleVisibility='off'; end
    end
end
yline(1,'k:'); xline(par.inv.burnin,'k:');
xlabel('Iteration'); ylabel('gamma'); title('Noise multipliers'); grid on;
legend(hh,par.inv.datatypes,'FontSize',8,'Location','best');
nexttile; hold on;
for c=1:numel(chains), semilogy(sum(chains{c}.E2,2)); end
xline(par.inv.burnin,'k:'); xlabel('Iteration'); ylabel('Sum base chi2'); title('Fit before noise scaling'); grid on;
nexttile; bar(summary.split_rhat); yline(1.01,'r--');
xticks(1:numel(summary.metric_names)); xticklabels(summary.metric_names); xtickangle(60);
ylabel('Classical split R-hat'); title('Fixed-metric mixing diagnostic'); grid on;
sgtitle('RJ chain diagnostics: counts must mix as well as velocities and noise');
exportgraphics(f,fullfile(out,'diagnostics.png'),'Resolution',170);

%% =================== SHALLOW STRUCTURE AT TWO DEPTH SCALES ===================
f=figure('Visible',visibility,'Color','w','Position',[100 100 900 600]); figures(end+1)=f;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for limit=[.2 2]
    nexttile; hold on;
    z=summary.depth_km; q=summary.vs_quantiles;
    fill([q(1,:) fliplr(q(3,:))],[z' fliplr(z')],[.88 .92 .98], ...
        'EdgeColor','none','DisplayName','Retained 5-95%');
    plot_layer_model(TRUEmodel.model_matrix,'k',2,'True');
    plot_layer_model(initial.model.matrix,gray,1.1,'Start, chain 1');
    plot_layer_model(best.model.matrix,orange,1.5,fit_label);
    set(gca,'YDir','reverse'); ylim([0 limit]); xlim([0 4]); grid on;
    xlabel('Vs (km/s)'); ylabel('Depth (km)'); title(sprintf('Upper %.0f m',1000*limit));
    legend('Location','southwest');
end
exportgraphics(f,fullfile(out,'shallow_model.png'),'Resolution',170);

%% =================== VELOCITY MIXING AT FIXED DEPTHS ===================
f=figure('Visible',visibility,'Color','w','Position',[100 100 900 600]); figures(end+1)=f;
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for limit=[2 max(summary.depth_km)]
    nexttile; hold on;
    r=summary.vs_split_rhat(:); z=summary.depth_km(:);
    finite_r=r; finite_r(~isfinite(finite_r))=NaN;
    plot(finite_r,z,'Color',blue,'LineWidth',1.3); xline(1.01,'r--');
    set(gca,'YDir','reverse'); ylim([0 limit]); grid on;
    xlabel('Classical split R-hat of Vs'); ylabel('Depth (km)');
    title(sprintf('Upper %.2g km',limit));
    selected=z<=limit;
    text(.03,.03,sprintf('Infinite: %d depths; unavailable: %d', ...
        sum(isinf(r(selected))),sum(isnan(r(selected)))), ...
        'Units','normalized','FontSize',8,'VerticalAlignment','bottom');
end
sgtitle('Velocity mixing: Inf = distinct constant segments; NaN = undefined or unavailable');
exportgraphics(f,fullfile(out,'velocity_mixing.png'),'Resolution',170);
write_thbi_report(par,chains,summary);
if ~par.res.show_figures, close(figures); end

%% =================== FULL THBI-STYLE RETAINED-MODEL DIAGNOSTICS ===================
% Five-argument calls retain the existing plot-only API. The DATA run and
% postprocessing scripts supply c2/c3/c4 outputs for the extended suite.
if nargin>=8
    plot_HEATMAP_ALLMODELS(suite_of_models,final_model,par);
    plot_HEATMAP_ALLMODELS(suite_of_models,final_model,par,.3);
    plot_SUITE_of_MODELS(suite_of_models,final_model,par);
    plot_TRADEOFF(posterior,par);
    plot_PARAMETERISATION(final_model.representative_model,par);
    plot_PERTURBATIONS(summary,par);
    plot_LAYER_TRENDS(chains,par);
    for c=1:numel(chains), plot_progress_chain(chains{c},par,c); end
    if ~isempty(summary.prior)
        plot_PRIORvsPOSTERIOR(summary.prior,posterior,final_model,par);
    end

    %% Signed RF wiggles: identical time grid and original amplitudes
    for d=5:nd
        name=par.inv.datatypes{d}; v=data.(name);
        f=figure('Visible',visibility,'Color','w','Position',[100 100 900 650]);
        tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
        traces={v.clean,v.obs,best.predata.(name)};
        titles={'True','Noisy synthetic',fit_label};
        for j=1:3
            ax=nexttile; plot_RF(ax,traces{j},v.x);
            title(ax,[titles{j} ' ' labels{d}]); ylim(ax,[v.x(1) v.x(end)]); grid(ax,'on');
            ylabel(ax,'Conversion lag (s)'); xlabel(ax,units{d});
        end
        exportgraphics(f,fullfile(out,[name '_wiggles.png']),'Resolution',170);
        if ~par.res.show_figures, close(f); end
    end
end
end
