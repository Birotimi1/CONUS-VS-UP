function c2_PLOT_RESULTS(par, data, TRUEmodel, chains, summary)
% MATLAB-native, reproducible scientific figures. Fits and mixing separate.

%% Figure settings and chosen states
out=par.res.dir; best=summary.bestfit; initial=chains{1}.initial;
blue=[.12 .38 .75]; orange=[.87 .29 .08]; gray=[.5 .5 .5];
visibility='off';
if par.res.show_figures, visibility='on'; end
%% Figure 1: layered Vs structure and all five observable fits
f=figure('Visible',visibility,'Color','w','Position',[50 50 1400 860]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
z=summary.depth_km; q=summary.vs_quantiles;
fill([q(1,:) fliplr(q(3,:))],[z' fliplr(z')],[.88 .92 .98], ...
    'EdgeColor','none','DisplayName','Retained 5-95% (inspect mixing)');
plot_layer_model(TRUEmodel.model_matrix,'k',1.8,'True');
plot_layer_model(initial.model.matrix,gray,1.1,'Start, chain 1');
plot_layer_model(best.model.matrix,orange,1.5,'Best visited fit');
set(gca,'YDir','reverse'); ylim([0 85]); xlim([0 5.1]); grid on;
xlabel('Vs (km/s)'); ylabel('Depth (km)'); title('Fixed-layer Earth');
legend('Location','southwest','FontSize',8);
labels={'HVSR','Rayleigh ellipticity','Rayleigh phase velocity','Rayleigh group velocity','P receiver function'};
units={'H/V','|H/V|','km/s','km/s','R/Z amplitude'};
for k=1:5
    d=par.inv.datatypes{k}; v=data.(d);
    nexttile; hold on;
    plot(v.x,v.obs,'.','Color',gray,'MarkerSize',6,'DisplayName','Noisy synthetic');
    plot(v.x,v.clean,'k-','LineWidth',1.5,'DisplayName','True response');
    plot(v.x,initial.predata.(d),'--','Color',blue,'LineWidth',1.0,'DisplayName','Start, chain 1');
    plot(v.x,best.predata.(d),'-','Color',orange,'LineWidth',1.4,'DisplayName','Best visited fit');
    if k==1, set(gca,'XScale','log'); xlabel('Frequency (Hz)');
    elseif k==5, xlabel('Time after P (s)');
    else, xlabel('Period (s)'); end
    ylabel(units{k}); title(sprintf('%s | base chi2/N = %.2f',labels{k},summary.base_reduced_chi2(k)));
    grid on; xlim([v.x(1) v.x(end)]);
    if k==1, legend('Location','best','FontSize',8); end
end
sgtitle(sprintf('Five-observable thin-model HBI | %s: %d chains x %d steps | convergence not established', ...
    par.inv.profile,par.inv.nchains,par.inv.niter),'FontWeight','bold');
exportgraphics(f,fullfile(out,'joint_fit.png'),'Resolution',170);
exportgraphics(f,fullfile(out,'joint_fit.pdf'),'ContentType','vector');

%% Figure 2: chain histories, noise multipliers, and convergence diagnostics
g=figure('Visible',visibility,'Color','w','Position',[50 50 1400 850]);
if ~par.res.show_figures, cleanup2=onCleanup(@() close(g)); end %#ok<NASGU>
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for c=1:numel(chains), plot(chains{c}.logL); end
xline(par.inv.burnin,'k:'); xlabel('Iteration'); ylabel('log likelihood'); title('Chain histories'); grid on;
nexttile; hold on;
for c=1:numel(chains), semilogy(sum(chains{c}.E2,2)); end
xline(par.inv.burnin,'k:'); xlabel('Iteration'); ylabel('Sum of base chi2'); title('Fit independent of inferred noise'); grid on;
nexttile; hold on;
for c=1:numel(chains), plot(sum(chains{c}.theta(:,1:7),2)); end
yline(TRUEmodel.moho_depth,'k--'); xline(par.inv.burnin,'k:');
xlabel('Iteration'); ylabel('Depth of interface 7 (km)'); title('Moho-labelled interface'); grid on;
nexttile; hold on;
for c=1:numel(chains), plot(sum(chains{c}.theta(:,1:2),2)*1000); end
yline(TRUEmodel.sed_thickness*1000,'k--'); xline(par.inv.burnin,'k:');
xlabel('Iteration'); ylabel('First two layers (m)'); title('Soft-cap thickness'); grid on;
nexttile; hold on;
colors=lines(5);
noise_handles=gobjects(5,1);
for k=1:5
    for c=1:numel(chains)
        style='-'; if c>1, style='--'; end
        p=plot(exp(chains{c}.eta(:,k)),'Color',colors(k,:),'LineStyle',style);
        if c==1, noise_handles(k)=p; else, p.HandleVisibility='off'; end
    end
end
yline(1,'k:'); xline(par.inv.burnin,'k:'); xlabel('Iteration'); ylabel('Noise multiplier gamma');
title('Noise multipliers (line styles distinguish chains)'); grid on;
legend(noise_handles,par.inv.datatypes,'Location','best','FontSize',8);
nexttile; bar(summary.split_rhat); yline(1.01,'r--');
xticks(1:24); xticklabels(summary.parameter_names); xtickangle(70);
ylabel('Classical split R-hat'); title('Mixing diagnostic; not a convergence certificate'); grid on;
sgtitle('HBI diagnostics: inspect before interpreting retained-sample uncertainty');
exportgraphics(g,fullfile(out,'diagnostics.png'),'Resolution',170);
exportgraphics(g,fullfile(out,'diagnostics.pdf'),'ContentType','vector');

%% Figure 3: shallow section showing the thin near-surface layers
h=figure('Visible',visibility,'Color','w','Position',[50 50 600 650]);
if ~par.res.show_figures, cleanup3=onCleanup(@() close(h)); end %#ok<NASGU>
hold on; plot_layer_model(TRUEmodel.model_matrix,'k',2,'True');
plot_layer_model(initial.model.matrix,gray,1.2,'Start, chain 1');
plot_layer_model(best.model.matrix,orange,1.8,'Best visited fit');
set(gca,'YDir','reverse'); ylim([0 .3]); xlim([0 3.3]); grid on;
xlabel('Vs (km/s)'); ylabel('Depth (km)'); title('Shallow thin-model structure'); legend('Location','southwest');
exportgraphics(h,fullfile(out,'shallow_model.png'),'Resolution',170);
write_hbi_report(out,par,chains,summary);
end
