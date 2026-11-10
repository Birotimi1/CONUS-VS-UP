function plot_reference_synthetics(truth,data,par)
% Show actual synthetic inputs before inversion; no posterior is implied.
visible='off'; if par.res.show_figures, visible='on'; end
f=figure('Visible',visible,'Color','w','Position',[50 50 1450 900]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for limit=[.1 250]
    nexttile; plot_layer_model(truth.model_matrix,'k',1.8,'ECM1 + AK135, THBI materials');
    set(gca,'YDir','reverse'); ylim([0 limit]); grid on;
    xlabel('Vs (km/s)'); ylabel('Depth (km)'); title(sprintf('Known structure: upper %g km',limit));
end
labels={'HVSR','Rayleigh ellipticity','Rayleigh phase velocity','Rayleigh group velocity'};
units={'H/V','|H/V|','Phase velocity (km/s)','Group velocity (km/s)'};
for k=1:4
    name=par.inv.datatypes{k}; v=data.(name); nexttile;
    plot(v.x,v.clean,'k-','LineWidth',1.3); hold on; plot(v.x,v.obs,'.','Color',[.55 .55 .55]);
    if k==1, set(gca,'XScale','log'); xlabel('Frequency (Hz)'); else, xlabel('Period (s)'); end
    grid on; title(labels{k}); ylabel(units{k}); xlim([v.x(1) v.x(end)]);
end
nexttile; plot(data.rf.waves.t,data.rf.waves.rf); xlim([-5 40]); grid on;
xlabel('Ps conversion lag (s)'); ylabel('R/Z'); title('Observed synthetic Ps waves for H-kappa');
nexttile; hk=data.rf.true_stack;
imagesc(hk.H_km,hk.kappa,hk.energy/hk.maximum_energy); axis xy; colorbar; hold on;
plot(truth.moho_depth,truth.zoned_parameters.vpvs,'wx','LineWidth',2,'MarkerSize',10);
xlabel('Moho depth H (km)'); ylabel('Vp/Vs'); title('Ps stack at true mean crustal Vs');
nexttile; plot(data.rf_sp.x,data.rf_sp.clean,'k-','LineWidth',1.3); hold on;
plot(data.rf_sp.x,data.rf_sp.obs,'.','Color',[.55 .55 .55]); grid on;
xlabel('Conversion depth (km)'); ylabel('P daughter amplitude'); title('Sp CCP depth stack');
sgtitle(truth.description,'Interpreter','none');
exportgraphics(f,fullfile(par.res.dir,'reference_synthetic_inputs.png'),'Resolution',180);
exportgraphics(f,fullfile(par.res.dir,'reference_synthetic_inputs.pdf'),'ContentType','vector');
if ~par.res.show_figures, close(f); end
end
