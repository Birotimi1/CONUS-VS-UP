function plot_PRIORvsPOSTERIOR(prior, posterior, final_model, par)
% Adapted THBI prior/posterior comparison: conditioned counts and envelopes.
assert(~isempty(prior),'thbi:plot','No empirical prior was saved for this run.');
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1200 800]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for z=1:3
    nexttile;
    assert(isequal(prior.count_values{z},posterior.count_values{z}),'thbi:plot','Prior count support differs.');
    bar(posterior.count_values{z},[prior.count_probability{z}(:) posterior.count_probability{z}(:)]);
    xlabel('Finite physical layers'); ylabel('Probability'); title(par.mod.zone_names{z}); grid on;
    legend({'Conditioned prior','Retained posterior'},'Location','best');
end
for limit=[.3 2 max(final_model.Z)]
    nexttile; hold on;
    z=prior.depth_km(:); q=prior.vs_quantiles;
    fill([q(1,:) fliplr(q(3,:))],[z' fliplr(z')],[.87 .87 .87], ...
        'EdgeColor','none','DisplayName','Prior 5-95%');
    z=final_model.Z(:); q=final_model.VSquantiles';
    fill([q(1,:) fliplr(q(3,:))],[z' fliplr(z')],[.65 .8 .95], ...
        'EdgeColor','none','FaceAlpha',.6,'DisplayName','Posterior 5-95%');
    plot(q(2,:),z,'k-','DisplayName','Posterior median');
    set(gca,'YDir','reverse'); ylim([0 limit]); xlabel('Vs (km/s)'); ylabel('Depth (km)'); grid on;
    title(sprintf('Upper %g km',limit)); legend('Location','best');
end
exportgraphics(f,fullfile(par.res.dir,'prior_vs_posterior.png'),'Resolution',170);
end
