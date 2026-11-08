function plot_TRADEOFF(posterior, par)
% Adapted THBI tradeoff panels for quantities shared across different K.
% No layer-index correlation or sampled Vp/Vs/anisotropy is implied.
x={posterior.zsed,posterior.zmoh,posterior.zmoh,posterior.gamma(:,1)};
y={posterior.zone_mean_vs(:,1),posterior.zone_mean_vs(:,2),posterior.fdVSmoh,posterior.gamma(:,5)};
xl={'Sediment base (km)','Moho depth (km)','Moho depth (km)','HVSR noise multiplier'};
yl={'Thickness-weighted sediment Vs (km/s)','Thickness-weighted crust Vs (km/s)', ...
    'Moho Vs jump (% of interface mean)','Ps RF noise multiplier'};
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1000 800]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
display_rows=unique(round(linspace(1,posterior.Nstored,min(10000,posterior.Nstored))));
for k=1:4
    ax=nexttile;
    scatter(ax,x{k}(display_rows),y{k}(display_rows),5,posterior.chain(display_rows),'filled');
    xlabel(ax,xl{k}); ylabel(ax,yl{k}); grid(ax,'on');
    dx=x{k}-mean(x{k}); dy=y{k}-mean(y{k}); denom=sqrt(sum(dx.^2)*sum(dy.^2));
    r=NaN; if denom>0, r=sum(dx.*dy)/denom; end
    title(ax,sprintf('All-state correlation r = %.3f; color = chain',r));
end
exportgraphics(f,fullfile(par.res.dir,'parameter_tradeoffs.png'),'Resolution',170);
end
