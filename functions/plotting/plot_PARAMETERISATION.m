function plot_PARAMETERISATION(model, par)
% Physical-layer counterpart of Zach's spline-basis parameterization plot.
[N,K]=count_N_layers_in_mod(model,par);
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1000 700]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
boundaries=[0 model.zsed model.zmoh par.mod.maxz]; colors=lines(3);
for limit=[.3 max(par.res.depth_km)]
    ax=nexttile; hold(ax,'on');
    for z=1:3
        zone=model.(par.mod.zone_fields{z});
        edges=boundaries(z)+[0;zone.f(:);1]*diff(boundaries(z:z+1));
        xx=reshape([zone.VS(:) zone.VS(:)]',[],1);
        zz=reshape([edges(1:end-1) edges(2:end)]',[],1);
        plot(ax,xx,zz,'o-','Color',colors(z,:),'MarkerSize',3,'DisplayName',par.mod.zone_names{z});
    end
    plot(ax,[model.halfspaceVS model.halfspaceVS],[par.mod.maxz max(par.res.depth_km)], ...
        'k--','DisplayName','Half-space');
    yline(ax,model.zsed,':','Sediment base','HandleVisibility','off');
    yline(ax,model.zmoh,':','Moho','HandleVisibility','off');
    set(ax,'YDir','reverse'); ylim(ax,[0 limit]); grid(ax,'on');
    xlabel(ax,'Vs (km/s)'); ylabel(ax,'Depth (km)'); legend(ax,'Location','best');
    title(ax,sprintf('%d finite layers: K = [%d %d %d]',N,K));
end
exportgraphics(f,fullfile(par.res.dir,'physical_parameterization.png'),'Resolution',170);
end
