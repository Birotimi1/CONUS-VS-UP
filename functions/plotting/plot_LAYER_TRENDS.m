function plot_LAYER_TRENDS(chains, par)
% Layer-interface counterpart of THBI plot_KNOT_TRENDS; no spline knots.
% Render at most 500 states per chain; count traces use every saved iteration.
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1200 700]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(2,numel(chains),'TileSpacing','compact','Padding','compact');
for c=1:numel(chains)
    chain=chains{c}; n=chain.last_iter;
    selected=unique(round(linspace(1,n,min(n,500))));
    interfaces=NaN(sum(par.mod.kmax)-1,numel(selected));
    for j=1:numel(selected)
        model=model_from_chain(chain,selected(j),par);
        depth=cumsum(model.matrix(1:end-2,1)); % excludes fixed half-space boundary
        interfaces(1:numel(depth),j)=depth;
    end
    ax=nexttile(c); hold(ax,'on');
    plot(ax,repmat(selected,size(interfaces,1),1),interfaces,'.','Color',[.4 .4 .4]);
    plot(ax,1:n,chain.boundaries(1:n,1),'b-',1:n,chain.boundaries(1:n,2),'r-');
    xline(ax,par.inv.burnin,'k:'); set(ax,'YDir','reverse'); ylim(ax,[0 par.mod.maxz]);
    xlabel(ax,'Iteration'); ylabel(ax,'Interface depth (km)'); title(ax,sprintf('Chain %d: blue basement, red Moho',c));
    ax=nexttile(numel(chains)+c);
    stairs(ax,1:n,chain.K(1:n,:)); xline(ax,par.inv.burnin,'k:');
    xlabel(ax,'Iteration'); ylabel(ax,'Physical layer count'); grid(ax,'on');
    legend(ax,par.mod.zone_names,'Location','best');
end
exportgraphics(f,fullfile(par.res.dir,'layer_interface_trends.png'),'Resolution',170);
end
