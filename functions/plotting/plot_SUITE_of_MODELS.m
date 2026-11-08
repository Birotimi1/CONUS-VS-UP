function plot_SUITE_of_MODELS(suite, final_model, par)
% THBI-style model-suite display; deterministic display subset by chain.
% The subset controls rendering only. All retained states enter summaries.
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1000 700]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ids=unique(suite.chain); colors=lines(numel(ids));
for field={'VS','VP'}
    name=field{1}; ax=nexttile; hold(ax,'on');
    for c=1:numel(ids)
        selected=find(suite.chain==ids(c));
        selected=selected(unique(round(linspace(1,numel(selected),min(200,numel(selected))))));
        plot(ax,suite.(name)(:,selected),suite.Z,'Color',.7+.3*colors(c,:),'LineWidth',.3);
    end
    plot(ax,final_model.([name 'av']),suite.Z,'k-','LineWidth',1.8);
    set(ax,'YDir','reverse'); grid(ax,'on'); ylim(ax,[0 max(suite.Z)]);
    xlabel(ax,[name ' (km/s)']); ylabel(ax,'Depth (km)');
    title(ax,'Up to 200 retained profiles per chain; black = median');
end
exportgraphics(f,fullfile(par.res.dir,'model_suite.png'),'Resolution',170);
end
