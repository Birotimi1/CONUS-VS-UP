function probability = plot_HEATMAP_ALLMODELS(suite, final_model, par, depth_limit)
% THBI-style velocity distributions on the physical common-depth grid.
% Histograms use all retained states, not just accepted or best-fit models.
if nargin<4, depth_limit=max(suite.Z); end
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1000 700]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
probability=struct();
for field={'VS','VP'}
    name=field{1}; values=suite.(name);
    lo=min(values(:)); hi=max(values(:)); padding=max(.02,.01*(hi-lo));
    edges=linspace(lo-padding,hi+padding,121); centers=(edges(1:end-1)+edges(2:end))/2;
    mass=zeros(numel(suite.Z),numel(centers));
    for j=1:numel(suite.Z)
        mass(j,:)=histcounts(values(j,:),edges)/size(values,2);
    end
    probability.(name)=struct('edges',edges,'depth_km',suite.Z,'mass',mass);
    ax=nexttile; hold(ax,'on');
    surf(ax,centers,suite.Z,log10(max(mass,1e-4)),'EdgeColor','none'); view(ax,2);
    plot3(ax,final_model.([name 'av']),suite.Z,ones(size(suite.Z)), ...
        'k-','LineWidth',1.5);
    set(ax,'YDir','reverse'); ylim(ax,[0 depth_limit]); xlim(ax,[edges(1) edges(end)]);
    caxis(ax,[-4 0]); colormap(ax,parula); cb=colorbar(ax);
    ylabel(cb,'log10(probability per velocity bin)');
    xlabel(ax,[name ' (km/s)']); ylabel(ax,'Depth (km)');
    title(ax,[name ': all retained models; black = pointwise median']);
end
exportgraphics(f,fullfile(par.res.dir,sprintf('posterior_heatmap_%gkm.png',depth_limit)),'Resolution',170);
end
