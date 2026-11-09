function plot_RF_BANDS(forward, par, prefix)
% Plot full-band and short/intermediate/long-period views of each RF.
% All amplitudes retain their original units and sign. No peak normalization.
visibility='off'; if par.res.show_figures, visibility='on'; end
names={'rf','rf_sp'}; names=names(isfield(forward.pred,names));
for j=numel(names):-1:1
    name=names{j};
    if isfield(forward.pred.(name),'waveform')
        forward.pred.(name)=forward.pred.(name).waveform;
    end
    if ~isfield(forward.pred.(name),'diagnostic_period_bands_s'), names(j)=[]; end
end
if isempty(names), return; end
f=figure('Visible',visibility,'Color','w','Position',[100 100 1250 900]);
nrows=1+size(forward.pred.(names{1}).diagnostic_period_bands_s,1);
tiledlayout(nrows,numel(names),'TileSpacing','compact','Padding','compact');
for row=1:nrows
    for d=1:numel(names)
        name=names{d}; rf=forward.pred.(name); cfg=par.forc.(name);
        bands=[rf.period_band_s;rf.diagnostic_period_bands_s];
        assert(size(bands,1)==nrows,'thbi:plot','Ps/Sp diagnostic band counts must agree.');
        values=[rf.val,rf.band_val]; ax=nexttile;
        plot(ax,rf.t,values(:,row),'k-','LineWidth',1.2); grid(ax,'on');
        window=cfg.plot_window;
        if row==2 % Zoom short periods so their conversion arrivals are visible.
            window=[-10 40]; if strcmp(name,'rf_sp'), window=[-40 10]; end
        end
        xlim(ax,window); xline(ax,0,':');
        label='Ps R/Z'; if strcmp(name,'rf_sp'), label='Sp Z/R'; end
        if strcmp(name,'rf_sp') && isfield(par,'rf') && strcmp(par.rf.method,'hk_ccp')
            label='Sp P daughter before migration';
        end
        title(ax,sprintf('%s | %.3g-%.3g s periods',label,bands(row,:)));
        xlabel(ax,'Conversion lag (s)'); ylabel(ax,'Signed RF amplitude');
    end
end
sgtitle('RF frequency bands: diagnostic views of the same waveform');
exportgraphics(f,fullfile(par.res.dir,[prefix '_rf_bands.png']),'Resolution',180);
exportgraphics(f,fullfile(par.res.dir,[prefix '_rf_bands.pdf']),'ContentType','vector');
if ~par.res.show_figures, close(f); end
end
