function plot_ANISOTROPY(chains,truth,par)
% Explicit THBI-average vs Vsv/Vsh convention and retained anisotropy traces.
if ~isfield(par.mod,'xi_bounds'), return; end
out=fullfile(par.res.dir,'figures'); if ~exist(out,'dir'), mkdir(out); end
visible='off'; if par.res.show_figures, visible='on'; end
f=figure('Visible',visible,'Color','w','Position',[100 100 1150 450]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
nexttile; hold on;
for c=1:numel(chains), plot(chains{c}.xi(:,2),'DisplayName',sprintf('Chain %d',c)); end
xline(par.inv.burnin,'k:'); yline(truth.zoned_parameters.xi(2),'k--','True');
xlabel('Iteration'); ylabel('Crustal xi = (Vsh/Vsv)^2'); grid on; legend('show');
nexttile; hold on;
keep=par.inv.burnin+1:par.inv.niter;
for c=1:numel(chains), histogram(chains{c}.xi(keep,2),'Normalization','pdf','DisplayStyle','stairs'); end
xline(truth.zoned_parameters.xi(2),'k--'); xlabel('Crustal xi'); ylabel('Retained density'); grid on;
nexttile; hold on; matrix=truth.model_matrix; e=vti_elastic(matrix);
plot_layer_model(matrix,'k',1.5,'THBI average Vs');
sv=matrix; sv(:,3)=e.Vsv; plot_layer_model(sv(:,1:4),[0 .35 .8],1.2,'Vsv');
sh=matrix; sh(:,3)=e.Vsh; plot_layer_model(sh(:,1:4),[.8 .1 .1],1.2,'Vsh');
set(gca,'YDir','reverse'); ylim([0 80]); xlabel('Velocity (km/s)'); ylabel('Depth (km)');
grid on; legend('show'); title('Assigned synthetic anisotropy');
exportgraphics(f,fullfile(out,'anisotropy.png'),'Resolution',180);
if ~par.res.show_figures, close(f); end
end
