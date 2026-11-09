function plot_progress_chain(chain, par, ichain)
% Inspect saved complete or partial history; no solver or sampler execution.
% Adapted THBI progress view for our acceptance and forward-failure fields.
n=chain.last_iter; it=(1:n)';
assert(n>=1,'thbi:plot','No completed iterations to plot.');
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1100 750]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile; plot(it,chain.logL(1:n)); xlabel('Iteration'); ylabel('log likelihood'); grid on;
nexttile; plot(it,cumsum(chain.accepted(1:n))./it); xlabel('Iteration'); ylabel('Cumulative acceptance fraction'); ylim([0 1]); grid on;
nexttile; stairs(it,cumsum(chain.forward_failed(1:n))); xlabel('Iteration'); ylabel('Cumulative forward failures'); grid on;
nexttile; semilogy(it,exp(chain.eta(1:n,:))); xlabel('Iteration'); ylabel('Noise multiplier'); legend(par.inv.datatypes,'Location','best'); grid on;
sgtitle(sprintf('Chain %d: %d completed iterations',ichain,n));
exportgraphics(f,fullfile(par.res.dir,sprintf('progress_chain_%02d.png',ichain)),'Resolution',170);
end
