function plot_PERTURBATIONS(summary, par)
% THBI-style attempted/prior-valid/accepted move diagnostics after burn-in.
visibility='off'; if par.res.show_figures, visibility='on'; end
f=figure('Visible',visibility,'Color','w','Position',[80 80 1100 700]);
if ~par.res.show_figures, cleanup=onCleanup(@() close(f)); end %#ok<NASGU>
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
nexttile;
has_second=isfield(summary,'dr_move_attempts');
values=[sum(summary.move_attempts,1)',sum(summary.move_accepts,1)'];
labels={'First attempted','First accepted'};
if has_second
    values=[values,sum(summary.dr_move_attempts,1)',sum(summary.dr_move_accepts,1)'];
    labels=[labels,{'Second attempted','Second path accepted'}];
end
bar(values);
xticks(1:6); xticklabels(par.proposal.names); xtickangle(25);
ylabel('Retained-iteration proposals'); legend(labels,'Location','best'); grid on;
for z=1:3
    nexttile;
    attempts=sum(summary.zone_move_attempts(:,4:5,z),1);
    accepts=sum(summary.zone_move_accepts(:,4:5,z),1);
    values=[attempts(:) accepts(:)];
    if has_second
        da=sum(summary.dr_zone_move_attempts(:,4:5,z),1);
        dc=sum(summary.dr_zone_move_accepts(:,4:5,z),1);
        values(:,3:4)=[da(:),dc(:)];
    end
    bar(values);
    xticks(1:2); xticklabels({'Birth','Death'}); ylabel('Proposals'); grid on;
    title(par.mod.zone_names{z}); legend(labels,'Location','best');
end
exportgraphics(f,fullfile(par.res.dir,'proposal_acceptance.png'),'Resolution',170);
end
