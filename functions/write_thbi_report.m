function write_thbi_report(par, chains, s)
fid=fopen(fullfile(par.res.dir,'inversion_report.txt'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Zoned transdimensional five-observable HBI; %s\n',par.inv.profile);
fprintf(fid,'MATLAB/runtime %s (%s)\n',version,computer);
fprintf(fid,'%d chains x %d steps; burn-in %d; all repeated states retained.\n',par.inv.nchains,par.inv.niter,par.inv.burnin);
fprintf(fid,'%s\n%s\n',s.warning,s.rhat_kind);
if isfield(par.proposal,'dr') && isfield(par.proposal.dr,'mode')
    fprintf(fid,'Delayed rejection: enabled=%d; mode=%s\n',par.proposal.dr.enabled,par.proposal.dr.mode);
end
for c=1:numel(chains)
    q=chains{c};
    fprintf(fid,'Chain %d: start fit energy %.6g; best visited fit energy %.6g; forward %d; failed %d; seconds %.2f\n', ...
        c,model_fit_score(q.initial.misfit.E2',par),model_fit_score(q.bestfit.misfit.E2',par),q.nforward,q.failures,q.elapsed_seconds);
    if isfield(q,'bestfit_retained') && ~isempty(q.bestfit_retained)
        fprintf(fid,'  Best retained fit energy %.6g at iteration %d\n', ...
            model_fit_score(q.bestfit_retained.misfit.E2',par),q.bestfit_retained_iteration);
    end
    for j=1:numel(par.proposal.names)
        fprintf(fid,'  %s retained FIRST-stage moves: %d accepted / %d prior-valid / %d attempted\n', ...
            par.proposal.names{j},s.move_accepts(c,j),s.move_valid(c,j),s.move_attempts(c,j));
    end
    for z=1:3
        fprintf(fid,'  %s FIRST-stage birth accepted/attempted %d/%d; death accepted/attempted %d/%d\n', ...
            par.mod.zone_names{z},s.zone_move_accepts(c,4,z),s.zone_move_attempts(c,4,z), ...
            s.zone_move_accepts(c,5,z),s.zone_move_attempts(c,5,z));
    end
    fprintf(fid,'  Retained iterations with native forward failures: %d\n',s.retained_forward_failures(c));
    if isfield(s,'dr_attempts')
        fprintf(fid,'  Retained second-stage DR: %d accepted / %d attempted\n',s.dr_accepts(c),s.dr_attempts(c));
        if isfield(s,'dr_move_attempts')
            for j=1:numel(par.proposal.names)
                fprintf(fid,'    Second proposal %s: %d accepted paths / %d prior-valid / %d attempted\n', ...
                    par.proposal.names{j},s.dr_move_accepts(c,j),s.dr_move_valid(c,j),s.dr_move_attempts(c,j));
            end
        end
    end
end
fprintf(fid,'%d saved response draws per chain. %s\n',s.nresponse_per_chain,s.response_note);
if any(s.retained_forward_failures)
    fprintf(fid,'Inspect failed models before interpreting the posterior: numerical failures were rejected, which can remove valid-prior regions.\n');
end
fprintf(fid,'%s\n',s.bestfit_selection);
if isfield(s,'bestfit_iteration') && isfinite(s.bestfit_iteration)
    fprintf(fid,'Fit-plot/final selection: chain %d iteration %d\n',s.bestfit_chain,s.bestfit_iteration);
end
for d=1:numel(par.inv.datatypes)
    fprintf(fid,'%s base residual energy/effective DOF: %.8g\n',par.inv.datatypes{d},s.base_reduced_chi2(d));
end
for z=1:3
    fprintf(fid,'%s K, retained probability\n',par.mod.zone_names{z});
    fprintf(fid,'%d, %.8g\n',[s.count_values{z};s.count_probability{z}]);
end
fprintf(fid,'Metric, q05, median, q95, classical_split_Rhat\n');
for j=1:numel(s.metric_names)
    fprintf(fid,'%s, %.9g, %.9g, %.9g, %.9g\n',s.metric_names{j},s.metric_quantiles(:,j),s.split_rhat(j));
end
if isfield(s,'rank_diagnostics')
    fprintf(fid,'Metric, rank_folded_Rhat, bulk_ESS, tail_ESS\n');
    for j=1:numel(s.metric_names)
        fprintf(fid,'%s, %.9g, %.9g, %.9g\n',s.metric_names{j},s.rank_diagnostics.rhat(j), ...
            s.rank_diagnostics.bulk_ess(j),s.rank_diagnostics.tail_ess(j));
    end
end
fprintf(fid,'Depth_km, Vs_q05, Vs_median, Vs_q95, classical_split_Rhat\n');
fprintf(fid,'%.9g, %.9g, %.9g, %.9g, %.9g\n', ...
    [s.depth_km(:)';s.vs_quantiles;s.vs_split_rhat(:)']);
end
