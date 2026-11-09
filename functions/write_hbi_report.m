function write_hbi_report(out,par,chains,s)
fid=fopen(fullfile(out,'inversion_report.txt'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'Fixed-layer five-observable HBI; %s\n',par.inv.profile);
fprintf(fid,'Runtime: %s (%s)\n',version,computer);
fprintf(fid,'%d chains x %d steps; burn-in %d; all repeated states retained.\n', ...
    par.inv.nchains,par.inv.niter,par.inv.burnin);
fprintf(fid,'%s\n%s\n',s.warning,s.rhat_kind);
fprintf(fid,'Maximum split R-hat: %.6g\n',max(s.split_rhat));
for c=1:numel(chains)
    q=chains{c};
    fprintf(fid,'Chain %d: start chi2 %.6g; best visited chi2 %.6g; acceptance %.4f; forward %d; failed %d; seconds %.2f\n', ...
        c,sum(q.initial.misfit.E2),sum(q.bestfit.misfit.E2),mean(q.accepted),q.nforward,q.failures,q.elapsed_seconds);
end
fprintf(fid,'Best fit uses minimum sum of C0-whitened energies, includes burn-in; it is not a posterior mean.\n');
for k=1:5
    fprintf(fid,'%s: base chi2/N %.6g; gamma at best fit %.6g\n', ...
        par.inv.datatypes{k},s.base_reduced_chi2(k),s.bestfit.model.datahparm.gamma(k));
end
fprintf(fid,'Parameter, q05, median, q95, classical_split_Rhat\n');
for j=1:numel(s.parameter_names)
    fprintf(fid,'%s, %.9g, %.9g, %.9g, %.9g\n',s.parameter_names{j},s.parameter_quantiles(:,j),s.split_rhat(j));
end
end
