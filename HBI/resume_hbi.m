function resume_hbi(result_dir)
% Resume the SAME configuration and data, including RNG and proposal tuning.
% Deliberate changes to data/settings require a new run directory.
a0_STARTUP_BAYES();
load(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','runtime');
assert(strcmp(par.res.dir,result_dir),'Use the original absolute result directory.');
chains=run_all_chains(par,data);
summary=c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary);
fid=fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w');
assert(fid>=0,'Cannot write latest-run pointer.');
fprintf(fid,'%s\n',result_dir); fclose(fid);
end
