function resume_thbi(result_dir)
% Continue the same RJ run, with unchanged code, data and settings.
a0_STARTUP_BAYES();
load(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','runtime');
assert(strcmp(par.res.dir,result_dir),'Use the original absolute result directory.');
chains=run_all_chains(par,data);
summary=c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');
[posterior,retained]=c2_BUILD_POSTERIOR(chains,par);
suite_of_models=c3_BUILD_MODEL_SUITE(retained,par);
final_model=c4_FINAL_MODEL(posterior,suite_of_models,retained,par);
save(fullfile(result_dir,'postprocessing.mat'),'posterior','suite_of_models','final_model','-v7.3');
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary,posterior,suite_of_models,final_model);
fid=fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w'); assert(fid>=0);
fprintf(fid,'%s\n',result_dir); fclose(fid);
end
