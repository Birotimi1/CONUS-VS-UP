function result_dir = run_synthetic_hbi(profile, overrides)
% New synthetic run; resume an existing one with resume_hbi(result_dir).
% Pass a configured par structure, or a profile plus optional inv overrides.
% Example: run_synthetic_hbi('pilot',struct('parallel',true)).
if nargin<1, profile='pilot'; end
if nargin<2, overrides=struct(); end
a0_STARTUP_BAYES();
if isstruct(profile)
    par=profile; profile=par.inv.profile;
else
    par=a1_DEFINE_PARAMETERS(profile);
end
allowed={'parallel','niter','burnin','nchains'};
names=fieldnames(overrides);
for k=1:numel(names)
    assert(any(strcmp(names{k},allowed)),'Unsupported run override.');
    par.inv.(names{k})=overrides.(names{k});
end
assert(par.inv.niter>par.inv.burnin && par.inv.burnin>=0 && ...
    par.inv.nchains>=1 && all(mod([par.inv.niter par.inv.burnin par.inv.nchains],1)==0), ...
    'Invalid chain lengths.');
if ~exist(par.res.root,'dir'), mkdir(par.res.root); end
result_dir=fullfile(par.res.root,['ECM1_thin_HBI_' profile '_' datestr(now,'yyyymmdd_HHMMSS')]);
if exist(result_dir,'dir'), error('hbi:exists','Result directory already exists.'); end
mkdir(result_dir); par.res.dir=result_dir;
fprintf('Fixed-layer five-observable HBI: %s\n%s\n',profile,result_dir);
[data,TRUEmodel,forward]=a2_LOAD_DATA(par);
par.forc=forward.config; % freeze resolved executable paths/settings for every chain
runtime.version=version; runtime.computer=computer;
save(fullfile(result_dir,'inputs.mat'),'par','data','TRUEmodel','forward','runtime','-v7');
chains=run_all_chains(par,data);
summary=c1_PROCESS_RESULTS(chains,data,par);
save(fullfile(result_dir,'results.mat'),'par','data','TRUEmodel','chains','summary','runtime','-v7');
c2_PLOT_RESULTS(par,data,TRUEmodel,chains,summary);
fid=fopen(fullfile(par.res.root,'LATEST_RUN.txt'),'w');
assert(fid>=0,'Cannot write latest-run pointer.');
fprintf(fid,'%s\n',result_dir); fclose(fid);
fprintf('Completed. Review fits AND diagnostics: %s\n',result_dir);
end
