function chains = run_all_chains(par, data)
% Independent process workers, or serial fallback explicitly selected in par.
root=fileparts(mfilename('fullpath'));
chains = cell(par.inv.nchains,1);
if par.inv.parallel
    if exist('OCTAVE_VERSION','builtin') || ~license('test','Distrib_Computing_Toolbox')
        error('hbi:parallel','Parallel mode requires MATLAB Parallel Computing Toolbox.');
    end
    pool = gcp('nocreate');
    if isempty(pool)
        cluster = parcluster('local');
        cluster.JobStorageLocation = fullfile(par.res.dir,'parallel_jobs');
        if ~exist(cluster.JobStorageLocation,'dir'), mkdir(cluster.JobStorageLocation); end
        pool = parpool(cluster,min(par.inv.nchains,cluster.NumWorkers)); %#ok<NASGU>
    elseif isa(pool,'parallel.ThreadPool')
        error('hbi:parallel','Use process workers: native solvers change working directories.');
    end
    parfor k = 1:par.inv.nchains
        addpath(root);
        a0_STARTUP_BAYES(); % select HBI stages and shared helpers on each worker
        chains{k} = run_one_chain(par,data,k);
    end
else
    for k = 1:par.inv.nchains, chains{k} = run_one_chain(par,data,k); end
end
end
