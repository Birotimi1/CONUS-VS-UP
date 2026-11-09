function preserve_forward_failure(destination, model, data, par, problem, workdir, context)
% Preserve the first failure at this destination before temporary files vanish.
% Saving diagnostics must never replace the original MATLAB exception.
if exist(destination,'dir'), return; end
try
    mkdir(destination);
    if exist(workdir,'dir')
        copyfile(workdir,fullfile(destination,'native_files'));
    end
    failure.identifier=problem.identifier;
    failure.message=problem.message;
    failure.report=getReport(problem,'extended','hyperlinks','off');
    failure.stack=problem.stack;
    failure.context=context;
    failure.rng=rng;
    save(fullfile(destination,'failure.mat'),'model','data','par','failure','-v7');
    fid=fopen(fullfile(destination,'failure.txt'),'w');
    if fid>=0
        fprintf(fid,'%s\n',failure.report);
        fclose(fid);
    end
catch saving_problem
    warning('thbi:failure_save','Could not preserve failure in %s: %s',destination,saving_problem.message);
end
end
