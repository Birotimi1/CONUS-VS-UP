function snapshot = source_snapshot()
% Exact source text guards checkpoint continuation against edited algorithms.
% Snapshot THBI stages/parameters/tests and the shared project helpers.
% Restrict the roots explicitly: generated DATA files must never enter this
% guard after the helper moves from THBI/functions to repository/functions.
% Uses MATLAB file IO only; no compiler, Java or external hashing dependency.
root=fileparts(fileparts(mfilename('fullpath')));
files=[dir(fullfile(root,'THBI','**','*.m')); ...
    dir(fullfile(root,'functions','**','*.m'))];
paths=arrayfun(@(f) fullfile(f.folder,f.name),files,'UniformOutput',false);
[paths,~]=sort(paths);
snapshot=cell(numel(paths),2);
for j=1:numel(paths)
    snapshot{j,1}=paths{j}(numel(root)+2:end);
    snapshot{j,2}=fileread(paths{j});
end
end
