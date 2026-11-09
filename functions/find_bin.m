function directory = find_bin(explicit, env_name, candidates, executable)
directory = explicit;
if isempty(directory), directory = getenv(env_name); end
if ~isempty(directory)
    if exist(fullfile(directory, executable), 'file'), return; end
    error('fiveobs:dependency', '%s is not in configured directory %s.', executable, directory);
end
% Prefer project-local binaries, then PATH, then the supplied shared project.
search = [candidates(1) strsplit(getenv('PATH'), ':') candidates(2:end)];
for k = 1:numel(search)
    if ~isempty(search{k}) && exist(fullfile(search{k}, executable), 'file')
        directory = search{k}; return;
    end
end
error('fiveobs:dependency', 'Cannot find %s. Set %s or its cfg bin directory.', executable, env_name);
end
