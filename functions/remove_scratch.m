function remove_scratch(workdir)
if exist(workdir,'dir'), rmdir(workdir,'s'); end
end
