function save_checkpoint(filename, checkpoint)
% Write then rename within the same directory; never partially overwrite.
temporary = [filename '.writing.mat'];
save(temporary,'checkpoint','-v7');
[ok,message] = movefile(temporary,filename,'f');
if ~ok, error('hbi:checkpoint','%s',message); end
end
