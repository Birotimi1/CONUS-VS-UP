function write_synthetic_text(filename, value)
fid = fopen(filename, 'w');
if fid < 0, error('fiveobs:io', 'Cannot write %s.', filename); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', value);
end
