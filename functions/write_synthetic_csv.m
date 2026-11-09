function write_synthetic_csv(filename, header, values)
fid = fopen(filename, 'w');
if fid < 0, error('fiveobs:io', 'Cannot write %s.', filename); end
closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s\n', header);
fprintf(fid, [repmat('%.12g,',1,size(values,2)-1) '%.12g\n'], values');
end
