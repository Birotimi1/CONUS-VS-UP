function quoted = quote_shell_path(value)
% POSIX shell quoting, including spaces, apostrophes, and metacharacters.
q = char(39);
quoted = [q strrep(char(value), q, [q char(34) q char(34) q]) q];
end
