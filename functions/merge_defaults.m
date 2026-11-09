function a = merge_defaults(a, b)
names = fieldnames(b);
for k = 1:numel(names)
    if ~isfield(a, names{k}), a.(names{k}) = b.(names{k}); end
end
end
