function logP = calc_log_prior(model, par)
% Fixed-layer compatibility name; the full prior stays in HBI/a3.
% Never change the active algorithm's search path inside a shared helper.
repo = fileparts(fileparts(mfilename('fullpath')));
for name = {'a3_CALC_LOG_PRIOR','a1_TEST_CONDITIONS'}
    expected = fullfile(repo,'HBI',[name{1} '.m']);
    assert(strcmp(which(name{1}),expected),'hbi:path', ...
        'Run the HBI startup before calling calc_log_prior; expected %s.',expected);
end
logP = a3_CALC_LOG_PRIOR(model,par);
end
