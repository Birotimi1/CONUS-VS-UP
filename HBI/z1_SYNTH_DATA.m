function [data, TRUEmodel, forward] = z1_SYNTH_DATA(par)
% Older entry point retained for existing callers/tests.
% The complete true-model/data/covariance/noise code is in a2_LOAD_DATA.m.
[data,TRUEmodel,forward] = a2_LOAD_DATA(par);
end
