function model = make_layer_model(theta, par, eta)
% Fixed-layer compatibility name; edit the full material-property block in HBI/b1.
% Never change the active algorithm's search path inside a shared helper.
if nargin < 3, eta = par.mod.noise.log_center; end
repo = fileparts(fileparts(mfilename('fullpath')));
expected = fullfile(repo,'HBI','b1_INITIATE_MODEL.m');
assert(strcmp(which('b1_INITIATE_MODEL'),expected),'hbi:path', ...
    'Run the HBI startup before calling make_layer_model; expected %s.',expected);
model = b1_INITIATE_MODEL(par,theta,eta);
end
