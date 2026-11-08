function par = bayes_inv_parms(profile)
% Compatibility name. Edit HBI/a1_DEFINE_PARAMETERS.m for the actual settings.
if nargin == 0, profile = 'pilot'; end
addpath(fileparts(fileparts(mfilename('fullpath'))));
par = a1_DEFINE_PARAMETERS(profile);
end
