function test_receiver_function_cps(cps_bin, test_root)
% Physical checks, not an inversion/recovery test. Requires real hrftn96.
if nargin < 2, test_root = tempname; end
mkdir(test_root);
cfg = struct('rayp_s_km', .06, 'gaussian_alpha', 2.5, ...
    'dt', .025, 'npts', 4096, 'delay', 10);

% A homogeneous medium has no delayed conversion or reverberation.
homogeneous = [0 6.3 3.6 2.7];
h = forward_receiver_function_cps(homogeneous, cfg, cps_bin, fullfile(test_root, 'homogeneous'));
late = h.t > 3 & h.t < 25;
assert(max(abs(h.val(late))) < 1e-4 * max(abs(h.val)), 'Spurious late RF energy in homogeneous medium.');

% Primary Moho Ps arrival must agree with the independently calculated lag.
crust = [35 6.3 3.6 2.7; 0 8.0 4.5 3.3];
r = forward_receiver_function_cps(crust, cfg, cps_bin, fullfile(test_root, 'moho35'));
expected = 35 * (sqrt(1/3.6^2 - .06^2) - sqrt(1/6.3^2 - .06^2));
window = find(abs(r.t - expected) < .6);
[peak, imax] = max(r.val(window));
measured = r.t(window(imax));
assert(peak > 0, 'Expected positive Ps conversion at a positive Moho velocity jump.');
assert(abs(measured - expected) <= .15, 'Moho Ps arrival is inconsistent with ray theory.');

% Artificial subdivision of a homogeneous layer must not change the RF.
split = [12 6.3 3.6 2.7; 23 6.3 3.6 2.7; 0 8.0 4.5 3.3];
s = forward_receiver_function_cps(split, cfg, cps_bin, fullfile(test_root, 'split_crust'));
assert(max(abs(s.val-r.val)) < 1e-4 * max(abs(r.val)), 'RF changed when an identical layer was subdivided.');

% Shared refinement must preserve layer properties at shallow/deep depths,
% including the final half-space and the physical interface at 35 km.
source = [crust(1,:); 10 8 4.5 3.3; crust(2,:)];
refined = refine_layer_model(source, 2);
depths = [0 .001 12 34.999 35.001 45 1000];
for depth = depths
    z0 = [0; cumsum(source(1:end-1,1))];
    z1 = [0; cumsum(refined(1:end-1,1))];
    row0 = find(z0 <= depth, 1, 'last');
    row1 = find(z1 <= depth, 1, 'last');
    assert(isequal(source(row0,2:4), refined(row1,2:4)), 'Numerical refinement changed the physical model.');
end
refined_rf = forward_receiver_function_cps(refined, cfg, cps_bin, fullfile(test_root, 'refined_crust'));
assert(max(abs(refined_rf.val-r.val)) < 1e-4 * max(abs(r.val)), 'Shared refinement changed the RF.');

% Reject units/configuration mistakes before invoking the solver.
bad = cfg; bad.rayp_s_km = 6.0;
failed = false;
try
    forward_receiver_function_cps(crust, bad, cps_bin, fullfile(test_root, 'invalid'));
catch problem
    failed = strcmp(problem.identifier, 'fiveobs:rayp');
end
assert(failed, 'Ray parameter in s/degree was not rejected as an invalid s/km value.');
fprintf('RF checks passed: homogeneous, Ps lag %.3f s (expected %.3f s), layer subdivision, ray-parameter validation.\n', measured, expected);
end
