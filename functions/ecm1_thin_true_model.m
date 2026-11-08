function tm = ecm1_thin_true_model()
% Prescribed synthetic truth from the existing HBI ECM1 thin test, not a
% measured/recovered real-site model. h in km and Vs in km/s.
h = [0.010; 0.040; 0.080; 0.870; 8; 14; 17; 8; 22];
vs = [0.20; 0.55; 1.50; 2.80; 3.55; 3.75; 3.98; 4.30; 4.45; 4.54];
% Preserve the original HBI Brocher conversion and floor order exactly.
vp = 0.9409 + 2.0947*vs - 0.8206*vs.^2 + 0.2683*vs.^3 - 0.0251*vs.^4;
rho = 1.6612*vp - 0.4721*vp.^2 + 0.0671*vp.^3 - 0.0043*vp.^4 + 0.000106*vp.^5;
vp = max(vp, 1.5*vs);
rho = max(rho, 1.5);
model = [[h; 0], vp, vs, rho];
validate_layer_model(model);
assert(abs(sum(h)-70) < 1e-10 && abs(sum(h(1:7))-40) < 1e-10);
tm = struct('case_name', 'ecm1_thin', 'description', ...
    'ECM1 thin synthetic truth: 50 m soft cap, 40 km Moho, 70 km to half-space', ...
    'model_matrix', model, 'theta', [h; vs], 'sed_thickness', .05, ...
    'moho_depth', 40, 'moho_interface', 7, 'total_depth', 70, ...
    'noise_added', false, 'noise_type', 'none', ...
    'source_note', 'Assigned h and Vs from HBI ECM1 thin; Vp/density derived using the same empirical equations.');
end
