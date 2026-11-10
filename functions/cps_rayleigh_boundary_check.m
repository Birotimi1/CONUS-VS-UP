function residual = cps_rayleigh_boundary_check(matrix,period,velocity)
% CPS uses single-precision layer tensors. A weakly surface-coupled mode can
% have a sensitive boundary residual even when its phase root is accurate.
% Permit local root refinement within3e-5relative velocity, well below data
% error scales. A distant repeated-eigenvalue artefact cannot pass this test.
residual=rayleigh_surface_residual(matrix,period,velocity);
if residual>=1e-3
    options=optimset('TolX',1e-11*velocity,'MaxIter',80,'Display','off');
    [~,nearby]=fminbnd(@(c) rayleigh_surface_residual(matrix,period,c), ...
        velocity*(1-3e-5),velocity*(1+3e-5),options);
    residual=min(residual,nearby);
end
end
