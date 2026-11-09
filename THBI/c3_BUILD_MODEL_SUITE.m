function suite_of_models = c3_BUILD_MODEL_SUITE(retained, par)
% THBI-style model suite evaluated as piecewise-constant physical layers.
% Unlike spline interpolation, step evaluation preserves sharp interfaces.

%% =================== COMMON DEPTH BASIS AND SAMPLE IDENTITIES ===================
suite_of_models.Z=par.res.depth_km(:);
suite_of_models.chain=retained.chain;
suite_of_models.iteration=retained.iteration;
suite_of_models.K=retained.K;
n=retained.Nstored; nz=numel(suite_of_models.Z);
suite_of_models.VS=zeros(nz,n);
suite_of_models.VP=zeros(nz,n);
suite_of_models.rho=zeros(nz,n);

%% =================== EVERY RETAINED EARTH ON THE SAME DEPTH GRID ===================
for i=1:n
    model=model_from_chain(retained,i,par);
    properties=sample_layer_model(model.matrix,suite_of_models.Z);
    suite_of_models.VP(:,i)=properties(:,1);
    suite_of_models.VS(:,i)=properties(:,2);
    suite_of_models.rho(:,i)=properties(:,3);
end
suite_of_models.note='Columns are retained iterations, including repeats; no averaging of variable layer indices.';
end
