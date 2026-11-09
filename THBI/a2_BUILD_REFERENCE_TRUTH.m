function truth = a2_BUILD_REFERENCE_TRUTH(par)
% Published ECM1 column + explicitly labelled AK135 mantle, in layered form.
% Raw source properties are retained; matched synthetics use THBI materials.

%% =================== PUBLISHED SEDIMENT AND CRUST ===================
root=fullfile(par.forc.project_root,'DATA','reference_models');
raw=dlmread(fullfile(root,'ecm1_19178_layers.csv'),',',1,0);
ak=dlmread(fullfile(root,'ak135_mantle.csv'),',',1,0);
h=raw(:,2); vs=raw(:,4); zsed=sum(h(raw(:,1)==1)); moho=sum(h);
assert(abs(zsed-.01)<1e-10 && abs(moho-43.7731)<1e-8, ...
    'thbi:truth','Unexpected ECM1 cell/units; inspect the source manifest.');
assert(par.mod.maxz==250,'thbi:truth','This reference fixture is defined through 250 km.');
crust=raw(:,1)==2;
vpvs=sum(h(crust)./raw(crust,4))/sum(h(crust)./raw(crust,3));

%% =================== AK135 REFERENCE MANTLE TO 250 KM ===================
edges=[moho;77.5;120;165;210;250];
mid=(edges(1:end-1)+edges(2:end))/2;
mantle_vs=interp1(ak(:,1),ak(:,3),mid,'linear');
h=[h;diff(edges)]; vs=[vs;mantle_vs];
halfspace=interp1(ak(:,1),ak(:,3),250,'linear');
counts=[sum(raw(:,1)==1),sum(raw(:,1)==2),numel(mantle_vs)];
model=struct('zsed',zsed,'zmoh',moho,'halfspaceVS',halfspace,'vpvs',vpvs);
if isfield(par.mod,'xi_bounds'), model.xi=par.synth.xi; end
model.datahparm.eta=par.mod.noise.log_center;
zone_edges=[0 zsed moho par.mod.maxz]; first=1;
for z=1:3
    rows=first:first+counts(z)-1;
    f=cumsum(h(rows(1:end-1)))/diff(zone_edges(z:z+1));
    model.(par.mod.zone_fields{z})=struct('f',f(:),'VS',vs(rows));
    first=rows(end)+1;
end
[ok,why]=a1_TEST_CONDITIONS(model,par);
assert(ok,'thbi:truth','Published reference structure outside prior: %s',why);
model=make_mod_from_parms(model,par); validate_layer_model(model.matrix);
assert(abs(sum(model.matrix(:,1))-250)<1e-8 && ...
    max(abs(model.matrix(:,3)-[vs;halfspace]))<1e-12, ...
    'thbi:truth','Layered assembly changed the reference thickness/Vs.');

%% =================== SEPARATE SOURCE VALUES FROM MATCHED MATERIAL CLOSURE ===================
truth=struct('case_name','ecm1_19178_ak135_thbi_materials', ...
    'description','ECM1 cell 19178 (10 m sediment, 43.7731 km Moho) + AK135 reference mantle to 250 km; THBI Vp/rho closure', ...
    'model_matrix',model.matrix,'theta',[h;vs;halfspace], ...
    'sed_thickness',zsed,'sediment_base_depth',zsed, ...
    'moho_depth',moho,'moho_interface',sum(counts(1:2)), ...
    'total_depth',250,'zone_counts',counts,'zoned_parameters',model, ...
    'noise_added',false,'noise_type','none', ...
    'published_ecm1_layers',raw,'ak135_velocity_nodes',ak, ...
    'latitude',36.5,'longitude',-82.5,'ecm1_cell',19178, ...
    'source_note','ECM1 thickness/Vs retained. Vp/rho projected to THBI material relations for a matched synthetic. Mantle Vs is a midpoint-layer approximation to AK135; no inferred or imposed LVZ.');
truth.source_manifest=jsondecode(fileread(fullfile(root,'ecm1_19178_provenance.json')));
if isfield(model,'xi')
    truth.assigned_xi=model.xi;
    truth.description=[truth.description '; assigned crustal radial anisotropy xi=1.06'];
    truth.source_note=[truth.source_note ' Anisotropy is an assigned synthetic test, not measured ECM1 or AK135 anisotropy.'];
end
end
