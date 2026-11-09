function report = test_thbi_reference_native(par,data,truth,folder)
% USER-RUN MATLAB native preflight. Calls CPS/HVf; not a convergence test.
% Inspect forward files and measured sensitivities before launching chains.
if ~exist(folder,'dir'), mkdir(folder); end
base=truth.zoned_parameters;
[born,p]=b2_PERTURB_MODEL(base,par,struct('type','birth','zone',2,'index',1, ...
    'split_f',base.crustmparm.f(1)/2,'delta_vs',.12));
assert(p.valid,'thbi:preflight','Native split fixture outside support.');
[back,p]=b2_PERTURB_MODEL(born,par,struct('type','death','zone',2,'index',1));
assert(p.valid && max(abs(back.matrix(:)-base.matrix(:)))<1e-10);
deep=base; deep.mantmparm.VS(end)=deep.mantmparm.VS(end)+.05;
deep=make_mod_from_parms(deep,par); assert(a1_TEST_CONDITIONS(deep,par));
kappa=base; kappa.vpvs=kappa.vpvs+.01; kappa=make_mod_from_parms(kappa,par);
models={base,born,back,deep,kappa}; predictions=cell(size(models));
if isfield(base,'xi')
    anis=base; anis.xi(2)=anis.xi(2)+.02; anis=make_mod_from_parms(anis,par);
    models{end+1}=anis; predictions=cell(size(models));
end
for j=1:numel(models)
    work=fullfile(folder,sprintf('model_%02d',j));
    try
        [predictions{j},fwd]=b3_FORWARD_MODEL(models{j},data,par,work);
        m=b4_CALC_MISFIT(data,predictions{j},par);
        L=b5_CALC_LIKELIHOOD(m,data,models{j}.datahparm,par);
        assert(isfinite(L) && all(isfinite(m.E2)));
        for d=1:numel(par.inv.datatypes)
            name=par.inv.datatypes{d}; y=predictions{j}.(name);
            assert(numel(y)==data.(name).N && isreal(y) && all(isfinite(y)));
            assert(isequal(fwd.pred.(name).model_input,fwd.model));
        end
    catch problem
        preserve_forward_failure(fullfile(folder,'failure'),models{j},data,par, ...
            problem,work,struct('stage','reference native preflight','iteration',0,'model_index',j));
        rethrow(problem);
    end
end
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; y=predictions{1}.(name);
    report.(name).roundtrip=max(abs(predictions{3}.(name)-y));
    report.(name).split_change=max(abs(predictions{2}.(name)-y));
    report.(name).deep_change=max(abs(predictions{4}.(name)-y));
    report.(name).kappa_change=max(abs(predictions{5}.(name)-y));
    if numel(models)>5, report.(name).anisotropy_change=max(abs(predictions{6}.(name)-y)); end
    assert(report.(name).roundtrip<1e-6*(1+max(abs(y))), ...
        'thbi:preflight','Native split/merge return mismatch for %s.',name);
    if d>=2 && d<=4
        ix=data.(name).x>=60;
        report.(name).long_period_deep_change=max(abs(predictions{4}.(name)(ix)-y(ix)));
    end
end
% HK uses crustal average Vs/H/kappa, so changing mantle Vs alone must NOT
% manufacture deep sensitivity. Other products' sensitivities are reported.
assert(report.rf.deep_change<1e-12,'thbi:preflight','HK acquired spurious mantle dependence.');
assert(report.cph.split_change>1e-8 && report.ugr.split_change>1e-8, ...
    'thbi:preflight','Dispersion did not respond to a nonzero crustal split.');
report.note='Execution/sensitivity checks only; 100-s coverage does not establish resolution at 250 km.';
save(fullfile(folder,'native_preflight.mat'),'models','predictions','report','-v7');
fid=fopen(fullfile(folder,'native_preflight.txt'),'w'); assert(fid>=0);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s\nObservable: split, deep-mantle, kappa change; roundtrip\n',report.note);
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; q=report.(name);
    fprintf(fid,'%s: %.9g %.9g %.9g; %.9g\n',name,q.split_change,q.deep_change,q.kappa_change,q.roundtrip);
end
end
