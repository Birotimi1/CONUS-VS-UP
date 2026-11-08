function test_thbi_forward()
% Optional native-solver check: one common model before/after a birth/death.
% Run manually after the algebraic checks. This DOES call CPS and HVf.
par=bayes_inv_parms('pilot'); par.mod.start.jitter=false;
root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
par.res.dir=root;
[data,~,forward]=a2_LOAD_DATA(par); par.forc=forward.config;
model=b1_INITIATE_MODEL(par);
scratch=fullfile(root,'reused'); mkdir(scratch);
[a,fa]=b3_FORWARD_MODEL(model,data,par,scratch);
[born,birth]=b2_PERTURB_MODEL(model,par,struct('type','birth','zone',2,'index',1,'split_f',.12,'delta_vs',.08));
assert(birth.valid);
[b,fb]=b3_FORWARD_MODEL(born,data,par,scratch);
[back,death]=b2_PERTURB_MODEL(born,par,struct('type','death','zone',2,'index',1));
assert(death.valid);
[c,fc]=b3_FORWARD_MODEL(back,data,par,scratch);
for d=1:5
    name=par.inv.datatypes{d};
    assert(all(isfinite(b.(name))) && numel(b.(name))==data.(name).N);
    assert(max(abs(a.(name)-b.(name)))>1e-8,'Prediction did not respond to a birth: %s',name);
    assert(max(abs(a.(name)-c.(name)))<1e-6,'Native round-trip mismatch: %s',name);
    assert(isequal(fa.pred.(name).model_input,fa.model) && ...
        isequal(fb.pred.(name).model_input,fb.model) && ...
        isequal(fc.pred.(name).model_input,fc.model),'Solvers used different Earths.');
end
fprintf('PASS: five native forward solvers accept a variable-layer model and recover predictions after merge.\n');
end
