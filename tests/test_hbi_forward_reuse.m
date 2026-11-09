function test_hbi_forward_reuse()
% Detect stale native solver files when one chain reuses its scratch folder.
repo=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(repo,'HBI')); a0_STARTUP_BAYES();
par=bayes_inv_parms('pilot'); par.res.dir=tempname; mkdir(par.res.dir);
cleanup=onCleanup(@() rmdir(par.res.dir,'s')); %#ok<NASGU>
[data,truth,forward]=z1_SYNTH_DATA(par); par.forc=forward.config;
scratch=fullfile(par.res.dir,'reused'); mkdir(scratch);
model=make_layer_model(truth.theta,par);
a=b3_FORWARD_MODEL(model,data,par,scratch);
t=truth.theta; t(14)=t(14)+.08;
b=b3_FORWARD_MODEL(make_layer_model(t,par),data,par,scratch);
c=b3_FORWARD_MODEL(model,data,par,scratch);
for k=1:5
    d=par.inv.datatypes{k};
    assert(max(abs(a.(d)-c.(d)))<1e-10,'Stale or nondeterministic solver output.');
    fprintf('%s maximum change: %.9g\n',d,max(abs(a.(d)-b.(d))));
    assert(any(abs(a.(d)-b.(d))>1e-7),'%s did not respond to changed structure.',d);
    assert(max(abs(a.(d)-data.(d).clean))<1e-10,'Generation/inversion forward mismatch.');
end
fprintf('PASS: true -> perturbed -> true in reused scratch, all five native solvers.\n');
end
