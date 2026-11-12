% VALIDATE_THBI_BODY_WAVES -- independent anisotropic body-wave benchmark.
% Local MATLAB only. No scheduler, data generation or inversion is launched.

%% =================== LOCAL CONFIGURATION AND PRESERVED EVIDENCE ===================
repo=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(repo,'THBI'));
a0_STARTUP_BAYES(); par=bayes_inv_parms('reference100_vti');
par.forc.verbose=true; validate_configuration(par);
folder=fullfile(par.res.root,['BODY_REFERENCE_' datestr(now,'yyyymmdd_HHMMSS')]);
assert(~exist(folder,'dir'),'thbi:audit','Audit directory already exists.');
mkdir(folder); diary(fullfile(folder,'validation.log'));
copyfile([mfilename('fullpath') '.m'],fullfile(folder,'validation_launcher.m'));
fprintf('Independent VTI body-wave validation: %s\n',folder);

%% =================== ANALYTIC, LAYERED AND NATIVE BODY-WAVE CHECKS ===================
% The reference uses physical displacement/traction ODEs, generic eigensystems
% and bounded impedance recursion. It never calls HVf's polarization or Wang
% routines. Analytic half-spaces and a separate expm recurrence check it first.
body_report=test_thbi_body_wave_reference(par,folder,true);
assert(body_report.passed,'thbi:preflight','Independent body-wave validation failed.');
quadrature_report=test_thbi_body_quadrature(par,body_report,fullfile(folder,'quadrature'));
assert(quadrature_report.passed,'thbi:preflight','Body-wave quadrature validation failed.');

%% =================== COMPONENT-BY-COMPONENT SUMMARY ===================
for j=1:numel(body_report.cases)
    c=body_report.cases{j};
    fprintf('%s: max relative [PSV-horizontal SH-horizontal PSV-vertical] = %s\n', ...
        c.name,mat2str(max(c.relative_error,[],1),5));
end
fprintf('PASS. Saved reference code, matrices, native logs and results in %s\n',folder);
fprintf('This tests forward physics on specified fixtures, not inversion convergence or depth resolution.\n');
diary off;
