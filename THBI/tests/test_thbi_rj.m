function test_thbi_rj()
% Deterministic transform/density checks. Run manually; no native solvers.
par=bayes_inv_parms('pilot'); par.mod.start.jitter=false;
validate_configuration(par); model=b1_INITIATE_MODEL(par);

%% =================== SPLIT -> MERGE AND MERGE -> SPLIT ===================
for z=1:3
    field=par.mod.zone_fields{z}; e=[0;model.(field).f;1];
    s=(e(1)+e(2))/2;
    request=struct('type','birth','zone',z,'index',1,'split_f',s,'delta_vs',.02);
    [born,birth]=b2_PERTURB_MODEL(model,par,request);
    assert(birth.valid && born.K(z)==model.K(z)+1 && born.M==model.M+2);
    [back,death]=b2_PERTURB_MODEL(born,par,struct('type','death','zone',z,'index',1));
    assert(death.valid && isequal(back.K,model.K));
    assert(max(abs(back.matrix(:)-model.matrix(:)))<1e-10,'Split/merge did not recover the Earth.');
    assert(abs(birth.log_hastings+death.log_hastings)<1e-10,'Proposal density ratios are not inverse.');
    old=struct('logL',-30,'logP',a3_CALC_LOG_PRIOR(model,par));
    new=struct('logL',-29.3,'logP',a3_CALC_LOG_PRIOR(born,par));
    [~,ab]=b6_IFACCEPT(old,new,birth); [~,ad]=b6_IFACCEPT(new,old,death);
    assert(abs(ab+ad)<1e-10,'Full RJ log acceptance ratios are not inverse.');

    [merged,death]=b2_PERTURB_MODEL(model,par,struct('type','death','zone',z,'index',1));
    assert(death.valid);
    request=struct('type','birth','zone',z,'index',1,'split_f',death.split_f,'delta_vs',death.delta_vs);
    [restored,birth]=b2_PERTURB_MODEL(merged,par,request);
    assert(birth.valid && max(abs(restored.matrix(:)-model.matrix(:)))<1e-10);
    assert(abs(birth.log_hastings+death.log_hastings)<1e-10);
end

%% =================== INDEPENDENT NUMERICAL JACOBIAN ===================
% Differentiate (v,s,delta)->(v_left,s,v_right), not an implementation helper.
a=.17; b=.83; x=[3.5;.38;.09]; J=zeros(3);
for j=1:3
    xp=x; xm=x; xp(j)=xp(j)+1e-6; xm(j)=xm(j)-1e-6;
    rp=(xp(2)-a)/(b-a); rm=(xm(2)-a)/(b-a);
    yp=[xp(1)-(1-rp)*xp(3);xp(2);xp(1)+rp*xp(3)];
    ym=[xm(1)-(1-rm)*xm(3);xm(2);xm(1)+rm*xm(3)];
    J(:,j)=(yp-ym)/(2e-6);
end
assert(abs(abs(det(J))-1)<1e-7,'Split Jacobian is not unity.');

%% =================== ABSOLUTE PROPOSAL DENSITIES, NOT JUST RECIPROCITY ===================
unequal=par;
unequal.proposal.probability=[.20 .10 .10 .25 .10 .25];
unequal.proposal.zone_probability=[.2 .3 .5];
sparse=model; sparse.crustmparm.f=[.0015;.5];
sparse=make_mod_from_parms(sparse,unequal);
assert(a1_TEST_CONDITIONS(sparse,unequal));
% The first of three crustal layers is too thin to split: E=2, not K=3.
request=struct('type','birth','zone',2,'index',2,'split_f',.2,'delta_vs',.04);
[born,move]=b2_PERTURB_MODEL(sparse,unequal,request);
assert(move.valid);
epsilon=unequal.mod.min_fraction(2); sigma=unequal.proposal.birth_vs_sigma(2);
q_forward=.25*.3/2/(.5-.0015-2*epsilon) * ...
    exp(-.5*(.04/sigma)^2)/(sqrt(2*pi)*sigma);
q_reverse=.10*.3/3;
assert(abs(move.log_q_forward-log(q_forward))<1e-12);
assert(abs(move.log_q_reverse-log(q_reverse))<1e-12);
% Independent prior ratio for K=3 -> 4 in the unordered crustal zone.
prior_ratio=(3/4)*(6/(1-4*epsilon)^3)/(2/(1-3*epsilon)^2) / ...
    (unequal.mod.vsmax(2)-unequal.mod.vsmin(2));
assert(abs(a3_CALC_LOG_PRIOR(born,unequal)-a3_CALC_LOG_PRIOR(sparse,unequal)-log(prior_ratio))<1e-11);
old=struct('logL',-5,'logP',a3_CALC_LOG_PRIOR(sparse,unequal));
new=struct('logL',-3,'logP',a3_CALC_LOG_PRIOR(born,unequal));
[~,logalpha]=b6_IFACCEPT(old,new,move);
assert(abs(logalpha-(2+log(prior_ratio)+log(q_reverse/q_forward)))<1e-11);

%% =================== COUNT LIMITS ARE SELF-TRANSITIONS ===================
limited=par; limited.mod.kmax=model.K;
[~,move]=b2_PERTURB_MODEL(model,limited,struct('type','birth','zone',1));
assert(~move.valid && strcmp(move.type,'birth'));
limited=par; limited.mod.kmin=model.K;
[~,move]=b2_PERTURB_MODEL(model,limited,struct('type','death','zone',2));
assert(~move.valid && strcmp(move.type,'death'));
limited=par; limited.mod.min_fraction(1)=.24; limited.mod.kmax(1)=4;
tight=model; tight.sedmparm.f=[1/3;2/3]; tight=make_mod_from_parms(tight,limited);
assert(a1_TEST_CONDITIONS(tight,limited));
[~,move]=b2_PERTURB_MODEL(tight,limited,struct('type','birth','zone',1));
assert(~move.valid && ~isempty(strfind(move.reason,'No layer')));

%% =================== ORDERED-VELOCITY NORMALIZATION AND ASSEMBLY ===================
ordered_par=par; ordered_par.conditions.monotonic(1)=true;
[ordered,~]=a3_CALC_LOG_PRIOR(model,ordered_par);
unordered=par; unordered.conditions.monotonic(1)=false;
[ordinary,~]=a3_CALC_LOG_PRIOR(model,unordered);
assert(abs(ordered-ordinary-log(factorial(model.K(1))))<1e-11,'Missing ordered-uniform K! density.');
assert(abs(sum(model.matrix(:,1))-par.mod.maxz)<1e-10);
assert(model.matrix(end,1)==0 && size(model.matrix,1)==sum(model.K)+1);

%% =================== MOHO CONVENTION FROM THE REFERENCE THBI ===================
check=model; check.crustmparm.VS(end)=3; check.mantmparm.VS(1)=3.85;
assert(~a1_TEST_CONDITIONS(check,par),'Absolute 0.8 km/s Moho limit was omitted.');
relaxed=par; relaxed.conditions.max_moho_jump=Inf;
assert(a1_TEST_CONDITIONS(check,relaxed),'Mean-normalized Moho fraction was not used.');
check.mantmparm.VS(1)=3.9;
assert(~a1_TEST_CONDITIONS(check,relaxed),'25 percent fractional Moho limit was omitted.');
fprintf('PASS: RJ round trips, absolute densities, reciprocal ratios, Jacobian, count bounds, priors, assembly and Moho limits.\n');
end
