function test_thbi_extended60()
% Run manually in MATLAB: RJ mixture, DR balance, RF bands and deconvolution.
% Pure MATLAB algebra only. Native CPS/HVf checks are a separate test.
par=bayes_inv_parms('extended60'); validate_configuration(par);
saved_rng=rng; cleanup=onCleanup(@() rng(saved_rng)); %#ok<NASGU>
rng(791,'twister');

%% =================== MIXTURE BIRTH DENSITY AND INVERSE MERGE ===================
par.mod.start.jitter=false; model=b1_INITIATE_MODEL(par);
par.proposal.probability=[.20 .10 .10 .25 .10 .25];
par.proposal.zone_probability=[.2 .3 .5];
for z=1:3
    edges=[0;model.(par.mod.zone_fields{z}).f;1];
    width=edges(2)-2*par.mod.min_fraction(z);
    eligible=sum(diff(edges)>2*par.mod.min_fraction(z));
    for delta=[.02 .24]
        request=struct('type','birth','zone',z,'index',1,'split_f',edges(2)*.4,'delta_vs',delta);
        [born,birth]=b2_PERTURB_MODEL(model,par,request);
        assert(birth.valid && birth.dr_eligible,'thbi:test','Mixture birth failed support.');
        w=par.proposal.birth_broad_probability;
        s1=par.proposal.birth_vs_sigma(z); s2=par.proposal.birth_vs_sigma_broad(z);
        density=(1-w)*exp(-delta^2/(2*s1^2))/(sqrt(2*pi)*s1) + ...
            w*exp(-delta^2/(2*s2^2))/(sqrt(2*pi)*s2);
        expected=.25*par.proposal.zone_probability(z)/eligible/width*density;
        assert(abs(birth.log_q_forward-log(expected))<1e-12);
        [back,death]=b2_PERTURB_MODEL(born,par,struct('type','death','zone',z,'index',1));
        assert(death.valid && death.dr_eligible && isequal(back.K,model.K));
        assert(max(abs(back.matrix(:)-model.matrix(:)))<1e-10);
        assert(abs(birth.log_hastings+death.log_hastings)<1e-10);
    end
end

%% =================== DELAYED-REJECTION PATHWISE DETAILED BALANCE ===================
% Compare the FULL x->y(rejected)->z path flux with its z->y->x reverse,
% including both proposal stages. This catches likelihood-only retries.
sig=.9; small=.25*sig; triples=[0 2 .4;-.7 1.8 .2;1 3 -.3];
for row=1:size(triples,1)
    x=triples(row,1); y=triples(row,2); z=triples(row,3);
    for outside=[false true]
        lx=-x^2/2; ly=-y^2/2; lz=-z^2/2;
        if outside, ly=-Inf; end % first candidate outside prior support
        sx=state(x,lx); sy=state(y,ly); sz=state(z,lz);
        p=proposal(x,y,sig); sz.ptb=proposal(x,z,small);
        [~,ab]=delayed_rejection_damp(sx,sy,sz,p);
        reverse=proposal(z,y,sig); sx.ptb=proposal(z,x,small);
        [~,ba]=delayed_rejection_damp(sz,sy,sx,reverse);
        forward_flux=lx+logphi(y-x,sig)+log1p(-min(1,exp(ly-lx))) + ...
            logphi(z-x,small)+min(0,ab);
        reverse_flux=lz+logphi(y-z,sig)+log1p(-min(1,exp(ly-lz))) + ...
            logphi(x-z,small)+min(0,ba);
        assert(abs(forward_flux-reverse_flux)<1e-12,'thbi:test','DR detailed balance failed.');
    end
end
sx=state(1,-.5); sy=state(0,0); sz=state(2,-2); sz.ptb=proposal(1,2,small);
[accepted,loga]=delayed_rejection_damp(sx,sy,sz,proposal(1,0,sig));
assert(~accepted && isinf(loga) && loga<0,'thbi:test','Certain first-stage acceptance should preclude DR.');
bad=proposal(1,0,sig); bad.type='birth';
expect_error(@() delayed_rejection_damp(sx,sy,sz,bad),'thbi:dr');
% NaN proposal factors must stop, never be turned into an acceptance by min.
current=struct('logL',0,'logP',0); trial=current;
bad=struct('valid',true,'log_q_reverse',NaN,'log_q_forward',0,'log_abs_jacobian',0);
expect_error(@() b6_IFACCEPT(current,trial,bad),'thbi:accept');
expect_error(@() validate_layer_model([0 -6 3.5 2.7]),'fiveobs:model');
% A NaN SAC start time used to evade abs(header-request)>tolerance checks.
sac_file=[tempname '.sac']; remove_sac=onCleanup(@() delete(sac_file)); %#ok<NASGU>
floats=zeros(70,1); floats([1 6 41 45])=[.1 NaN 1 .1];
integers=zeros(40,1); integers([7 10 16 36])=[6 64 1 1];
fid=fopen(sac_file,'wb','ieee-le'); assert(fid>=0);
fwrite(fid,floats,'float32'); fwrite(fid,integers,'int32');
fwrite(fid,zeros(192,1),'uint8'); fwrite(fid,zeros(64,1),'float32'); fclose(fid);
expect_error(@() read_cps_sac(sac_file),'fiveobs:sac');

%% =================== 30-60 S BAND RETAINS LONG, SUPPRESSES SHORT PERIODS ===================
cfg=par.forc.rf; t=((0:cfg.npts-1)'-cfg.npts/2)*cfg.dt;
trace=struct('val',sin(2*pi*t/45)+sin(2*pi*t/5),'dt',cfg.dt,'t',t);
filtered=receiver_function_bands(trace,cfg);
middle=abs(t)<90; basis=[sin(2*pi*t/45),cos(2*pi*t/45),sin(2*pi*t/5),cos(2*pi*t/5)];
fit=basis(middle,:)\filtered.band_val(middle,3);
long_gain=1/sqrt(1+(45/60)^8)/sqrt(1+(30/45)^8);
assert(abs(hypot(fit(1),fit(2))-long_gain)<.04 && hypot(fit(3),fit(4))<.01, ...
    'thbi:test','The 30-60 s diagnostic did not isolate long periods.');
assert(isequal(filtered.period_band_s,[1.25 60]) && size(filtered.band_val,2)==3);
short=trace; short.val=short.val(1:1000); short.t=short.t(1:1000);
expect_error(@() receiver_function_bands(short,cfg),'fiveobs:rf_config');

%% =================== SIGNED NEGATIVE-LAG Sp AND PHYSICAL AMPLITUDE ===================
cfg=par.forc.rf_sp; N=cfg.npts;
parent=struct('val',zeros(N,1),'dt',cfg.dt,'t',(0:N-1)'*cfg.dt-cfg.delay);
daughter=parent; parent.val(2001)=1; daughter.val(1881)=-.4;
rf=waterlevel_receiver_function(parent,daughter,cfg);
[peak,ix]=min(rf.val);
assert(abs(rf.t(ix)+12)<cfg.dt && abs(peak+.4*cfg.gaussian_alpha/sqrt(pi))<1e-5, ...
    'thbi:test','Sp lag, sign or continuous-time Gaussian amplitude is incorrect.');
daughter.val=2*daughter.val; twice=waterlevel_receiver_function(parent,daughter,cfg);
assert(max(abs(twice.val-2*rf.val))<1e-12,'thbi:test','RF amplitude was normalized.');
parent.val(:)=0;
expect_error(@() waterlevel_receiver_function(parent,daughter,cfg),'fiveobs:rf');
fprintf('PASS: mixture RJ densities, DR path balance, RF period bands and signed Sp deconvolution.\n');
end

function s=state(value,logdensity)
s=struct('model',struct('K',[3 3 2]),'logP',0,'logL',logdensity,'value',value);
end
function p=proposal(x,y,sigma)
p=struct('type','velocity','zone',1,'index',1,'dr_eligible',true, ...
    'value_from',x,'value_to',y,'sigma',sigma);
end
function val=logphi(x,sigma)
val=-.5*(x/sigma)^2-log(sqrt(2*pi)*sigma);
end
function expect_error(callback,identifier)
try
    callback(); error('thbi:test','Expected error was not raised.');
catch problem
    assert(strcmp(problem.identifier,identifier),'thbi:test','Unexpected error: %s',problem.message);
end
end
