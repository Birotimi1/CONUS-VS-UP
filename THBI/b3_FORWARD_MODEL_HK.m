function hk = b3_FORWARD_MODEL_HK(model, observed, par)
% Ps H-kappa energy on OBSERVED receiver functions, recomputed for this Earth.
% Layered isotropic adaptation of Brennan's hk_forward_model/HKstack wrapper.
% The sampled crustal Vp/Vs is not inferred by fitting a synthetic R/Z vector.

%% =================== OBSERVED WAVES AND CRUSTAL AVERAGE ===================
cfg=par.rf.hk; waves=observed.waves;
t=waves.t(:); rf=waves.rf; rayp=waves.rayp_s_km(:);
assert(size(rf,1)==numel(t) && size(rf,2)==numel(rayp) && ...
    all(isfinite(rf(:))) && all(diff(t)>0),'thbi:hk','Invalid observed Ps waves.');
h=model.matrix(1:end-1,1); top=[0;cumsum(h(1:end-1))];
hc=max(0,min(h,model.zmoh-top)); crust=hc>0;
assert(abs(sum(hc)-model.zmoh)<1e-8,'thbi:hk','Crustal integration does not reach Moho.');
% Harmonic velocity average, as in Brennan hk_average_crust. We keep the
% current FLAT CPS geometry; spherical/flat equivalence remains an open item.
vs=model.matrix(1:end-1,3); vsav=sum(hc)/sum(hc(crust)./vs(crust));
H=cfg.H_km(:)'; K=cfg.kappa(:); weights=cfg.phase_weights/sum(cfg.phase_weights);
assert(model.zmoh>=min(H) && model.zmoh<=max(H) && ...
    model.vpvs>=min(K) && model.vpvs<=max(K),'thbi:hk','Trial H/kappa outside stack grid.');

%% =================== Ps, PpPs AND PpSs+PsPs STACK ===================
energy=zeros(numel(K),numel(H));
for j=1:numel(rayp)
    p=rayp(j); qs=sqrt(vsav^-2-p^2); qp=sqrt((vsav*K).^-2-p^2);
    assert(isreal(qs) && isreal(qp),'thbi:hk','Nonpropagating HK phase.');
    times={(qs-qp)*H,(qs+qp)*H,repmat(2*qs*H,numel(K),1)};
    for phase=1:3
        assert(min(times{phase}(:))>=t(1) && max(times{phase}(:))<=t(end), ...
            'thbi:hk','Observed Ps window must include all HK grid phase arrivals.');
        sign_phase=1; if phase==3, sign_phase=-1; end
        energy=energy+sign_phase*weights(phase)*interp1(t,rf(:,j),times{phase},'linear');
    end
end
% Match reference hk_forward_model: interpolate the GRID energy at the
% candidate, then divide by that same grid's maximum (no candidate clipping).
trial_energy=interp2(H,K,energy,model.zmoh,model.vpvs,'linear');
maximum=max(energy(:));
assert(isfinite(maximum) && maximum>0,'thbi:hk','Ps stack has no positive energy.');
hk=struct('H_km',H,'kappa',K,'energy',energy,'trial_energy',trial_energy, ...
    'maximum_energy',maximum,'E_by_Emax',trial_energy/maximum,'vs_average',vsav, ...
    'trial_H_km',model.zmoh,'trial_kappa',model.vpvs, ...
    'kind','Ps H-kappa stack; isotropic, harmonic crustal Vs; weights [0.5 0.3 0.2]');
end
