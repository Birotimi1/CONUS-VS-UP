function hk = b3_FORWARD_MODEL_HK(model, observed, par)
% Ps H-kappa energy on OBSERVED receiver functions, recomputed for this Earth.
% Layered adaptation of Brennan's hk_forward_model/HKstack wrapper.
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
% Reference hk_average_crust uses arithmetic thickness means for density
% and anisotropy, harmonic Vs. Keep the same averages including sediments.
elastic=vti_elastic(model.matrix); rows=find(crust);
xiav=sum(hc(rows).*elastic.xi(rows))/sum(hc);
phiav=sum(hc(rows).*elastic.phi(rows))/sum(hc);
etaav=sum(hc(rows).*elastic.eta(rows))/sum(hc);
rhoav=sum(hc(rows).*elastic.rho(rows))/sum(hc);
average=[zeros(size(K)),vsav*K,vsav*ones(size(K)),rhoav*ones(size(K)), ...
    xiav*ones(size(K)),phiav*ones(size(K)),etaav*ones(size(K))];
assert(model.zmoh>=min(H) && model.zmoh<=max(H) && ...
    model.vpvs>=min(K) && model.vpvs<=max(K),'thbi:hk','Trial H/kappa outside stack grid.');

%% =================== Ps, PpPs AND PpSs+PsPs STACK ===================
energy=zeros(numel(K),numel(H));
for j=1:numel(rayp)
    [qp,qs]=vti_vertical_slowness(average,rayp(j));
    times={(qs-qp)*H,(qs+qp)*H,2*qs*H};
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
    'xi_average',xiav,'phi_average',phiav,'eta_average',etaav, ...
    'kind','Ps H-kappa stack; flat VTI, harmonic crustal Vs; weights [0.5 0.3 0.2]');
end
