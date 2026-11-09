function rf = waterlevel_receiver_function(parent, daughter, cfg)
% Regularized daughter/parent RF from CPS component waveforms.
% For Sp: parent=R, daughter=Z; signed Z/R, negative conversion lags.
% No time reversal, polarity flip or peak-amplitude normalization.
P=parent.val(:); D=daughter.val(:); N=numel(P); dt=cfg.dt;
assert(N==cfg.npts && numel(D)==N && isreal(P) && isreal(D) && ...
    all(isfinite([P;D])) && numel(parent.t)==N && numel(daughter.t)==N && ...
    max(abs(parent.t(:)-daughter.t(:)))<2e-5 && ...
    abs(parent.dt-dt)<1e-6 && abs(daughter.dt-dt)<1e-6, ...
    'fiveobs:rf','Incompatible parent/daughter waveforms.');
assert(cfg.water_level>0 && cfg.water_level<1 && ...
    cfg.edge_taper_fraction>=0 && cfg.edge_taper_fraction<.5, ...
    'fiveobs:rf_config','Invalid water level or edge taper.');
nt=floor(cfg.edge_taper_fraction*N); taper=ones(N,1);
if nt>1
    edge=.5*(1-cos(pi*(0:nt-1)'/(nt-1)));
    taper(1:nt)=edge; taper(end-nt+1:end)=flipud(edge);
end
nfft=2^nextpow2(2*N); pf=fft(P.*taper,nfft); df=fft(D.*taper,nfft);
power=abs(pf).^2; maximum=max(power);
assert(isfinite(maximum) && maximum>0,'fiveobs:rf','Zero parent spectrum.');
denominator=max(power,cfg.water_level*maximum);
% Signed FFT frequencies ensure a real delayed output for noninteger shifts.
f=[0:nfft/2,(-nfft/2+1):-1]'/(nfft*dt);
gaussian=exp(-(pi*f/cfg.gaussian_alpha).^2); gaussian(nfft/2+1)=0;
shift=exp(-2i*pi*f*cfg.delay);
response=real(ifft(df.*conj(pf)./denominator.*gaussian.*shift))/dt;
rf=struct('val',response(1:N),'dt',dt,'t',(0:N-1)'*dt-cfg.delay, ...
    'gaussian_alpha',cfg.gaussian_alpha,'rayp_s_km',cfg.rayp_s_km, ...
    'water_level',cfg.water_level,'edge_taper_fraction',cfg.edge_taper_fraction);
rf.kind='Synthetic Sp RF: regularized signed Z/R, negative precursor lags';
end
