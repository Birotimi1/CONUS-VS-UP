function rf = receiver_function_bands(rf, cfg)
% Explicit RF period passband and diagnostic bands; no duplicate likelihoods.
% Zero-phase fourth-order Butterworth amplitude factors, applied by padded FFT.
% Native Gaussian filtering is retained. Periods are seconds, not time lags.
if ~isfield(cfg,'period_band_s'), return; end
bands=cfg.period_band_s;
if isfield(cfg,'diagnostic_period_bands_s'), bands=[bands;cfg.diagnostic_period_bands_s]; end
assert(size(bands,2)==2 && isreal(bands) && all(isfinite(bands(:))) && ...
    all(bands(:,1)>0 & bands(:,2)>bands(:,1)), 'fiveobs:rf_config','Invalid RF period bands.');
rf.val=rf.val(:); N=numel(rf.val); dt=rf.dt;
assert(isreal(rf.val) && all(isfinite(rf.val)) && isfinite(dt) && dt>0, ...
    'fiveobs:rf','Invalid RF waveform or sampling interval.');
assert(N*dt>=4*max(bands(:,2)), 'fiveobs:rf_config', ...
    'RF record must span at least four cycles of the longest configured period.');
nfft=2^nextpow2(2*N); f=(0:nfft-1)'/(nfft*dt); f=min(f,1/dt-f);
spectrum=fft(rf.val,nfft); filtered=zeros(N,size(bands,1));
for j=1:size(bands,1)
    flo=1/bands(j,2); fhi=1/bands(j,1);
    assert(fhi<.5/dt,'fiveobs:rf_config','RF band exceeds Nyquist frequency.');
    highpass=zeros(size(f)); positive=f>0;
    highpass(positive)=1./sqrt(1+(flo./f(positive)).^8);
    lowpass=1./sqrt(1+(f/fhi).^8);
    y=real(ifft(spectrum.*highpass.*lowpass)); filtered(:,j)=y(1:N);
end
rf.unbanded_val=rf.val; rf.val=filtered(:,1);
rf.period_band_s=bands(1,:); rf.diagnostic_period_bands_s=bands(2:end,:);
rf.band_val=filtered(:,2:end);
rf.band_note='Diagnostic filtered views of the SAME waveform; not independent observations.';
end
