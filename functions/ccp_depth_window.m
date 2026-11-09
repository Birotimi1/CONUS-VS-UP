function weight = ccp_depth_window(z, limits, taper_km)
% Continuous cosine version of Ge Jin/Z. Eilon flat_hanning_win.
% limits give the half-amplitude window edges; no Signal Toolbox dependency.
z=z(:); left=limits(1)-taper_km/2; right=limits(2)+taper_km/2;
assert(taper_km>0 && diff(limits)>taper_km,'thbi:ccp','Invalid CCP window.');
weight=ones(size(z)); weight(z<left | z>right)=0;
i=z>=left & z<left+taper_km; weight(i)=.5*(1-cos(pi*(z(i)-left)/taper_km));
i=z>right-taper_km & z<=right; weight(i)=.5*(1-cos(pi*(right-z(i))/taper_km));
end
