function logq = birth_contrast_logpdf(delta, par, zone)
% Exact density of the fixed narrow/broad Gaussian birth-contrast mixture.
% Use this SAME marginal density for births and their reverse deaths.
sigma=par.proposal.birth_vs_sigma(zone);
w=0;
if isfield(par.proposal,'birth_broad_probability'), w=par.proposal.birth_broad_probability; end
a=-.5*(delta/sigma)^2-log(sigma)-.5*log(2*pi);
if w==0, logq=a; return; end
broad=par.proposal.birth_vs_sigma_broad(zone);
b=-.5*(delta/broad)^2-log(broad)-.5*log(2*pi);
terms=[log1p(-w)+a,log(w)+b]; m=max(terms);
logq=m+log(sum(exp(terms-m)));
end
