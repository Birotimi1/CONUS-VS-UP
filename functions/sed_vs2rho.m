% Imported from the local Zach/Brennan THBI_MCMC/functions source; see docs/THBI_REFERENCE_IMPORTS.json.
function rho  = sed_vs2rho( Vs )
% rho  = sed_vs2rho( Vs )
%   Empirical scaling of Vs to rho for sedimentary rocks. 
%   Equations from Shen and Ritzwoller (JGR, 2016) equation (2), based
%   on results of Brocher (BSSA 2005).

rho = 1.227 + 1.53*Vs - 0.837*Vs.^2 + 0.207*Vs.^3 - 0.0166*Vs.^4;

end

