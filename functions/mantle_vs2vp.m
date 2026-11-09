% Imported from the local Zach/Brennan THBI_MCMC/functions source; see docs/THBI_REFERENCE_IMPORTS.json.
function Vp  = mantle_vs2vp( Vs,Zkm )
% Vp  = mantle_vs2vp( Vs,Zkm )
%   Scaling of Vs to Vp for mantle rocks, using Vp/Vs ratio at each depth
%   from AK135.


% akmod = ak135('depths',Zkm,'crust',false);

try
%     Vp = Vs.*(akmod.vp./akmod.vs);
    Vp = Vs.*ak135vpvs(Zkm);
catch e
    Vp=1.81*Vs;
    error(getReport(e))
end

end

function vpvs = ak135vpvs(Zkm)


% values are taken from ak135 function with continental crust
% http://rses.anu.edu.au/seismology/ak135/ak135f.html
% columns are Z, vpvs_nocrust
% a = [
%      0.000      1.6763
%     20.000      1.6763
%     20.0001     1.6883
%     35.000      1.6883
%     35.0001     1.7946
%     77.500      1.7918
%    120.000      1.7889
%    165.000      1.8130
%    210.000      1.8371
%    210.0001     1.8351
%    260.000      1.8404
%    310.000      1.8452
%    360.000      1.8498
%    410.000      1.8425];
%%% brb2022.07.22 Whoops, at some point i put crustal vpvs ratios in here.
%%% Zach had already removed them. Now I've taken them back out...
a = [
    0           1.7946 % Extend mantle vpvs ratio from the moho to the surface. This function is for MANTLE vpvs. If we have a shallow Moho, its vpvs should still be MANTLE vpvs and not CRUST vpvs. 
    35.000      1.7946
    77.500      1.7918
   120.000      1.7889
   165.000      1.8130
   210.000      1.8371 % Remove discontinuities. We don't want those from ak135, we would want any such discontinuities to come straight from our model. 
   260.000      1.8404
   310.000      1.8452
   360.000      1.8498
   410.000      1.8425];

vpvs = interp1(a(:,1),a(:,2),Zkm, 'linear'); % Cubic interpolation gives smooth results, so it doesn't artificially introduce a layer that our modelling didn't insert. 

% figure(1); clf; hold on; plot(vpvs)

end
