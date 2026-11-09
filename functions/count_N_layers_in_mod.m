function [Nfinite, K, Nincluding_halfspace] = count_N_layers_in_mod(model, par)
% Count physical layers in ONE zoned model, not numerical solver subdivisions.
% Adapted role from THBI count_N_layers_in_mod; the original layerises spline
% profiles and leaves its output unassigned. This API returns exact counts.
if nargin<2
    fields={'sedmparm','crustmparm','mantmparm'};
else
    fields=par.mod.zone_fields;
end
K=zeros(1,3);
for z=1:3
    zone=model.(fields{z});
    K(z)=numel(zone.VS);
    assert(K(z)>=1 && numel(zone.f)==K(z)-1,'thbi:layers', ...
        'Zone %s must have K velocities and K-1 interfaces.',fields{z});
end
Nfinite=sum(K);
Nincluding_halfspace=Nfinite+1;
if isfield(model,'matrix')
    assert(size(model.matrix,1)==Nincluding_halfspace && model.matrix(end,1)==0, ...
        'thbi:layers','Physical matrix and zone layer counts disagree.');
end
end
