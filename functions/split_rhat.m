function r = split_rhat(x)
% Classical split R-hat for fixed-dimensional summaries, not layer-index values.
% Input: retained iterations x metrics x chains. This is not rank-normalized.
n=floor(size(x,1)/2); m=size(x,3);
if n<2 || m<2, r=NaN(1,size(x,2)); return; end
y=cat(3,x(1:n,:,:),x(end-n+1:end,:,:));
W=mean(var(y,0,1),3); B=n*var(mean(y,1),0,3);
r=sqrt(((n-1)/n*W+B/n)./W);
r(W==0 & B>0)=Inf; % different constant chains: no within-chain mixing
r(W==0 & B==0)=NaN; % identical constants: variance ratio is undefined
end
