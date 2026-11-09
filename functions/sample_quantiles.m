function q = sample_quantiles(x, p)
% Empirical quantiles along rows, with linear interpolation (no toolbox).
x=sort(x,1); n=size(x,1); q=zeros(numel(p),size(x,2));
for k=1:numel(p)
    r=1+(n-1)*p(k); lo=floor(r); hi=ceil(r);
    q(k,:)=x(lo,:)+(r-lo)*(x(hi,:)-x(lo,:));
end
end
