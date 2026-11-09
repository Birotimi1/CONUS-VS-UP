function diagnostic = mcmc_rank_diagnostics(x)
% Rank-normalized/folded split R-hat and bulk/tail ESS, without toolboxes.
% x: iterations x fixed metrics x chains. Retain repeated states.
% Equations: Vehtari et al. (2021), doi:10.1214/20-BA1221;
% https://mc-stan.org/posterior/reference/rhat.html and ess_bulk.html.
% Geyer positive, monotone paired-autocorrelation truncation; FFT covariance.
n=floor(size(x,1)/2); nm=size(x,2); nc=size(x,3);
diagnostic=struct('rhat',NaN(1,nm),'bulk_ess',NaN(1,nm),'tail_ess',NaN(1,nm));
if n<4 || nc<2, return; end
for k=1:nm
    original=reshape(x(:,k,:),size(x,1),nc);
    if any(~isfinite(original(:))), continue; end
    y=[original(1:n,:),original(end-n+1:end,:)];
    ranked=normal_scores(y);
    folded=normal_scores(abs(y-median(original(:))));
    r=[rhat(ranked),rhat(folded)];
    % An exactly constant folded variable has no scale information. Preserve
    % the bulk diagnostic, while an entirely constant metric stays undefined.
    diagnostic.rhat(k)=max(r);
    diagnostic.bulk_ess(k)=ess(ranked);
    q=sample_quantiles(original(:),[.05 .95]);
    tails=[ess(double(y<=q(1))),ess(double(y>=q(2)))];
    if all(isfinite(tails)), diagnostic.tail_ess(k)=min(tails); end
end
end

function z=normal_scores(y)
[v,order]=sort(y(:)); ranks=zeros(size(v)); starts=[1;find(diff(v)~=0)+1];
ends=[starts(2:end)-1;numel(v)];
for j=1:numel(starts), ranks(starts(j):ends(j))=(starts(j)+ends(j))/2; end
u=zeros(size(v)); u(order)=(ranks-3/8)/(numel(v)+1/4);
z=reshape(sqrt(2)*erfinv(2*u-1),size(y));
end

function value=rhat(y)
n=size(y,1); W=mean(var(y,0,1)); between=var(mean(y,1),0,2);
if W==0
    value=NaN; if between>0, value=Inf; end
else
    value=sqrt(((n-1)/n*W+between)/W);
end
end

function value=ess(y)
[n,m]=size(y); W=mean(var(y,0,1)); value=NaN;
if W<=0, return; end
variance=(n-1)/n*W+var(mean(y,1),0,2);
y=bsxfun(@minus,y,mean(y,1));
spectrum=fft(y,2^nextpow2(2*n),1);
acov=real(ifft(abs(spectrum).^2,[],1)); acov=acov(1:n,:)/n;
rho=1-(W-mean(acov,2))/variance; rho(1)=1;
pairs=rho(1:2:2*floor(n/2))+rho(2:2:2*floor(n/2));
stop=find(pairs<=0,1); if ~isempty(stop), pairs=pairs(1:stop-1); end
for j=2:numel(pairs), pairs(j)=min(pairs(j),pairs(j-1)); end
tau=max(1/log10(n*m),-1+2*sum(pairs));
value=n*m/tau;
end
