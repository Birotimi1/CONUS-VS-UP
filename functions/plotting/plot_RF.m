function [h,hline,hpos,hneg] = plot_RF(ax, RF, tt, level, xshift, iffill, pcol, ncol)
% THBI-style vertical RF wiggle, with explicit axes and no external maxab.
% Adapted from THBI plotting/plot_RF.m; finite, increasing-time input only.
if nargin<4 || isempty(level), level=0; end
if nargin<5 || isempty(xshift), xshift=0; end
if nargin<6 || isempty(iffill), iffill=true; end
if nargin<7 || isempty(pcol), pcol=[.8 .15 .1]; end
if nargin<8 || isempty(ncol), ncol=[.1 .3 .8]; end
RF=RF(:); tt=tt(:);
assert(numel(RF)==numel(tt) && numel(tt)>1 && all(diff(tt)>0) && ...
    all(isfinite(RF)) && all(isfinite(tt)),'thbi:plot','RF must be finite with strictly increasing times.');
assert(isscalar(level) && isfinite(level) && level>=0,'thbi:plot','Invalid RF fill threshold.');
hold(ax,'on'); hpos=[]; hneg=[];
if iffill
    positive=RF; positive(positive<level)=0;
    negative=RF; negative(negative>-level)=0;
    hpos=patch(ax,[xshift;positive+xshift;xshift],[tt(1);tt;tt(end)],pcol, ...
        'EdgeColor','none','FaceAlpha',.5);
    hneg=patch(ax,[xshift;negative+xshift;xshift],[tt(1);tt;tt(end)],ncol, ...
        'EdgeColor','none','FaceAlpha',.5);
end
hline=plot(ax,RF+xshift,tt,'k-','LineWidth',1.2);
set(ax,'YDir','reverse'); xlabel(ax,'R/Z amplitude'); ylabel(ax,'Time after P (s)');
h={hline,hpos,hneg};
end
