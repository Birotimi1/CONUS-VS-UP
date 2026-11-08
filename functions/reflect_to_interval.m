function x = reflect_to_interval(x, lo, hi)
% Symmetric reflection of a Gaussian random walk at fixed conditional bounds.
if lo==hi, x=lo; return; end
assert(isfinite(lo) && isfinite(hi) && lo<hi,'thbi:bounds','Invalid interval.');
w=hi-lo; q=mod(x-lo,2*w);
if q>w, q=2*w-q; end
x=min(hi,max(lo,lo+q));
end
