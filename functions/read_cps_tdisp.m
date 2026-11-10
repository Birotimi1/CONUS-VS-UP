function out = read_cps_tdisp(filename)
% CPS330 tio.f gtsmdt/gtshed/gtsval: sequential unformatted TI dispersion.
% Read original double phase velocities; ASCII tables round to5significant
% digits and cannot distinguish a trapped mode from the half-space boundary.
fid=fopen(filename,'rb','ieee-le'); assert(fid>=0,'fiveobs:cps','Cannot open %s.',filename);
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
name=record(fid); ipar=record(fid); fpar=record(fid);
assert(numel(name)==80 && numel(ipar)==80 && numel(fpar)==80,'fiveobs:cps','Unknown CPS TI header.');
n=decode(record(fid),'int32'); assert(isscalar(n) && n>=1 && n<=200);
elastic=decode(record(fid),'single'); assert(numel(elastic)==7*double(n));
elastic=reshape(double(elastic),7,[]); nper=decode(record(fid),'int32');
assert(isscalar(nper) && nper>=1 && nper<=2049);
out=struct('period',zeros(0,1),'wave',zeros(0,1),'velocity',{{}}, ...
    'halfspace_shear',sqrt(elastic(5:6,end)/elastic(7,end))');
while true
    head=record(fid); if isempty(head), break; end
    assert(numel(head)==16,'fiveobs:cps','Malformed CPS dispersion header.');
    w=double(decode(head(1:4),'int32')); nm=double(decode(head(5:8),'int32'));
    T=decode(head(9:16),'double'); if w<0, break; end
    assert(any(w==[1 2]) && nm>=0 && nm<=2000 && isfinite(T) && T>0, ...
        'fiveobs:cps','Invalid CPS wave/mode count/period.');
    c=[];
    if nm>0
        c=decode(record(fid),'double');
        assert(numel(c)==nm && all(isfinite(c) & c>0),'fiveobs:cps','Invalid CPS roots.');
    end
    out.period(end+1,1)=T; out.wave(end+1,1)=w; out.velocity{end+1,1}=c(:);
end
assert(numel(out.period)==nper,'fiveobs:cps','Incomplete CPS frequency records.');
end

function bytes=record(fid)
n=fread(fid,1,'int32=>double');
if isempty(n), bytes=uint8([]); return; end
assert(n>=0 && n<1e8,'fiveobs:cps','Invalid Fortran record length.');
bytes=fread(fid,n,'uint8=>uint8'); tail=fread(fid,1,'int32=>double');
assert(numel(bytes)==n && isequal(n,tail),'fiveobs:cps','Truncated Fortran record.');
end

function value=decode(bytes,kind)
value=typecast(bytes(:),kind);
[~,~,endian]=computer;
if endian=='B', value=swapbytes(value); end
end
