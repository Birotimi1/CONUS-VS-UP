function signature = check_native_vti_build(cfg)
% Validate compiled HVf source/binary pair, then record actual solver hashes.
% Call once before native preflight/chain start, not in the proposal loop.
root=fileparts(cfg.hvf_vti_bin); previous=pwd;
cleanup=onCleanup(@() cd(previous)); %#ok<NASGU>
assert(exist(fullfile(root,'bin','source.sha256'),'file')==2 && ...
    exist(fullfile(root,'bin','executable.sha256'),'file')==2, ...
    'thbi:build','Build native/HVf_VTI/build.sh before the native VTI test.');
cd(root);
[has_sha256sum,~]=system('command -v sha256sum >/dev/null 2>&1');
if has_sha256sum==0, hash_command='sha256sum'; else, hash_command='shasum -a 256'; end
[status,message]=system([hash_command ' -c bin/source.sha256 bin/executable.sha256 2>&1']);
assert(status==0,'thbi:build','HVf binary/source mismatch: rebuild native/HVf_VTI/build.sh. %s',message);
programs={'surf96','hrftn96','tprep96','tdisp96','tregn96','tlegn96','tdpegn96','trftn96'};
signature=struct('hvf_source',fileread('bin/source.sha256'), ...
    'hvf_binary',fileread('bin/executable.sha256'),'cps',cell(1,1));
signature.cps=cell(numel(programs),2);
for j=1:numel(programs)
    file=fullfile(cfg.cps_bin,programs{j});
    [status,message]=system([hash_command ' ' quote_shell_path(file)]);
    assert(status==0,'thbi:build','Cannot hash CPS program %s.',file);
    signature.cps(j,:)={programs{j},strtrim(message)};
end
end
