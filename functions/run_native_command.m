function [status,message] = run_native_command(command,cfg,label)
% Execute an unchanged native command with live progress and a hard time limit.
% Timeouts abort the job: they must never be converted into MCMC rejections.
if nargin<2, cfg=struct(); end
if nargin<3, label='native forward'; end
cfg=merge_defaults(cfg,struct('native_heartbeat_seconds',60, ...
    'native_timeout_seconds',7200,'verbose',false));
interval=cfg.native_heartbeat_seconds; limit=cfg.native_timeout_seconds;
assert(all(isfinite([interval limit])) && all([interval limit]>=0) && ...
    all(mod([interval limit],1)==0),'thbi:native_config','Use integer nonnegative native timing settings.');
root=fileparts(fileparts(mfilename('fullpath')));
script=fullfile(root,'scripts','native_watchdog.sh');
assert(exist(script,'file')==2,'thbi:dependency','Missing native process monitor.');
wrapper=sprintf('bash %s %d %d %s %s',quote_shell_path(script),interval,limit, ...
    quote_shell_path(label),quote_shell_path(command));
if cfg.verbose, fprintf('[native] START %s | %s\n',label,pwd); end
started=tic;
[status,message]=system(wrapper,'-echo');
if cfg.verbose, fprintf('[native] END %s | %.1fs | exit %d\n',label,toc(started),status); end
if status==124
    error('thbi:native_timeout','%s exceeded %ds in %s. Resume the last checkpoint after investigating; no model was rejected for time cost.',label,limit,pwd);
elseif status==130 || status==143
    error('thbi:native_interrupted','%s interrupted in %s.',label,pwd);
end
end
