function test_thbi_failures(profile)
% Error-propagation and scratch preservation checks with a cheap forward map.
% No native executables are invoked.
if nargin<1, profile='pilot'; end
par=bayes_inv_parms(profile); par.mod.start.jitter=false;
par.inv.niter=8; par.inv.burnin=2; par.inv.max_consecutive_failures=2;
par.proposal.probability=[0 0 1 0 0 0]; % boundary always gives a structural trial
root=tempname; mkdir(root); cleanup=onCleanup(@() rmdir(root,'s')); %#ok<NASGU>
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d};
    data.(name)=struct('x',(1:3)','obs',[1;2;3],'sig',ones(3,1));
end
data=a2_PREPARE_DATA(data,par);

%% =================== FAILED INITIAL MODEL IS PRESERVED ===================
par.res.dir=fullfile(root,'initial_failure');
try
    run_one_chain(par,data,1,[],@initial_failure);
    error('thbi:test','Initial error was hidden.');
catch problem
    assert(strcmp(problem.identifier,'thbi:injected_initial_error'));
end
a=load(fullfile(par.res.dir,'chain_01','initial_forward_failure','failure.mat'));
assert(a.failure.context.iteration==0 && isfield(a.model,'matrix'));

%% =================== EXPECTED NATIVE FAILURE STOPS AT THE FAILURE LIMIT ===================
par.res.dir=fullfile(root,'native_failure');
try
    run_one_chain(par,data,1,[],@native_failure_after_start);
    error('thbi:test','Expected failure limit was not reached.');
catch problem
    assert(strcmp(problem.identifier,'thbi:forward_failure'));
end
chain_dir=fullfile(par.res.dir,'chain_01');
a=load(fullfile(chain_dir,'first_forward_failure','failure.mat'));
b=load(fullfile(chain_dir,'fatal_forward_failure','failure.mat'));
assert(a.failure.context.iteration==1 && b.failure.context.iteration==2);
assert(exist(fullfile(chain_dir,'first_forward_failure','native_files','started.txt'),'file')==2);

%% =================== PROGRAMMING ERROR PROPAGATES ON THE FIRST PROPOSAL ===================
par.res.dir=fullfile(root,'programming_error');
try
    run_one_chain(par,data,1,[],@programming_failure_after_start);
    error('thbi:test','Programming error was hidden as a rejected model.');
catch problem
    assert(strcmp(problem.identifier,'thbi:injected_programming_error'));
end
a=load(fullfile(par.res.dir,'chain_01','fatal_forward_failure','failure.mat'));
assert(a.failure.context.iteration==1);
assert(strcmp(a.failure.identifier,'thbi:injected_programming_error'));
%% =================== TIMEOUT ABORTS WITHOUT REJECTING A COSTLY EARTH ===================
par.res.dir=fullfile(root,'timeout');
try
    run_one_chain(par,data,1,[],@timeout_after_start);
    error('thbi:test','Timeout was hidden as a rejected model.');
catch problem
    assert(strcmp(problem.identifier,'thbi:native_timeout'));
end
saved=load(fullfile(par.res.dir,'chain_01','checkpoint.mat'),'checkpoint');
assert(saved.checkpoint.iter==0 && saved.checkpoint.chain.failures==0);
%% =================== SECOND-STAGE NATIVE FAILURE IS PRESERVED, NOT ACCEPTED ===================
if par.proposal.dr.enabled
    par.res.dir=fullfile(root,'second_stage_failure'); par.inv.max_consecutive_failures=1;
    try
        run_one_chain(par,data,1,[],@second_stage_failure);
        error('thbi:test','Second-stage failure was hidden.');
    catch problem
        assert(strcmp(problem.identifier,'thbi:forward_failure'));
    end
    a=load(fullfile(par.res.dir,'chain_01','fatal_forward_failure','failure.mat'));
    assert(a.failure.context.iteration==1 && ...
        strcmp(a.failure.context.stage,'delayed-rejection proposal'));
end
fprintf('PASS: programming errors propagate; failed Earth, context and scratch files survive cleanup.\n');
end

function predicted=initial_failure(varargin) %#ok<INUSD,STOUT>
error('thbi:injected_initial_error','Deliberately injected initial failure.');
end

function predicted=native_failure_after_start(model,data,par,workdir)
predicted=controlled_forward(model,data,par,workdir,'thbi:forward_failure');
end

function predicted=programming_failure_after_start(model,data,par,workdir)
predicted=controlled_forward(model,data,par,workdir,'thbi:injected_programming_error');
end

function predicted=controlled_forward(~,data,par,workdir,identifier)
marker=fullfile(workdir,'started.txt');
if exist(marker,'file'), error(identifier,'Deliberately injected test failure.'); end
fid=fopen(marker,'w'); assert(fid>=0); fprintf(fid,'Initial forward completed.\n'); fclose(fid);
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; predicted.(name)=zeros(data.(name).N,1);
end
end

function predicted=second_stage_failure(~,data,par,workdir)
marker=fullfile(workdir,'call_number.txt'); n=0;
if exist(marker,'file'), n=str2double(fileread(marker)); end
n=n+1; fid=fopen(marker,'w'); assert(fid>=0); fprintf(fid,'%d',n); fclose(fid);
if n>=3, error('thbi:forward_failure','Deliberate second-stage solver failure.'); end
for d=1:numel(par.inv.datatypes)
    name=par.inv.datatypes{d}; predicted.(name)=1e6*(n==2)*ones(data.(name).N,1);
end
end

function predicted=timeout_after_start(model,data,par,workdir)
predicted=controlled_forward(model,data,par,workdir,'thbi:native_timeout');
end
