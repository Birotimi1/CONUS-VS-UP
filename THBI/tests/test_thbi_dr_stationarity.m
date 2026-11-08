function report = test_thbi_dr_stationarity(dr_mode)
% Exact finite-state variable-dimension test of the production acceptors.
% Two zones each have 1 or 2 binary-valued layers: 36 Earth-like states.
% Explicitly enumerate every first/second path and build the WHOLE transition
% matrix. This checks invariant mass, beyond reciprocal path algebra.
% This discrete test validates acceptance; continuous physical maps are
% separately checked by test_thbi_layered_dr and test_thbi_prior_counts.
par=bayes_inv_parms('extended60');
if nargin>0, par.proposal.dr.mode=dr_mode; end
par.proposal.dr.same_zone_weight=[.2 .8 .5]; % only used in targeted comparison mode
par.proposal.probability=[.6 0 0 .25 .15 0];
par.proposal.zone_probability=[.2 .3 .5];
saved=rng; cleanup=onCleanup(@() rng(saved)); %#ok<NASGU>
bits={0,1,[0;0],[0;1],[1;0],[1;1]}; N=36;
models=cell(N,2); states=cell(N,1); logtarget=zeros(N,1);
for a=1:6
    for b=1:6
        i=(a-1)*6+b; models(i,:)={bits{a},bits{b}};
        v=bits{a}; u=bits{b};
        logtarget(i)=-.7*numel(v)-.3*numel(u)+.8*sum(v)-.4*sum(u)-.6*(mean(v)-mean(u))^2;
        states{i}=struct('model',struct('K',[numel(v) numel(u) 1]), ...
            'logL',logtarget(i),'logP',0);
    end
end
pi=exp(logtarget-max(logtarget)); pi=pi/sum(pi);

%% =================== EXPLICIT BASE PROPOSAL MATRIX Q1 ===================
Q=zeros(N); move=zeros(N); zone=zeros(N);
for x=1:N
    for z=1:2
        v=models{x,z}; n=numel(v); pz=par.proposal.zone_probability(z);
        for j=1:n
            new=models(x,:); flipped=v; flipped(j)=1-flipped(j); new{z}=flipped;
            y=find_state(models,new);
            Q(x,y)=.6*pz/n; move(x,y)=1; zone(x,y)=z;
        end
        if n==1
            for bit=0:1
                new=models(x,:); new{z}=[v;bit]; y=find_state(models,new);
                probability=.7; if bit==1, probability=.3; end
                Q(x,y)=.25*pz*probability; move(x,y)=4; zone(x,y)=z;
            end
        else
            new=models(x,:); new{z}=v(1); y=find_state(models,new);
            Q(x,y)=.15*pz; move(x,y)=5; zone(x,y)=z;
        end
    end
end
% Unrepresented moves (third zone, impossible counts) are rejected at once.
assert(all(sum(Q,2)<1));
P=zeros(N); wrong_first_rule=zeros(N); wrong_targeting=zeros(N); paths=0;
for x=1:N
    for y=find(Q(x,:)>0)
        p1=proposal(x,y,Q,move,zone,states,par,NaN);
        [~,loga1]=b6_IFACCEPT(states{x},states{y},p1,false);
        a1=exp(min(0,loga1)); mass=Q(x,y)*a1;
        P(x,y)=P(x,y)+mass; wrong_first_rule(x,y)=wrong_first_rule(x,y)+mass;
        wrong_targeting(x,y)=wrong_targeting(x,y)+mass;
        if a1==1, continue; end
        for z=find(Q(y,:)>0)
            p2=proposal(y,z,Q,move,zone,states,par,p1.zone);
            second=states{z}; second.ptb=p2;
            [~,loga2,t]=delayed_rejection_layered(states{x},states{y},second,p1,par);
            q2=Q(y,z)*exp(p2.log_targeting_forward);
            pathmass=Q(x,y)*(1-a1)*q2;
            P(x,z)=P(x,z)+pathmass*exp(min(0,loga2));
            % Negative controls must visibly fail the SAME invariant target.
            % (1) Ordinary first-stage rule applied to the compound move,
            %     omitting both rejection factors, as in a simple retry.
            reverse2=proposal(y,x,Q,move,zone,states,par,p2.zone);
            naive=logtarget(z)-logtarget(x)+log(Q(y,x)/Q(x,y))+log(Q(z,y)/Q(y,z)) + ...
                reverse2.log_targeting_forward-p2.log_targeting_forward;
            wrong_first_rule(x,z)=wrong_first_rule(x,z)+pathmass*exp(min(0,naive));
            % (2) Full rejection factors but missing Q2 targeting correction.
            wrong_targeting(x,z)=wrong_targeting(x,z)+ ...
                pathmass*exp(min(0,loga2-t.log_targeting_ratio));
            paths=paths+1;
        end
    end
    P(x,x)=P(x,x)+1-sum(P(x,:));
    wrong_first_rule(x,x)=wrong_first_rule(x,x)+1-sum(wrong_first_rule(x,:));
    wrong_targeting(x,x)=wrong_targeting(x,x)+1-sum(wrong_targeting(x,:));
end
flux=diag(pi)*P;
report=struct('dr_mode',par.proposal.dr.mode,'states',N,'second_paths',paths, ...
    'stationarity_error',max(abs(pi'*P-pi')), ...
    'balance_error',max(abs(flux(:)-reshape(flux',[],1))), ...
    'naive_retry_error',max(abs(pi'*wrong_first_rule-pi')), ...
    'missing_targeting_error',max(abs(pi'*wrong_targeting-pi')));
assert(all(P(:)>=-1e-14) && max(abs(sum(P,2)-1))<1e-14);
assert(report.stationarity_error<1e-12 && report.balance_error<1e-12, ...
    'thbi:test','Full delayed-rejection transition matrix does not preserve the target.');
assert(report.naive_retry_error>1e-6,'thbi:test','Negative control did not detect omitted rejection factors.');
if strcmp(par.proposal.dr.mode,'targeted_layered')
    assert(report.missing_targeting_error>1e-6,'thbi:test','Targeting negative control failed.');
else
    assert(report.missing_targeting_error<1e-12,'thbi:test','Random kernel acquired a targeting correction.');
end
fprintf('PASS: exact %d-state %s RJ transition matrix, %d second paths.\n',N,par.proposal.dr.mode,paths);
disp(report);
end

function i=find_state(models,new)
i=0;
for j=1:size(models,1)
    if isequal(models(j,:),new), i=j; return; end
end
assert(i>0,'thbi:test','Unknown enumerated state.');
end

function p=proposal(x,y,Q,moves,zones,states,par,previous)
assert(Q(x,y)>0 && Q(y,x)>0);
z=zones(x,y);
if strcmp(par.proposal.dr.mode,'layered_path'), previous=NaN; end
p=struct('type',par.proposal.names{moves(x,y)},'valid',true,'step_scale',1, ...
    'dr_mode',par.proposal.dr.mode,'zone',z,'previous_zone',previous, ...
    'K_before',states{x}.model.K,'K_after',states{y}.model.K, ...
    'log_q_forward',log(Q(x,y)),'log_q_reverse',log(Q(y,x)), ...
    'log_abs_jacobian',0,'log_targeting_forward',0);
if isfinite(previous) && previous>0
    w=par.proposal.dr.same_zone_weight(previous);
    weights=(1-w)*par.proposal.zone_probability; weights(previous)=weights(previous)+w;
    p.log_targeting_forward=log(weights(z)/par.proposal.zone_probability(z));
end
end
