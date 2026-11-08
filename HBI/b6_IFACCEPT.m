function accept = b6_IFACCEPT(logPost, trial_logPost)
% Log-space Metropolis-Hastings for symmetric reflected proposals, T=1.
%% Reject invalid trials before drawing an acceptance variate
if ~isfinite(trial_logPost)
    accept = false;
    return;
end

%% Compare a uniform draw with the posterior ratio, both in log space
logRatio = trial_logPost-logPost;
logUniform = log(rand);
accept = logUniform < min(0,logRatio);
end
