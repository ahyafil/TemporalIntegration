function [Y,attended] = snapshot_rnd(X,pi,beta)
%Y = snapshot_rnd(X) simulated the snapshot model in deterministic mode
%(i.e. with no sensory noise). 
% X is the nTrial-by-nSample matrix of sensory evidence, with nTrial the
% number of trials and nSample the number of samples in each stimulus. Pad
% with nans if sequence is shorter than nSample for some trials. 
% Y is the vector of simulated responses (of length nTrial).
% The model randomly picks a sample on each trial and outputs a response
% according to the sign of X for te corresponding sample (1 if positive, 0
% if negative).
% 
% snapshot_rnd(X,pi) allows for non-uniform draws for snapshot. pi is the
% vector of mixture coefficients, i.e. the vector of probabilities of length nSample 
% indicating the probability of using each sample (values must sum to one). 
% pi can also be a vector of length nSample+2 to include lapses: 
% pi(nSample+1) corresponds to the probability of a lapse with response 0, 
% pi(nSample+2) corresponds to the probability of a lapse with response 1.
%
% snapshot_rnd(X,pi, beta) to simulate the stochastic (non-deterministic) version of the
% snapshot model, i.e. given that the attended sample is s, response 1 is
% selected with probability 1/(1+exp(-beta(i)*X(t,s)). The sample evidence is
% weighted by the corresponding weight in beta and passed through a logistic function.
%
% [Y,attended] = snapshot_rnd(...)
% attended is the vector for the attended sample in each trial
%
% See also snapshot_fit, snapshot_test, extremadetection_rnd

[nTrial,nSample] = size(X); % number of trials and samples


if nargin<2
pi = ones(1,nSample)/nSample; % uniformly distributed between samples
else
    assert(all(pi>=0) && abs(sum(pi)-1)<1e-6,'pi must be a vector of non-negative values summing to one.')
end
assert(length(pi)==nSample || length(pi)==nSample+2, 'pi must be a vector of length nSample or nSample+2');

is_deterministic = nargin<3;
assert(is_deterministic || length(beta)==nSample,'beta must be a vector of length nSample');

% 1. Select which sample (or lapse) is attended on each trial
attended = mnrnd(1,pi,nTrial); % nTrial x nSample matrix of 0s and 1s

%  2. Select the response for each trial conditioned on the attended frame
%  (ignoring lapses)
if is_deterministic
     % simply the sign of the stimulus
    resp = (X>0);
else
    % pass through logistic function with sample-dependent sensitivity
    p_resp = 1./(1+exp(-beta.*X)); 

    % generate response according to probability
    resp = binornd(1,p_resp); 
end

% add responses for lapses
if length(pi)>nSample % if lapses are included
resp = [resp zeros(nTrial,1) ones(nTrial,1)]; % systematically 0 for one, systematically 1 for the other
end

% 3. Select response based on attended sample on each trial
Y = sum(attended.*resp,2); % response according to attended
end