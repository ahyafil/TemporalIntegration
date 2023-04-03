function Y = extremadetection_rnd(Stimulus,mu, sigma, UseLastSample, Lapse)
% extremadetection_rnd simulates the extrema detection model
%
% Y = extremadetection_rnd(X, mu, sigma)
% - X is the nTrial x nSample matrix of stimulus evidence from the dynamic sequential stimuli
% (rows: trials; columns: sample positions)
% - mu is the absolute value of decision boundaries (scalar, or a vector of length nSample if the threshold varies across the sequence)
% - sigma is the standard deviation of gaussian sensory noise (scalar)
% - Y is a binary vector of simulated choices (1 if hits the positive threshold,
% 0 if hits the negative threshold)
%
%  Y = extremadetection_rnd(X, mu, sigma, UseLastSample) where
% UseLastSample is a boolean (default: false). If true, in trials where the boundary hasn't been reached at
%  the end of the sequence, the choice is made based on the sign of the
%  last sample. If false, the choice is completely random.
%
% Y = extremadetection_rnd(X,mu, sigma, UseLastSample, Lapse) to add lapse
% terms. Lapse should be a vector of two values for probability of lapses
% associated with the negative and positive boundaries (i.e. left and right
% lapse probabilities)

if nargin<4
    UseLastSample = false;
end
if nargin<5
    Lapse = [];
else
    assert(isvector(Lapse) && length(Lapse)==2 && all(Lapse>=0) && sum(Lapse)<=1,...
        'Lapse should be a vector of two non-negative values summing no more than one')
end

% extract number of trials and samples
[nTrial, nSamp] = size(Stimulus);

assert(isscalar(mu)||length(mu)==nSamp, 'mu should be a scalar or a vector of length nSample');
assert(isscalar(sigma), 'sigma should be scalar');

if isscalar(mu)
    mu = mu*ones(1,nSamp);
end

% initialize response
Y = zeros(nTrial,1);

HasNotReachedBoundary = true(nTrial,1); % boolean: true if still hasn't reached the boundary

for s=1:nSamp % for each sample

    % instantaneous sensory evidence (stimulus corrupted by gaussian noise)
    Xnoisy = Stimulus(:,s) + sigma*randn(nTrial,1); 

    % trials where this reaches negative threshold
    LeftResp = HasNotReachedBoundary & Xnoisy<-mu(s); % mask for those trials
    Y(LeftResp) = 0; % set response
    HasNotReachedBoundary(LeftResp) = false; % has reached boundary

    % trials where this reaches the positive threshold
    RightResp = HasNotReachedBoundary & Xnoisy>mu(s); % reaches positive threshold
    Y(RightResp) = 1;
    HasNotReachedBoundary(RightResp) = false;
end

%% if threshold is not reached after all the sequence is finished
if UseLastSample % use last sample
    LeftResp = HasNotReachedBoundary & Xnoisy<0; % reaches negative threshold
    Y(LeftResp) = 0;
    HasNotReachedBoundary(LeftResp) = false;

    RightResp = HasNotReachedBoundary & Xnoisy>0; % reaches positive threshold
    Y(RightResp) = 1;
else % else, decide randomly

    RdmChoice = rand(nTrial,1)>.5;
    Y(HasNotReachedBoundary & RdmChoice) = 1;
    Y(HasNotReachedBoundary & ~RdmChoice) = 0;
end

%% add lapses
if ~isempty(Lapse)
    % override responses if used lapses
    ppLapse = rand(nTrial,1); % draw value from uniform [0 1]for each trial: lowest values are classified are lapses
    Y(ppLapse <= Lapse(1)) = 0; % left lapses
    Y(ppLapse > Lapse(1) & ppLapse<= sum(Lapse)) = 1; % right lapses
end
end