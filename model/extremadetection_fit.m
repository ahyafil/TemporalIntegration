function [T,sigma,lapse, S] = extremadetection_fit(X, Y, param)
%[T sigma] =  extremadetection_fit(X, Y) fits an extrema detection
%model to behavioral data, where X is the design matrix and Y are binary
%choices. A sensory gaussian noise assumes.
%
% Use nan values in sample evidence matrix X if the length of sequence
% varies across trials.
%
% extremadetection_fit(X, Y, param) for additional parameters. param is a
% structure with the following possible fields:
% - 'LastSample': boolean to indicate whether responses when neither
% boundary is reached at the end of the sample sequence are based on the
% last sample evidence (true) or randomly chosen (false, default)
% - 'Lapse': whether to include lapses into response (default: true)
% - 'FixedLapse': provide values of lapse rate parameters (vector of two)
% if these parameters are not fitted to data
% 'FixedThreshold' (default: true) whether there is a single fixed
% threshold across samples or varies across samples
% - 'nInit': number of initial points for maximum-likelihood optimization
% (default: 10)
% -'TolFun': optimality tolerance (default: 1e-9)
% - 'maxIter': maximum number of iterations (default: 1000)
%
% Output:
% T is the estimated mean threshold parameter, sigma is the estimated sensory noise.
%
% [T sigma lapse] =  extremadetection_fit(...)
% lapse are the estimated lapse rates for each response
%
% [T sigma lapse S] =  extremadetection_fit(...)
% S is a structure with various model metrics ('LLH','AIC','AICc','BIC','PosteriorCovariance'), standard
% error over parameters ('T_se','sigma_se','lapse_se') and probability of
% response for each trial ('pModel')
%
% See also extremadetection_rnd, extremadetection_test, snapshot_fit

%% 1. PROCESS ARGUMENTS

% default values
with_lapse = 1;
FixedLapse = [];
FixedThreshold = true;
nInit = 10;  % number of initial points
maxIter = 1000; % maximum number of iterations
TolFun = 1e-9; % function tolerance (when to stop optimization algorithm)
LastSample = false; % in case threshold is not reached during the stim, whether response is made based on last sample and not just chance

[nTrial, nSamp] = size(X); % number of observations and samples

% process parameters
if nargin<3
    param = struct;
end

fnames = fieldnames(param);
for f = 1:length(fnames)
    switch fnames{f}
        case 'Lapse'
            with_lapse = logical(param.Lapse);
        case 'FixedLapse'
            FixedLapse = param.FixedLapse;
            with_lapse = 1;
        case 'FixedThreshold'
            FixedThreshold = param.FixedThreshold;
        case 'LastSample'
            LastSample = param.LastSample;
        case 'nInit'
            nInit = param.nInit;
        case 'maxIter'
            maxIter = param.maxIter;
        case 'TolFun'
            TolFun = param.TolFun;
        otherwise
            error('incorrect field for param: %s', f);

    end
end


nLapse = 2*with_lapse; % number of lapse parameters
FreeLapse = with_lapse && isempty(FixedLapse); % whether there are free lapse parameters
nLapsePar = nLapse*FreeLapse; % number of free lapse parameters

if FixedThreshold
    nThresholdPar =1;
else
nThresholdPar = nSamp; %one threshold parameter per sample
end

if ~isvector(Y) && any((Y~=0) & (Y~=1))
    error('Y should be a vector of 0s and 1s');
end

%% 2. DEFINE INITIAL VALUES OF PARAMETERS



% initial values for threshold and sensory noise
par_ini(1,:) = [nanmean(abs(X(:)))*ones(1,nThresholdPar) nanstd(X(:))]; % first initialization: T set to mean absolute value, sigma to standard deviation
par_ini(2:nInit,:) = nanmax(X(:))*rand(nInit-1,nThresholdPar+1); % following initializations

lb = zeros(1,nThresholdPar+1); % lower bounds for parameters
ub = Inf(1,nThresholdPar+1);

if FreeLapse
    lapse_idx=nThresholdPar+(2:3);
    par_ini(1,lapse_idx) = [.01 .01]; % initial lapse value
    rnd_lapse  = drchrnd(ones(1,3),nInit-1); % draw from uniform Dirichlet distrib
    par_ini(2:nInit,lapse_idx) = rnd_lapse(:,1:2);
    lb(lapse_idx) = [0 0];
    ub(lapse_idx) = [1 1];
end

%% 3. FIND ML PARAMETERS

% Define objective function: neg-LLH function
nLLH_fun = @(x) negLogLikelihood(X, Y, x,FixedLapse, FixedThreshold, LastSample);

% options for optimization
options = optimoptions('fmincon', 'TolFun',TolFun, 'MaxIter',maxIter);

% find ML parameters
[par_hat,nLLH, exitflag, output, ~,~,Hessian] = fmincon_multinit(nLLH_fun, par_ini,[],[],[],[],lb,ub,[],options);

LLH = -nLLH; % log-likelihood

% recover estimated values of each parameter
T = par_hat(1:nThresholdPar); %estimated threshold
sigma = par_hat(nThresholdPar+1); % estimated sensory noise
if FreeLapse
    lapse = par_hat(lapse_idx); % initial lapse value
elseif ~isempty(FixedLapse)
    lapse = FixedLapse;
else
    lapse = [];
end

%% 4. COMPUTE LAPLACE APPROXIMATION AND MODEL METRICS

% use Laplace approximation to compute standard error of parameters
PostCov = inv(Hessian); % posterior covariance (Laplace approximation)
se = sqrt(diag(PostCov));
S.T_se = se(1:nThresholdPar);
S.sigma_se = se(nThresholdPar+1);

if FreeLapse
    S.lapse_se = se(lapse_idx);
else
    S.lapse_se = zeros(length(FixedLapse),1);
end

% likelihood of response 1 (to simulate model)
pResp = pModelResponse(X, par_hat,FixedLapse, FixedThreshold, LastSample);
pModel = pResp(:,2);


%% 5. MODEL METRICS
S.LLH = LLH;

nFreeParameters = 1+nThresholdPar + nLapsePar; % number of free parameters

S.nFreeParameters = nFreeParameters;
S.BIC = nFreeParameters*log(nTrial) -2*LLH; % Bayes Information Criterion
S.AIC = 2*nFreeParameters - 2*LLH; % Akaike Information Criterior
S.AICc = S.AIC + 2*nFreeParameters*(nFreeParameters+1)/(nTrial-nFreeParameters-1); % AIC corrected for sample size

% probability of response =1 for each trial, according to fitted model
S.pModel = pModel;

S.PosteriorCovariance = PostCov;

end


%% compute probability of response 1 for each trial, given parameters
function pResp = pModelResponse(X, pars,FixedLapse, FixedThreshold, LastSample)

[nTrial, nSamp] = size(X); % number of osbervations and samples

%process parameters
if FixedThreshold
    nThresholdPar=1;
    T = pars(1)*ones(1,nSamp);
else
nThresholdPar=nSamp;
T = pars(1:nThresholdPar); % threshold 
end

sigma = pars(nThresholdPar+1); % sensory noise
if length(pars)>nThresholdPar+1 % free lapse parameters
    lapse = pars(nThresholdPar+(2:3));
elseif ~isempty(FixedLapse) % fixed lapse parameters
    lapse = FixedLapse;
else % no lapses
    lapse = [];
end


% normalize evidence and threshold by sensory noise
X = X/sigma;
T = T/sigma;

% marginalized probability for response 0 (i.e. associated with
% negative threshold)
pResp0 = zeros(nTrial,1);

% probability that it has not responded yet
pNoResp = ones(nTrial,1);

% find number of samples per trial
nSamp_Trial = nSamp*ones(nTrial,1); % by default full length of sequence
for t=1:nTrial

    % find last non-nan values in sequence
    last_nonan = find(~isnan(X(t,:)),1,'last');
    if ~isempty(last_nonan)
        nSamp_Trial(t) = last_nonan;
    end
end

% process each sample sequentially
for s=1:nSamp

    % threshold vector
    TT = T(s)*ones(nTrial,1);
    if LastSample % last sample rule:
        TT(nSamp_Trial==s) = 0; % if last sample, threshold is collapsed to 0
    end
    NotFinished = nSamp_Trial>=s; % all trials for which sequence of samples is not extinguished yet

    % probability of each response at this sample if has not reached a response yet
    pRespIfNoResp = [normcdf(-X(NotFinished,s)-TT(NotFinished))   normcdf(X(NotFinished,s)-TT(NotFinished))];

    % add to marginalized probability for each response
    pResp0(NotFinished) = pResp0(NotFinished) + pNoResp(NotFinished).*pRespIfNoResp(:,1);

    % update probability of no response until this frame
    pNoResp(NotFinished) = pNoResp(NotFinished) .* (1-sum(pRespIfNoResp,2));

end

% assign probability of no response when sequence is finished equally
% between both options
if ~LastSample
    pResp0 = pResp0 + 0.5*pNoResp;
end

% add lapses
if ~isempty(lapse)
    pResp0 = lapse(1) + (1-sum(lapse))*pResp0;
end

% first column: prob of responses 0, second column: prob of responses 1
pResp = [pResp0 1-pResp0];
end

%% compute neg-loglikelihood of data given parameters
function nLLH = negLogLikelihood(X, Y, pars, FixedLapse, FixedThreshold, LastSample)

% compute probability for each response at each trial
pResp = pModelResponse(X, pars,FixedLapse, FixedThreshold, LastSample);

[nTrial, ~] = size(X); % number of osbervations and samples

% compute likelihood for observed response
pY = zeros(nTrial,1);
pY(Y==1) = pResp(Y==1,2);
pY(Y==0) = pResp(Y==0,1);

% neg-log-likelihood
nLLH = -sum(log(pY));
end