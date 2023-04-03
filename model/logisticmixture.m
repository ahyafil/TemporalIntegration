function [w, mixt, S] = logisticmixture(X, Y, Ridge, nMixture, x0,param)
% function  [w, mixt, S] = logisticmixture(X, Y, Ridge,nMixture, x0 [,param])
%
% Mixture of Logistic regression models.
% w is a matrix nRegressor x nMixture
%
% add help!


%Design matrix A, targets Y,
% Use probitreg(a,y,Ridge) to use MAP approximation
% where prior for weights is multivariate Gaussian with precision sqrt(Ridge) :
% p(w) = sqrt(prod(Ridge)/(2*pi)^npar) * exp( -sum(Ridge.*w(i)^2)/2 );
% prior for lapse is flat Dirichlet:
% p(lapse) = Dir(lapse|alphas =1) = 1/2

% add help

if (nargin < 6)
    param = [];
end

[nObs, nReg,nD3] = size(X); % number of observations and regressors
samereg = ismatrix(X); % whether same set of regressors for each component or not (then, design matrix for each component)

if ~samereg && nD3 ~= nMixture
    error('if X is of dimension 3, the number of elements along dimension 3 must match the number of components');
end

if ~isvector(Y)
    error('Y must be a vector');
elseif length(Y) ~= nObs
    error('number of trials in X and Y do not match');
end
Y = Y(:);

% add constant vector (for bias)
with_bias = ~isfield(param, 'bias') || param.bias;
if with_bias
    X = cat(2,ones(nObs,1,size(X,3)), X);
    nReg = nReg + 1;
end

% whether we add lapses parameters
if ~isfield(param, 'lapse')
    param.lapse = true;
end

nMixtureTot = nMixture + 2*param.lapse; % lapses added two mixtures
nPar = nMixture*nReg + nMixtureTot; % number of parameters: each set of weight + for mixtures


% empty vector : return all with nan
if nObs ==0
    S.resid = nan(0,1);
    S.w = nan(1,nReg-with_bias);
    S.bias =  nan(1,1*with_bias);
    S.mixture = nan(1,nMixtureTot);
    S.exitflag = 0;
    S.niter = 0; % total number of iterations
    S.covb = nan(nPar);
    S.se = nan(1,nPar);
    S.T = nan(1,nPar);
    S.p = nan(1,nPar);
    S.w_se = nan(1, nReg-with_bias);
    S.bias_se = nan(1,1*with_bias);
    S.mixture_se = nan(1,nlapse);
    S.w_p = nan(1, nReg-with_bias);
    S.bias_p = nan(1,1*with_bias);
    S.mixture_p = nan(1,nlapse);
    S.joint = [[S.bias, S.w, S.lapse]; S.se]';
    S.LLH_ML = 0;
    S.LLH = 0;
    S.xent = nan; % cross-entropy estimate
    S.certainty = nan; % exponential of cross-entropy estimate (certainty of the model)
    S.x_all = nan(nPar, 1);
    S.LLH_all = [];
    S.exitflag_all = [];
    return;
end

% whether use fixed weights (only fitting mixing coeffs)
if ~isfield(param, 'fixedweights')
    param.fixedweights = [];
end

% by default : no prior
if nargin<3
    Ridge =  zeros(nReg,1);
end
if isfield(param, 'alpha') % Dirichlet prior
    alpha = param.alpha;
else
    alpha = ones(1,nMixtureTot) ; % uninformative Dirichlet prior
end
alpha = alpha(:);
logB = sum(gammaln(alpha)) - gammaln(sum(alpha)); % log of normalizing factor of Dirichlet (multivariate beta function)

if (length(Ridge) == 1)
    Ridge = Ridge*ones(1,nReg);
elseif length(Ridge(:)) == nReg
    Ridge = Ridge(:)';
else
    error('Ridge weight vector should be length 1 or %d', nReg);
end
if any(Ridge<0)
    error('Ridge values must be positive or null');
end

% number of initial points
if nMixture ==1
    nInit = 1; % with just one component, LLH is convex so always converges to global maxima
elseif isfield(param,'ninit')
    nInit = param.ninit;
else
    if isempty(param.fixedweights)
        nInit = 10;
    else % fixed set of weights
        nInit = 1;
    end
end

% initial value for parameter
if ((nargin < 5) || (isempty(x0)))
    if samereg % if the regressors are the same, we need to assign different weights at the beginning
        amax = max(abs(max(X,[],3)),[],1); % maximum value for each factor
        maxw = 8e-3./ amax / nReg; % maximum initial values for weight: make sure that no probit reaches 0 or 1
        maxw(amax==0) = 0; % in case there is null factor
        w = bsxfun(@times, maxw', rand(nReg,nMixture)); % draw random weights
    else
        w = zeros(nReg, nMixture);
    end
    mixt = alpha/sum(alpha); % maximum for prior
    x0 = [w(:);mixt];

else % provided
    if isvector(x0)
        x0 = x0(:);
    end
end

if nInit>1 && size(x0,2)<2
    amax = max(abs(max(X,[],3)),[],1); % maximum value for each factor
    maxw = 8./ amax / nReg; % maximum initial values for weight: make sure that no probit reaches 0 or 1
    maxw(amax==0) = 0; % in case there is null factor
    w = bsxfun(@times, maxw', rand(nReg,nMixture,nInit-1)); % draw random weights
    w = reshape(w,nReg*nMixture, nInit-1);

    if nMixtureTot>1
        mixt = gamrnd(repmat(alpha',nInit-1,1),1,nInit-1,nMixtureTot); % initial lapse value: sample from prior Dirichlet distribution
        mixt = mixt ./ repmat(sum(mixt,2),1,nMixtureTot); % normalize and exclude last term (proba for nolapse)
        %   lapse = .5*rand(nlapse, ninit-1); % all lapse values between 0 and .5
    else
        mixt = ones( nInit-1,1);
    end
    x0(:,2:nInit) = [w;mixt'];
end


% exclude data with nan
exclVar = any(any(isnan(X),1),3);
exclPar = [repmat(exclVar,1,nMixture) false(1,nMixtureTot)];

X(:,exclVar,:) = [];
Ridge(exclVar) = [];
x0(exclPar,:) = [];
nReg_withnan = nReg;
nPar_withnan = nPar;
nReg = nReg - sum( exclVar);
nPar = nPar - sum(exclVar);

illcond = (nPar>nObs);
if illcond
    warnstat1 = warning('OFF', 'MATLAB:nearlySingularMatrix');
    warnstat2 = warning('OFF', 'MATLAB:singularMatrix');

    warning(message('stats:glmfit:IllConditioned'));
end


% log of Ridge scaling factor (ensures that prior sums to 1)
withRidge = Ridge>0;
k_ridge = 1;
%k_ridge = log ( (prod(Ridge(withRidge))/(2*pi)^sum(withRidge)) )/2;

% Ridge diagonal matrix
Ridgemat =  spdiags(Ridge',0,nReg,nReg);
%Ridge_w =  Ridge(1:nvar);
%Ridgemat_w = Ridgemat(1:nvar,1:nvar);

if ~isfield(param, 'maxiter')
    param.maxiter = 1000;
end

if ~isfield(param, 'verbose')
    param.verbose = 0;
end

% tolerance criterion for convergence (on LLH)
if ~isfield(param, 'epsilon')
    param.epsilon =  1e-3;
end

if ~isfield(param, 'maxprint')
    param.maxprint = 5;
end

crossvalid = 0;
if isfield(param, 'crossvalidation') && ~isempty(param.crossvalidation)
    crossvalid = 1;
    if iscell(param.crossvalidation) % permutationas are already provided as nperm x 2 cell (first col: train set, 2nd: train set)
        generateperm = 0;
        allperm =   param.crossvalidation; % all permutations
        nperm = size(allperm,1); % number of permutations
        ntrain = length(allperm{1,1});
        ntest = length(allperm{1,2});

    else % we wil draw the permutations
        ntrain = param.crossvalidation; % number of trials in training set
        generateperm = 1;

        if ntrain == -1
            ntrain = nObs-1; % leave one out
        elseif ntrain<1 % defined as percentage
            ntrain = round(ntrain*n);
        end
        if isfield(param, 'ntest') % defined number of test trials
            ntest = param.ntest;
        else % by default all trials not in train set
            ntest = nObs-ntrain;
        end
        if isfield(param, 'nperm')
            nperm = param.nperm;
        else
            nperm = 100; % default number of permutations
        end
    end
end

yy = [1-Y(:) Y(:) ]; % first column for response 0, second for response 1 (allows also to enter y not as binary but as responsability (real value between 0 and 1, for EM algo)
yy_s = [Y(:)-1 Y(:) ]; % signed version

i_w = 1:nReg*nMixture; % index for variables, excluding bias
i_mixt = nReg*nMixture + (1:nMixtureTot); % index for mixture variables


%%  training on all dataset
[w,mixt,exitflag,logpost,LLH_MAP,S]  = fitting_multinit(X, nInit, x0,i_w,i_mixt,...
    Y,yy, k_ridge,Ridge, Ridgemat, alpha, nMixture, logB, param);
if  exitflag<=0
    warning('multilogistic:notconverged', 'Failed to converge');
end

%% cross validation
if crossvalid

    validationscore = zeros(1,nperm); % score for each (LLH per trial)
    accuracy = zeros(1,nperm); % proportion correct
    exitflag_CV = zeros(1,nperm); % exit flag (converged or not)
    w_CV = zeros(nReg,nMixture, nperm);
    mixt_CV = zeros(nMixture,nperm);

    parfor p=1:nperm % for each permutation
        if generateperm % generate permutation
            this_ntrain = ntrain;
            trainset = randperm(n,this_ntrain); % training set
            notrainset= setdiff(1:n, trainset); % trials not in training set
            validationset = notrainset(randperm(n-ntrain,ntest)); % test set
        else
            traintest = allperm(p,:);
            trainset = traintest{1}; %allperm{p,1};
            validationset = traintest{2};%allperm{p,2};
            this_ntrain = length(trainset);
        end

        %fit on training set
        [this_w,this_mixt,exitflag_CV(p)] = fitting_multinit(X(trainset,:,:), nInit, x0, i_w,i_mixt,...
            Y(trainset),yy(trainset,:),  k_ridge,Ridge, Ridgemat, alpha, nMixture, logB, param);

        %compute score on testing set (mean log-likelihood per trial)
        [validationscore(p),accuracy(p)] = loglike(X(validationset,:,:), yy(validationset,:),this_w, this_mixt, param.lapse);
        validationscore(p) = validationscore(p)/ length(validationset); % normalize by number of trials

        w_CV(:,:,p) = this_w;
        mixt_CV(:,p) = this_mixt;
    end


    n_nonconverged = sum(exitflag_CV<=0);
    if  n_nonconverged>0
        warning('multilogistic:notconverged', 'Failed to converge for %d/%d permutations', n_nonconverged, nperm);
    end

    S.w_CV = w_CV;
    S.mixt_CV = mixt_CV;

    S.validationscore = mean(validationscore);
    S.validationscore_all = validationscore;
    S.accuracy_validation = mean(accuracy);
    S.accuracy_all = accuracy;
    S.exitflag_CV = exitflag_CV;
    S.converged_CV = sum(exitflag_CV>0); % number of permutations with convergence achieved
end

% residuals
A = activation(X,w);
Pm = 1 ./ (1 + exp(-A)); % logistic function (prob of response 1 according to model)
if param.lapse % add predictions of lapse-0 and lapse-1 model
    Pm = [Pm zeros(nObs,1) ones(nObs,1)];
end
p1 = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
lhh = [1-p1 p1]; % probability of each of the responses

zz1 = (Pm .* mixt')./p1; % responsability, assuming response 1
zz0 = ((1-Pm).*mixt')./(1-p1); % responsability, assuming response 0
S.responsibility = yy(:,1).*zz0 + yy(:,2).*zz1; % responsability % posterior probability for lapse given current parameters
S.Y = p1;
S.resid = Y-p1;


% weights and mixtures

S.w = nan(nMixture,nReg_withnan);
S.w(:,~exclVar) = w';

S.mixture = mixt';

S.n = nObs;
S.exitflag = exitflag;
S.niter = sum(S.all_niter); % total number of iterations

%% Laplace approximation
if exitflag>0  % if converged, compute full Hessian (observed information matrix)
    H = zeros(nPar); % pre-allocate Hessian matrix

    % compute Hessian for weight terms within same model

    z_nl = yy_s ./ lhh;
    z_nl(yy==0) = 0; % to avoid nan values for 0/0
    srat1 = sum(z_nl,2);

    rat2 = yy ./ lhh.^2;
    rat2(yy==0) = 0; % avoid nan values
    srat2 = sum(rat2,2);

    lh = sum(yy .* lhh,2); %likelihood of observed response

    grad = zeros(nObs,nPar); % gradient of p(right) for each trial
    grad(:,nReg*nMixture+(1:nMixtureTot)) = Pm; % gradient w.r.t mixture coefficients

    yy_ss = sum(yy_s,2); % -1 if 0 resp, +1 if 1 resp

    for m=1:nMixture
        deriv = Pm(:,m) .* (1-Pm(:,m));

        this_iw = m*nReg+(1-nReg:0); % index of wieght parameters
        this_im = nMixture*nReg + m; % index of mixture parameter

        if samereg
            grad(:,this_iw) = mixt(m)*deriv.*yy_ss.*X; % gradient w.r.t weights
        else
            grad(:,this_iw) = mixt(m)*deriv.*yy_ss.*X(:,:,m); % gradient w.r.t weights
        end

        % compute Hessian for weight terms
        hh = yy_ss.*deriv.*(1-2*Pm(:,m))./lh;
        if samereg
            H_ww =  mixt(m) * X' *  spdiags(hh, 0, nObs, nObs) * X; % Hessian of LLH
        else
            H_ww =  mixt(m) * X(:,:,m)' *  spdiags(hh, 0, nObs, nObs) * X(:,:,m); % Hessian of LLH
        end

        % add Ridge term (Bishop - Pattern recognition and Machine learning - eq 4.143)
        H_ww = H_ww - k_ridge * Ridgemat;

        if samereg
            H_wm = (srat1.*deriv)'*X; % hessian w.r.t to weight and mixture coeff
        else
            H_wm = (srat1.*deriv)'*X(:,:,m); % hessian w.r.t to weight and mixture coeff
        end


        H_mm = - (alpha(m)-1)/mixt(m)^2; % hessian due to Dirichlet prior

        % fill in big hessian matrix
        H(this_iw, this_iw) = H_ww;
        H(this_iw, this_im) = H_wm';
        H(this_im, this_iw) = H_wm;
        H(this_im, this_im) = H_mm;

    end

    % lapse parameters
    for m=nMixture+1:nMixtureTot
        this_im = nMixture*nReg + m; % index of mixture parameter
        H(this_im, this_im) = - (alpha(m)-1)/mixt(m)^2; % hessian due to Dirichlet prior
    end


    H = H - grad'*spdiags(srat2,0,nObs,nObs)*grad; % part of hessian due to double derivation

    % project to free basis (removing constraint that sum of mixtures is 1)
    P = zeros(nMixtureTot);
    for i=1:nMixtureTot-1
        P(i,:) = [ones(1,i) -i zeros(1,nMixtureTot-i-1)]/sqrt(i*(i+1)); % coordinates for i-th basis vector of free space
    end
    P = blkdiag(eye(nReg*nMixture), P); % add weight parameters

    Hfree = P*H*P'; % Hessian in free basis
    Hfree = (Hfree+Hfree')/2;

    %Covariance matrix
    S.covb = nan(nPar_withnan);
    invH = inv(-Hfree);

    if any(diag(invH)<0)
        warning('information matrix is not definite positive... tolerance criterion might be set too high');
    end
    if ~isempty(param.fixedweights) % no covariance for fixed weight
        invH(1:nReg*nMixture,:) = 0;
        invH(:,1:nReg*nMixture) = 0;
    end
    S.covb(~exclPar,~exclPar) = P'*invH* P; % see covariance under constraint Seber & Wild Appendix E


    % standard error of estimates
    S.se = sqrt(diag(S.covb))';

    % T-statistic for the weights
    S.T = nan(1,nPar_withnan);
    S.T(~exclPar) = [w(:)' mixt'] ./ S.se(~exclPar);

    % p-value for significance of each coefficient
    S.p = 1-chi2cdf(S.T.^2,1)';

else % did not converge: cannot estimate observed information matrix

    warning('probitreg:notconverged', 'Failed to converge');
    H = nan(nPar);

    S.covb = nan(nPar_withnan);
    S.se = nan(1,nPar_withnan);
    S.T = nan(1,nPar_withnan);
    S.p = nan(1,nPar_withnan);
end

i_w = 1:nReg_withnan*nMixture; % index for variables, excluding bias
i_mixt = nMixture*nReg_withnan+(1:nMixtureTot); % index for mixture variables

S.w_se = reshape(S.se(i_w),nReg_withnan,nMixture);
S.mixt_se = S.se(i_mixt);
S.w_p = reshape(S.p(i_w),nReg_withnan,nMixture);
S.mixt_p = S.p(i_mixt);
S.joint = [[ S.w(:)', S.mixture]; S.se]';

% LLH at Max Likelihood / Max A Posteriori parameters
S.LLH_MAP = LLH_MAP;

% model evidence using Laplace approximation (Bishop - eq 4.137)  -
% requires that a prior has been defined
if exitflag>0 && all(Ridge>0)
    S.LLH = logpost - Ridge*w.^2/2 + nPar/2*log(2*pi) - log(det(-H))/2;
    S.LLH = S.LLH-nMixtureTot; % Dirichlet prior (correct formula!!)
    S.LLH = nan;
else
    S.LLH = nan;
end

S.nObs = nObs;
S.BIC = nPar*log(S.nObs) -2*LLH_MAP; % Bayes Information Criterion
S.AIC = 2*nPar - 2*LLH_MAP; % Akaike Information Criterior
S.AICc = S.AIC + 2*nPar*(nPar+1)/(nObs-nPar-1); % AIC corrected for sample size

S.xent = S.LLH/nObs; % cross-entropy estimate
S.certainty = exp(S.xent); % exponential of cross-entropy estimate (certainty of the model)

%% test set
if isfield(param, 'testset')
    [S.testscore, S.accuracy_test] = loglike(X(param.testset,:,:), yy(param.testset,:),w, mixt, param.lapse) ;
    S.testscore = S.testscore/ length(param.testset); % normalize by number of trials
end

S.n_init = nInit; % number of initial points
S.x_all = zeros(nPar_withnan, nInit);
S.x_all(~exclPar,:) = S.all_x;

if illcond
    warning(warnstat1.state, 'MATLAB:nearlySingularMatrix');
    warning(warnstat2.state, 'MATLAB:singularMatrix');
end

end


%% compute predicting accuracy and LLH for given parameter set
function [testscore,accuracy] = loglike(X, yy,w, mixt,lapse)

A = activation(X,w); % activation: multiply design matrix by weight matrix
Pm = 1 ./ (1 + exp(-A)); % logistic function (prob of response 1 without lapse)
if lapse % add predictions of lapse-0 and lapse-1 model
    nObs = size(Pm,1);
    Pm = [Pm zeros(nObs,1) ones(nObs,1)];
end


Y = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
lhh = [1-Y Y]; % probability of each of the responses

%% log-likelihood of observed data
ylh = yy .* log(lhh);
ylh(yy==0) = 0;
testscore = sum(ylh(:));

lh_o = sum(yy .*lhh,2); % proba of observed responses
accuracy = mean(lh_o>.5); % accuracy
end


%% fitting algo (multiple starting points)
function [w,mixt,exitflag,logpost,LLH_MAP,S]  = fitting_multinit(X, nInit, x0, i_w,i_mixt,...
    Y,yy, k_ridge,Ridge, Ridgemat, alpha, nmixture, logB, param)

%initialize variables
nPar = size(x0,1);
nReg = length(i_w)/nmixture;
S.all_x = zeros(nPar,nInit);
S.all_logpost = zeros(1,nInit);
S.all_LLH_MAP = zeros(1,nInit);
S.all_LLH_iter = cell(1,nInit);
S.all_exitflag = zeros(1,nInit);
S.all_niter =  zeros(1,nInit);

for n = 1:nInit

    w = reshape(x0(i_w,n),nReg, nmixture); % set of weights for each mixture
    mixt = x0(i_mixt,n); % mixture parameters


    % find best-fitting parameters from starting point
    [w,mixt,converged,iter,logpost,LLH,logpost_iter] = fitting(X, w, mixt,...
        Y,yy, k_ridge,Ridge, Ridgemat, alpha, nmixture, logB, param);

    % put that in big structure
    S.all_x(:,n) = [w(:);mixt];
    S.all_exitflag(n) = converged;
    S.all_niter(n) = iter;
    S.all_logpost(n) = logpost;
    S.all_LLH_MAP(n) = LLH;
    S.all_LLH_iter{n} = logpost_iter;
end

% find optimization procedure that yielded the best overall LLH
[logpost,imax] = max(S.all_logpost);
LLH_MAP = S.all_LLH_MAP(imax);
x = S.all_x(:,imax);
w = reshape(x(i_w),nReg,nmixture); % weights matrix
mixt = x(i_mixt); % lapse
exitflag = S.all_exitflag(imax);
S.logpost_iter = S.all_LLH_iter{imax};

end

%% fitting algorithm (single starting point)
function  [w,mixt,converged,iter,logpost,LLH,logpost_iter] = fitting(X, w, mixt,Y, yy, k_ridge, Ridge, Ridgemat, alpha, nMixture, logB, param)

SameReg = ismatrix(X); % whether all components have same regressors or component-specific regressors
nObs =size(X,1);
nMixtureTot = nMixture + 2*param.lapse;

%%deal with fixed weights first (simply one pass)
if ~isempty(param.fixedweights)
    w = param.fixedweights;
    A = activation(X,w);
    Pm = 1 ./ (1 + exp(-A)); % logistic function (prob of response 1 without lapse)
    if param.lapse % add predictions of lapse-0 and lapse-1 model
        Pm = [Pm zeros(nObs,1) ones(nObs,1)];
    end

    % compute likelihood of each trial for given parameters
    lh = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
    zz1 = (Pm .* mixt')./lh; % responsability, assuming response 1
    zz0 = ((1-Pm).*mixt')./(1-lh); % responsability, assuming response 0
    z_l = yy(:,1).*zz0 + yy(:,2).*zz1; % responsability

    mixt = (sum(z_l,1)' + (alpha-1))/(sum(z_l(:))+sum(alpha)-nMixtureTot);
    mixt = min(mixt,1); % keep lapses between 0 and 1 ( only appears if alpha<0)
    mixt = max(mixt,0);

    lh = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
    lhh = [1-lh lh]; % probability of each of the responses

    %% log-likelihood of observed data
    ylh = yy .* log(lhh);
    ylh(yy==0) = 0;
    LLH = sum(ylh(:));
    logpost = LLH - k_ridge * Ridge*sum(w.^2,2)/2;         % weight prior
    logpost = logpost + sum( (alpha-1).*log(mixt)) - logB;% Dirichlet prior

    converged = 1;
    iter = 1;
    logpost_iter = logpost;
    return;
end



A = activation(X,w);
Pm = 1 ./ (1 + exp(-A)); % logistic function (prob of response 1 without lapse)
if param.lapse % add predictions of lapse-0 and lapse-1 model
    Pm = [Pm zeros(nObs,1) ones(nObs,1)];
end


%oldproby = -ones(size(y));
oldLLH = -Inf;
logpost_iter = zeros(1,param.maxiter);
still =1;
iter = 0;

%% ITERATION PROCEDURE
while still

    for m=1:nMixture % for each model

        %% E-step 1 : estimate responsability of each model for each trial

        % compute likelihood of each trial for given parameters
        lh = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
        zz1 = (Pm .* mixt')./lh; % responsability, assuming response 1
        zz0 = ((1-Pm).*mixt')./(1-lh); % responsability, assuming response 0
        z_l = yy(:,1).*zz0 + yy(:,2).*zz1; % responsability

        % compute likelihood of each trial for given parameters
        % lh = expy * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
        % lhh = [1-lh lh ]; % probability of each of the responses (first column: resp 0, second column: resp1)

        %% M-step 1: update weights only
        if SameReg
            grad_w = X' * (z_l(:,m) .* (Pm(:,m)-Y)); % compute gradient over weights
        else
            grad_w = X(:,:,m)' * (z_l(:,m) .* (Pm(:,m)-Y)); % compute gradient over weights
        end
        grad_w = grad_w +  k_ridge * Ridge' .* w(:,m); % add Ridge term

        RR = z_l(:,m).* Pm(:,m) .* (1-Pm(:,m));
        if SameReg
            H_ww = X' * spdiags(RR,0, nObs, nObs) *X; % Hessian of Q
        else
            H_ww = X(:,:,m)' * spdiags(RR,0, nObs, nObs) *X(:,:,m); % Hessian of Q

        end


        % add Ridge term (Bishop - Pattern recognition and Machine learning - eq 4.143)
        H_ww = H_ww + k_ridge * Ridgemat;

        %Newton-Raphson update on weights
        w(:,m) = w(:,m) - H_ww\grad_w;

        %% E-step again 2: estimate posterior probability for each mixture for each trial

        A = activation(X,w);
        Pm = 1 ./ (1 + exp(-A)); % logistic function (prob of response 1 without lapse)
        if param.lapse % add predictions of lapse-0 and lapse-1 model
            Pm = [Pm zeros(nObs,1) ones(nObs,1)];
        end

        if nMixtureTot >1


            % compute likelihood of each trial for given parameters
            lh = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
            zz1 = (Pm .* mixt')./lh; % responsability, assuming response 1
            zz0 = ((1-Pm).*mixt')./(1-lh); % responsability, assuming response 0
            z_l = yy(:,1).*zz0 + yy(:,2).*zz1; % responsability

            %% M-step 2: update mixture parameters only (and compute LLH for new parameters)

            mixt = (sum(z_l,1)' + (alpha-1))/(sum(z_l(:))+sum(alpha)-nMixture);
            if mixt<0
                warning('ede');
            end
            mixt = min(mixt,1); % keep lapses between 0 and 1 ( only appears if alpha<0)
            mixt = max(mixt,0);

        end
    end

    % compute LLH for new parameters
    lh = Pm * mixt; % likelihood of response 1: sum up likelihood from each model, weighted
    lhh = [1-lh lh]; % probability of each of the responses

    %% log-likelihood of observed data
    ylh = yy .* log(lhh);
    ylh(yy==0) = 0;
    LLH = sum(ylh(:));
    logpost = LLH - k_ridge * Ridge*sum(w.^2,2)/2;         % weight prior
    logpost = logpost + sum( (alpha-1).*log(mixt)) - logB;% Dirichlet prior

    % check whether should iterate more
    iter = iter + 1; % update iteration counter;

    % has converged if improvement in LLH is smaller than epsilon
    converged = abs(oldLLH-logpost)<param.epsilon;
    oldLLH = logpost;
    logpost_iter(iter) = logpost;

    still = (iter<param.maxiter) && ~converged && ~any(isnan(w(:))) && ~any(isnan(mixt));

end

% sort models by decreasing mixing coefficients
if SameReg
    [mixt(1:nMixture), ord] = sort(mixt(1:nMixture), 'descend');
    w = w(:,ord);
end

logpost_iter(iter+1:end) = [];

end

%% compute activation from design matrix and weights
function A = activation(X,w)
if ismatrix(X) % same regressors for each component
    A = X * w; % multiply design matrix by weights
else % component-specific regressors
    A = sum(X .* permute(w,[3 1 2]),2); % multiply component-specific weight and sum over regressors
    A = permute(A,[1 3 2]);
end
end

