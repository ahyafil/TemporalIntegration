% function S = probitlapse(X, Y)
%
% Probit regression with lapses.
% X is the predictor/design matrix, i.e. a matrix with columns for
% factors and rows for observations
% Y is the target vector, containing binary values (0/1).
%
% S = probitlapse(X, Y, Ridge) uses ridge regression for weight estimation.
% Ridge is either a scalar or a vector defining the
% value of ridge regression for each weight. The corresponding prior is multivariate Gaussian with precision sqrt(Ridge) :
% p(w) = sqrt(prod(Ridge)/(2*pi)^npar) * exp( -sum(Ridge.*w(i)^2)/2 );
%
% S = probitlapse(X, Y, Ridge,nlapse) defines in nlapse the number of lapse
% parameters
% in the model:
% - 0: no lapse
% - 1: symmetrical lapse (same value for both responses)
% - 2: asymmetrical lapses (one parameter for each response)
% [default]
%
% S = probitlapse(a, y, Ridge,nlapse, x0) defines initial point(s) for the
% minimization procedure. If not provided, a random starting point for the
% minimization procedure will be drawn for each iteration.
% x0 is either a vector defining the starting value
% for each parameter, or a matrix defining starting values for each
% iteration in separate rows.
%
% S = probitlapse(X, Y, Ridge,nlapse, x0,P) allows to define in structure P properties
% of the fitting procedure. Possible fields in P include:
% - 'bias': value 0 or 1 (default). Whether a constant bias is including in
% the design matrix (the corresponding weight is the first).
% - 'alpha': a (nlapse+1) vector defining the parameters of Dirichlet prior for lapses.
% By default, alphas are selected to yield 3% lapse rate for each response.
% - 'ninit': number of different initializations of the minimization
% procedure (default:1 if model without lapse, as problem is convex; 100 otherwise)
% - 'maxiter': maximum number of iterations for each minimization procedure
% (default:10000)
% - 'epsilon': tolerance value for log-posterior to stop minimization
% procedure (default: 1e-6)
% - 'fixedlapse': provide fixed lapse parameter (scalar for symmetrical lapse) that will NOT be fitted
% - 'fixedregressor' : value of a fixed regressor with non-fitted weight
% set to 1
% - 'gating': fits weight for gating function
%
% OUTPUT:
% The output variable S is a structure with the following fields:
%     -'plapse': a 2 x nobs matrix representing the posterior probability
%     of a lapse for the corresponding response for each observed value.
%     - 'lh': the likelihood for response 1 according to the model with
%     best-fitting parameters
%      - 'resid': residual
%      - 'w': a vector with best fitting value for each weight
%     - 'bias': scalar with best-fitting value of the constant bias (void
%     if not included in the model)
%     - 'lapse': scalar or vector with best-fitting lapse parameters
%     -'n': number of observed values
%     -'exitflag': exitflag for minimization procedure for overall best
%     minimization iteration (should be positive if minimization suceeded)
%     - 'niter': number of minimization steps of the minimization procedure for
%     overall best minimization iteration
%     - 'covb':  covariance matrix of the parameter (inverse of Fisher
%     observed information matrix ?)
%     -'se': vector of standard errors on the values of parameters
%     - 'T': vector of T-statistics on the values of parameter
%     - 'p': p-value from Wald test for each parameter
%     - 'w_se','bias_se' and 'lapse_se': values of se divided into weight, bias and lapse parameters
%     - 'w_p','bias_p' and 'lapse_p': values of p divided into weight, bias and lapse parameters
%     - 'joint': a matrix with the estimated value, standard-error, and p-value for each parameter in first, second and third column resp.
%     - 'LLH_MAP': log-posterior value for Maximum A Priori parameters
%     - 'LLH': estimated log-evidence for the model (to improve...)
%     - 'BIC': Bayesian Information Criterion
%     - 'AIC': Akaike Information Criterion
%     - 'AICc': corrected Akaike Information Criterion
%     - 'xent': average cross-entropy per observed value
%     - 'certainty': exponential of cross-entropy estimate (certainty of the model)
%     - 'n_init': number of starting points of the minimization procedure
%     - 'LLH_iter': vector of lost-posterior at each minimization step,
%     for the overall best minimization iteration
%     - 'x_all': matrix with MAP parameters obtained for each minimization
%     parameter
%     - 'LLH_all': vector with the value of lost-posterior obtained for
%     each minimization iteration
%     - 'exitflag_all': vector with exitflag for each minimization step
%     - 'S.all_LLH_iter':  cell with the vectors of lost-posterior at each minimization step,
%     for each minimization iteration
%     - 'MinFvalCount': number of minimization iterations that achieved the
%     overall best lost-posterior value
%
%  Example:
%  X = randn(1000,4);
%  Z = X*[.3 -.1 .6 2]' + .3 + randn(1000,1);
%  Y = (Z>0);
%  Y(rand(1,1000)<.05) = 0;
%  Y(rand(1,1000)<.1) = 1;
%  S = probitlapse(X,Y,1)
%
% See also glmfit

function S = probitlapse(X, Y, Ridge,nlapse, x0,param)

if nargin==0
    S = '1.0'; % version
    return;
end

%% process arguments
if (nargin < 6)
    param = [];
end

[nobs, nreg] = size(X); % number of observed values and regressors

if ~isvector(Y)
    error('Y must be a vector');
elseif length(Y) ~= nobs
    error('number of trials in X and Y do not match');
end
Y = Y(:);
yy = [1-Y(:) Y(:) ]; % first column for response 0, second for response 1 (allows also to enter y not as binary but as responsability (real value between 0 and 1, for EM algo)
yy_s = [Y(:)-1 Y(:) ]; % signed version

% add constant vector (for bias)
with_bias = ~isfield(param, 'bias') || param.bias;
if with_bias
    X = [ones(nobs,1) X];
    nreg = nreg + 1;
end

% lapse parameter
if (nargin<4)
    nlapse = 2;
elseif ~isnumeric(nlapse) || length(nlapse)~=1 || ~any(nlapse==[0 1 2])
    error('nlapse must take value 0, 1 or 2');
end
if isfield(param,'fixedlapse')
    lapse = param.fixedlapse;
    do_lapse = 0;
else
    do_lapse = (nlapse>0);
end

% fixed regressors
if isfield(param,'fixedregressor')
    FR = param.fixedregressor;
else
    FR = zeros(nobs,1);
end

%gating
if isfield(param,'gating')
    if ~ismatrix(param.gating) || size(param.gating,1)~=nobs
        error('gating variable must be a matrix with same number of observations as Y');
    end
    ngating = 1+size(param.gating,2); % add constant regressor
else
    param.gating = [];
    ngating = 0;
end
if isempty(param.gating)
    npar = nreg + do_lapse*nlapse; % number of parameters: one for each variable + for lapse
else
    npar = nreg + 2*ngating;
end



% empty vector : return all with nan
if nobs ==0
    S = nullstruct(nreg,npar, with_bias,nlapse);
    return;
end

% by default : no prior
if nargin<3
    Ridge =  zeros(nreg,1);
end
if isfield(param, 'alpha') % Dirichlet prior
    alpha = param.alpha;
else
    hh_alpha = .03;
    switch nlapse
        case 0
            alpha = 0;
        case 1
            alpha = [1 (1/(2* hh_alpha)-1)];
        case 2
            alpha = [ones(1,nlapse) nlapse*(1/hh_alpha-1)]; %
    end
end
alpha = alpha(:);
logB = sum(gammaln(alpha)) - gammaln(sum(alpha)); % log of normalizing factor of Dirichlet (multivariate beta function)

if (length(Ridge) == 1)
    Ridge = Ridge*ones(1,nreg);
elseif length(Ridge(:)) == nreg
    Ridge = Ridge(:)';
else
    error('Ridge weight vector should be length 1 or %d', nreg);
end
if any(Ridge<0)
    error('Ridge values must be positive or null');
end
% log of Ridge scaling factor (ensures that prior sums to 1)
%withRidge = Ridge>0;
k_ridge = 1;
%k_ridge = log ( (prod(Ridge(withRidge))/(2*pi)^sum(withRidge)) )/2;

Ridge_gating = eye(2*ngating);

if (~isfield(param, 'maxiter'))
    param.maxiter = 10000;
end

if (~isfield(param, 'verbose'))
    param.verbose = 0;
end

% tolerance criterion for convergence (on LLH)
if (~isfield(param, 'epsilon'))
    param.epsilon =  1e-6;
end

% exclude data with nan
exclvar = any(isnan(X),1);
if isempty(param.gating)
    exclpar = [exclvar false(1,nlapse*do_lapse)];
else
    exclpar = [exclvar false(1,2*ngating)];
end

X(:,exclvar) = [];
Ridge(exclvar) = [];
nvar_withnan = nreg;
npar_withnan = npar;
nreg = nreg - sum( exclvar);
npar = npar - sum(exclvar);

% Ridge diagonal matrix
Ridgemat =  spdiags(Ridge',0,nreg,nreg);

%% generate initial value for parameters
% number of initial points
if nlapse ==0
    ninit = 1; % without lapse, LLH is convex so always converges to global maxima
elseif isfield(param,'ninit')
    ninit = param.ninit;
else
    ninit = 100;
end
if ((nargin < 5) || (isempty(x0)))
    x0 = zeros(nreg, 1);
    param2 = param;
    param2.gating = [];
    x0 = fitting(X, x0, [],yy, yy_s, Ridge, Ridgemat, [], alpha, 0, 0, FR,logB, param2); % start by fitting without lapse

    if do_lapse
        if ngating>0 % gating variables
            x0 = [x0; zeros(nlapse*ngating,1)]; % initialize gating weights at 0
        else % fixed mixing coeffs
            x0 = [x0;alpha(1:nlapse)/sum(alpha)]; % maximum for prior
        end
    end
else % provided
    if isvector(x0)
        x0 = x0(:);
    end
end

if ninit>1 && size(x0,2)<2
    amax = max(abs(X),[],1); % maximum value for each factor
    maxw = 1./ amax / nreg; % maximum initial values for weight: make sure that no probit reaches 0 or 1
    maxw(amax==0) = 0; % in case there is null factor
    w = bsxfun(@times, maxw', rand(nreg,ninit-1));

    if do_lapse
        if ngating>0
            amax = [1 max(abs(param.gating),[],1)]; % maximum value for each gating regressor
            maxw = 1./ amax / ngating; % maximum initial values for weight: make sure that no probit reaches 0 or 1
            maxw(amax==0) = 0; % in case there is null factor
            lapse = bsxfun(@times, maxw', rand(ngating,(ninit-1)*nlapse));
            lapse = reshape(lapse, ngating*nlapse, ninit-1)';
        else
            lapse = f_gamrnd(alpha(:)',ninit-1); % initial lapse value: sample from prior Dirichlet distribution (gamma distribution and then normalize)
            lapse = lapse(:,1:nlapse) ./ repmat(sum(lapse,2),1,nlapse); % normalize and exclude last term (proba for nolapse)
        end
        %   lapse = .5*rand(nlapse, ninit-1); % all lapse values between 0 and .5
    end
    if do_lapse
        x0(:,2:ninit) = [w;lapse'];
    else
        x0(:,2:ninit) = w;
    end
end




illcond = (npar>nobs);
if illcond
    warnstat1 = warning('OFF', 'MATLAB:nearlySingularMatrix');
    warnstat2 = warning('OFF', 'MATLAB:singularMatrix');

    warning('message','stats:glmfit:IllConditioned');
end

crossvalid = 0;
if isfield(param, 'crossvalidation')
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
            ntrain = nobs-1; % leave one out
        elseif ntrain<1 % defined as percentage
            ntrain = round(ntrain*n);
        end
        if isfield(param, 'ntest') % defined number of test trials
            ntest = param.ntest;
        else % by default all trials not in train set
            ntest = nobs-ntrain;
        end
        if isfield(param, 'nperm')
            nperm = param.nperm;
        else
            nperm = 100; % default number of permutations
        end
    end
end

i_bias = 1:1*with_bias; % index for bias (either one or empty)
i_w = 1+1*with_bias:nreg; % index for variables, excluding bias
if ngating>0
    i_lapse = nreg+(1:nlapse*ngating); % index for lapse variables
else
    i_lapse = nreg+(1:nlapse*do_lapse); % index for lapse variables
end

%% cross validation
if crossvalid

    testscore = zeros(1,nperm); % score for each (LLH per trial)
    accuracy = zeros(1,nperm); % proportion correct
    exitflag = zeros(1,nperm); % exit flag (converged or not)
    w_CV = zeros(nreg,nperm);
    lapse_CV = zeros(nlapse,nperm);

    for p=1:nperm % for each permutation
        if generateperm % generate permutation
            this_ntrain = ntrain;
            trainset = randperm(n,this_ntrain); % training set
            notrainset= setdiff(1:n, trainset); % trials not in training set
            testset = notrainset(randperm(n-ntrain,ntest)); % test set
        else
            traintest = allperm(p,:);
            trainset = traintest{1}; %allperm{p,1};
            testset = traintest{2};%allperm{p,2};
            this_ntrain = length(trainset);
        end

        %fit on training set
        [w,lapse,exitflag(p)] = fitting_multinit(X(trainset,:), ninit, x0, i_lapse,...
            yy(trainset,:), yy_s(trainset,:), Ridge, Ridgemat, Ridge_gating, alpha, nlapse, do_lapse, FR(trainset,:),logB, param);


        %compute score on testing set (mean log-likelihood per trial)
        [testscore(p),accuracy(p)] = loglike(X(testset,:),...
            yy(testset,:),FR(testset,:),w, lapse, param) ;
        testscore(p) = testscore(p)/ sum(w(testset)); % normalize by number of trials

        w_CV(:,p) = w;
        lapse_CV(:,p) = lapse;
    end


    n_nonconverged = sum(exitflag<=0);
    if  n_nonconverged>0
        warning('multilogistic:notconverged', 'Failed to converge for %d/%d permutations', n_nonconverged, nperm);
    end

    % mean over-cross validation sets
    w = mean(w_CV,2);
    lapse = mean(lapse_CV,2);
    S.w_CV = w_CV;
    S.lapse_CV = lapse_CV;


    S.testscore = mean(testscore);
    S.testscore_all = testscore;
    S.accuracy = mean(accuracy);
    S.accuracy_all = accuracy;
    S.exitflag = exitflag;
    S.converged = sum(exitflag>0); % number of permutations with convergence achieved



else % just training on all dataset

    [w,lapse,exitflag,logpost,LLH_MAP,S] = fitting_multinit(X, ninit, x0, i_lapse,...
        yy, yy_s, Ridge, Ridgemat,Ridge_gating, alpha, nlapse, do_lapse, FR,logB, param);
    if  exitflag<=0
        warning('multilogistic:notconverged', 'Failed to converge');
    end
end


%% fit weights and lapse parameters

% residuals
A = X * w + FR; % multiply design matrix by weights
Pm = (1+erf(A/sqrt(2)))/2; % probit function (prob of response 1 without lapse)
if isempty(param.gating)
    switch nlapse
        case 0
            pi_m = 1;
            lh = Pm;
            lhh = [1-lh lh]; % probability of each of the responses
            S.plapse =  nan(2,nobs); % posterior probability that trial was a lapse
        case 1
            pi_m = 1-2*lapse; % probability of non-lapse
            lh = lapse + pi_m*Pm; % likelihood of response 1
            lhh = [1-lh lh]; % probability of each of the responses
            S.plapse =  lapse*( yy ./ lhh)'; % posterior probability for lapse given current parameters
        case 2
            pi_m = 1-sum(lapse); % probability of non-lapse
            lh = lapse(2) + pi_m*Pm;
            lhh = [1-lh lh]; % probability of each of the responses
            S.plapse = bsxfun(@times, lapse', yy ./ lhh)'; % posterior probability for lapse given current parameters
    end
else % mixture of experts
    %  aa = exp(param.gating*lapse); % trial x side of activations
    aa = exp(lapse(1,:) + param.gating*lapse(2:end,:)); % trial x side of activations

    pi_all = aa./(1+sum(aa,2)); % multinomial logit
    pi_m = 1-sum(pi_all,2); % prior probability of no-laspe
    lh = pi_all(:,2) +  pi_m.*Pm;
    lhh = [1-lh lh]; % probability of each of the responses
    S.plapse = (yy .* pi_all ./ lhh)'; % posterior probability for lapse
end


S.lh = lh;
S.resid = Y-lh; % residual

% weights and lapses
S.w = nan(1,nvar_withnan-with_bias);
S.w(~exclvar(1+with_bias:end)) = w(i_w)';

S.bias =  w(i_bias)';
S.lapse = lapse';

S.n = nobs;
S.exitflag = exitflag;
S.niter = sum(S.all_niter); % total number of iterations

%% if converged, compute full Hessian (observed information matrix)
if exitflag>0

    % compute Hessian for weight terms
    z_nl = yy_s ./ lhh; % posterior probability for no-lapse given current parameters
    z_nl(yy==0) = 0; % to avoid nan values for 0/0

    gauss = exp(-A.^2/2)/sqrt(2*pi); % gaussian function
    srat1 = sum(z_nl,2);
    rat2 = yy ./ lhh.^2;
    rat2(yy==0) = 0; % avoid nan values
    srat2 = sum(rat2,2);

    hh = srat1.*A + pi_m .* srat2.*gauss;
    H_ww = - X' *  spdiags( pi_m.*gauss.*hh, 0, nobs, nobs) * X; % Hessian of LLH

    % add Ridge term (Bishop - Pattern recognition and Machine learning - eq 4.143)
    H_ww = H_ww - k_ridge * Ridgemat;

    % add hessian over lapses and lapse-weight
    if isempty(param.gating) % fixed mixing coefficients
        switch nlapse*do_lapse
            case 0
                H = H_ww;
            case 1
                H_ll = + sum((1-2*Pm) .^2 .* srat2) - (alpha(1)-1)/lapse^2 - 4*(alpha(2)-1)/(1-2*lapse)^2;  % second lapse param
                H_wl = (( rat2(:,2) - rat2(:,1)) .* gauss)' * X;   % hessian over first lapse param and wieght
                H = [H_ww H_wl';  H_wl H_ll]; % full Hessian

            case 2
                hh_pip = (alpha(3)-1)/pi_m^2; % hessian of prior over nolapse parameter
                H_ll = - sum(Pm .^2 .* srat2) - (alpha(1)-1)/lapse(1)^2 - hh_pip;  % second lapse param
                H_mm =  - sum((1-Pm) .^2 .* srat2) - (alpha(2)-1)/lapse(2)^2 - hh_pip; %.g_l(:).^2.*yy(:)); % second lapse param
                H_lm =  sum(Pm .*(1-Pm) .* srat2) - hh_pip; % first-second lapse param
                H_lapse = [H_ll H_lm; H_lm H_mm]; % Hessian over lapse parameters

                H_wl = ( (-lapse(2)*srat2 + rat2(:,1)) .* gauss)' * X;   % hessian over first lapse param and wieght
                H_wm = ( (lapse(1)*srat2 - rat2(:,2)) .* gauss)' * X;  % hessian over second lapse param and wieght
                H_wlapse = [H_wl;H_wm];

                H = [H_ww H_wlapse';  H_wlapse H_lapse]; % build full Hessian  matrix
        end
    else
        pi0 = pi_all(:,1); pi1 = pi_all(:,2);
        GG = [ones(nobs,1) param.gating];
        nlh = 1-lh;
        dd =  pi0 .* lh .* ( (2*pi0-1).*srat1 - pi0.*lh.*srat2 );
        H_ll =  GG' *  spdiags(dd, 0, nobs, nobs) * GG;  % second lapse param
        dd =  pi1 .* nlh .* ( (1-2*pi1).*srat1 - pi1.*nlh.*srat2 );
        H_mm =  GG' *  spdiags(dd, 0, nobs, nobs) * GG;  % second lapse param
        dd =  pi0 .* pi1 .* ( (1-2*lh).*srat1 + lh.*nlh.*srat2 );
        H_lm =  GG' *  spdiags(dd, 0, nobs, nobs) * GG;  % first-second lapse param
        H_lapse = [H_ll H_lm; H_lm H_mm]; % Hessian over lapse parameters

        dd = pi0.*pi_m.*gauss.*(1-Y)./nlh.^2;
        H_wl = GG' *  spdiags(dd, 0, nobs, nobs) * X;   % hessian over first lapse param and wieght
        dd = pi1.*pi_m.*gauss.*Y./lh.^2;
        H_wm = GG' *  spdiags(dd, 0, nobs, nobs) * X;   % hessian over second lapse param and wieght
        H_wlapse = [H_wl;H_wm];

        H = [H_ww H_wlapse';  H_wlapse H_lapse]; % build full Hessian  matrix
    end

    %Covariance matrix
    S.covb = nan(npar_withnan);
    invH = inv(-H);
    S.covb(~exclpar,~exclpar) =  invH;

    if any(diag(invH)<0)
        %warning('information matrix is not definite positive... tolerance criterion might be set too high');
    end

    % standard error of estimates
    S.se = sqrt(diag(S.covb))';

    % T-statistic for the weights
    S.T = nan(1,npar_withnan);
    S.T(~exclpar) = [w;lapse(:)]' ./ S.se(~exclpar);

    % p-value for significance of each coefficient (Wald T-test)
    %S.p = 2*normcdf(-abs(S.T))';
    S.p = erfc(abs(S.T)/sqrt(2))'; %use built-in error function to compute normal c.d.f

    % S.p = 1-chi2cdf(S.T.^2,1)';

else % did not converge: cannot estimate observed information matrix

    warning('probitreg:notconverged', 'Failed to converge');

    S.covb = nan(npar_withnan);
    S.se = nan(1,npar_withnan);
    S.T = nan(1,npar_withnan);
    S.p = nan(1,npar_withnan);
end

i_w = 1+1*with_bias:nvar_withnan; % index for variables, excluding bias
if ngating>0
    i_lapse = nvar_withnan+(1:nlapse*ngating); % index for lapse variables
else
    i_lapse = nvar_withnan+(1:nlapse*do_lapse); % index for lapse variables
end

S.w_se = S.se(i_w);
S.bias_se = S.se(i_bias);
S.lapse_se = S.se(i_lapse);
S.w_p = S.p(i_w);
S.bias_p = S.p(i_bias);
S.lapse_p = S.p(i_lapse);
if do_lapse
    S.joint = [[S.bias, S.w, S.lapse(:)']; S.se; S.p(:)']';
else
    S.joint = [[S.bias, S.w]; S.se; S.p(:)']';
end
if ngating>0
    S.lapse_se = reshape(S.lapse_se, ngating,2)'; % weight x response for gating
    S.lapse_p = reshape(S.lapse_p, ngating,2)'; % weight x response for gating
end

S.ridge = Ridge;

% LLH at Max Likelihood / Max A Posteriori parameters
S.LLH_MAP = LLH_MAP;

% model evidence using Laplace approximation (Rasmussen - eq 3.32)  -
% requires that a prior has been defined
if exitflag>0 && all(Ridge>0)
    [CH, flag] = chol(-H);
    if flag
        S.LLH = nan;
    else
        logdetH = 2*sum(log(diag(CH))); % log-determinant of -H (fast reliable way)
        logdetprior = -sum(log(Ridge));
        S.LLH = logpost - (logdetH+logdetprior)/2 ;
    end
else
    S.LLH = nan;
end
S.BIC = npar*log(nobs) -2*LLH_MAP; % Bayes Information Criterion
S.AIC = 2*npar - 2*LLH_MAP; % Akaike Information Criterior
S.AICc = S.AIC + 2*npar*(npar+1)/(nobs-npar-1); % AIC corrected for sample size

S.xent = S.LLH/nobs; % cross-entropy estimate
S.certainty = exp(S.xent); % exponential of cross-entropy estimate (certainty of the model)

S.n_init = ninit; % number of initial points
%S.LLH_iter = SS.logpost_iter;
S.x_all = zeros(npar_withnan, ninit);
S.x_all(~exclpar,:) = S.all_x;
S = rmfield(S,'all_x');


% compute how many initial points arrived at the (overall) minimum
% objective function value
Fvaldiff = bsxfun(@minus,logpost, S.all_logpost); % difference between final point for each starting point and overall final point
S.MinFvalCount = sum(Fvaldiff<param.epsilon);


if illcond
    warning(warnstat1.state, 'MATLAB:nearlySingularMatrix');
    warning(warnstat2.state, 'MATLAB:singularMatrix');
end

end

%% compute predicting accuracy and LLH for given parameter set
function         [testscore,accuracy] = loglike(X, yy,FR,w, lapse,nlapse, param)

A = X * w + FR; % multiply design matrix by weights
Pm = (1+erf(A/sqrt(2)))/2; % probit function (prob of response 1 without lapse)

% compute likelihood of each trial for given parameters
lhh = likeli(Pm, nlapse, lapse, param);

ylh = yy .* log(lhh);
ylh(yy==0) = 0;
testscore = sum(ylh(:));

lh_o = sum(yy .*lhh); % proba of observed responses
accuracy = mean(lh_o>.5); % accuracy
end

%% fitting algo (multiple starting points)
function [w,lapse,exitflag,logpost,LLH_MAP,S]  = fitting_multinit(X, ninit, x0, i_lapse,...
    yy, yy_s, Ridge, Ridgemat, Ridge_gating,alpha, nlapse, do_lapse, FR,logB, param)

npar = size(x0,1);
nreg = size(X,2); %npar - nlapse*do_lapse;

%initialize variables
S.all_x = zeros(npar,ninit);
S.all_logpost = zeros(1,ninit);
S.all_LLH_MAP = zeros(1,ninit);
S.all_LLH_iter = cell(1,ninit);
S.all_exitflag = zeros(1,ninit);
S.all_niter =  zeros(1,ninit);



for n = 1:ninit
    % fprintf('%d.',n);
    w = x0(1:nreg,n); % weights for each variable
    if do_lapse || nlapse==0
        lapse = x0(i_lapse,n); % lapse
    else
        lapse = [];
    end

    % find best-fitting parameters from starting point
    [w,lapse,converged,iter,logpost,LLH,logpost_iter] = fitting(X, w, lapse,...
        yy, yy_s, Ridge, Ridgemat, Ridge_gating,alpha, nlapse, do_lapse, FR,logB, param);

    % put that in big structure
    if isempty(param.gating)
        S.all_x(:,n) = [w;lapse(1:nlapse*do_lapse)];
    else
        S.all_x(:,n) = [w;lapse(:)];
    end
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
w = x(1:nreg); % weights for each variable
if do_lapse
    lapse = x(i_lapse); % lapse
    if ~isempty(param.gating)
        lapse = reshape(lapse, length(lapse)/2,2); % weight x response for gating
    end
else
    lapse = lapse';
end
exitflag = S.all_exitflag(imax);
S.logpost_iter = S.all_LLH_iter{imax};

end


%% fitting algorithm (single starting point)
function  [w,lapse,converged,iter,logpost,LLH,logpost_iter] = fitting(X, w, lapse,yy, yy_s, Ridge, Ridgemat, Ridge_gating, alpha, nlapse, do_lapse, FR,logB, param)

A = X * w + FR; % multiply design matrix by weights
Pm = (1+erf(A/sqrt(2)))/2; % probit function (prob of response 1 without lapse)

%oldproby = -ones(size(y));
old_logpost = -Inf;
logpost_iter = zeros(1,param.maxiter);
still =1;
iter = 0;
k_ridge = 1;
nobs = size(X,1);

G =[ones(nobs,1) param.gating];
ngating = size(G,2);

%% let's start by lapse update to see if it helps
if do_lapse

    if    ~isempty(param.gating)
        lapse = reshape(lapse,length(lapse)/2,2);
    end

    %% E-step again 2: estimate probability of lapse and no-lapse for each trial

    % compute likelihood for each trial
    [lhh,pi_all] = likeli(Pm, nlapse, lapse, param);

    z_l = bsxfun(@times, pi_all(:,1:2), yy ./ lhh); % posterior probability for lapse given current parameters
    %z_l = pi_all(:,1:2) .* yy ./ lhh; % posterior probability for lapse given current parameters


    %% M-step 2: update lapse probability only (and compute LLH for new parameters)
    if isempty(param.gating)
        switch nlapse
            case 1 % symmetrical lapse
                lapse = (alpha(1)-1 + sum(z_l(:)))/(nobs + sum(alpha)-2)/2;
            case 2 % asymmetrical lapse
                lapse = (sum(z_l,1)' + (alpha(1:2)-1))/(nobs+sum(alpha)-3);
                if lapse<0
                    warning('ede');
                end
                lapse = min(lapse,1); % keep lapses between 0 and 1 ( only appears if alpha<0)
                lapse = max(lapse,0);
        end
    else
        zz_l = [z_l 1-sum(z_l,2)]; % three columns: post proba for rsep 1, resp2 and model
        % lapse = mnrfit(param.gating,zz_l); % fit multiclass logistic model
        lapse = mnr_update(G, zz_l,lapse, ngating, nobs,k_ridge, Ridge_gating);
    end
end

% for EM convergence speed-up
thr_FHU = [.9 .9 1.1];
prev_lapse = nan(size(lapse));
FastJump = false;
LastFastJump = 0;
%  dLapse = nan(size(lapse));
% dLogpost = nan;



%% ITERATION PROCEDURE
while still

    %% E-step 1 : estimate probability of lapse and no-lapse for each trial

    % compute likelihood for each trial
    [lhh, pi_all] = likeli(Pm, nlapse, lapse, param);
    pi_m = pi_all(:,3);

    z_nl = yy_s ./ lhh; % posterior probability for no-lapse given current parameters
    z_nl(yy==0) = 0; % to avoid nan values for 0/0

    %% M-step 1: update weights only

    % compute gradient over weights
    gauss = exp(-A.^2/2)/sqrt(2*pi); % gaussian function
    srat1 = sum(z_nl,2);
    dd = pi_m .* srat1.*gauss;
    grad_w = dd' * X; % gradient of LLH with respect to weights

    % add Ridge term
    grad_w = grad_w -  k_ridge * Ridge .* w';

    % compute Hessian for weight terms
    rat2 = yy ./ lhh.^2;
    rat2(yy==0) = 0; % avoid nan values
    srat2 = sum(rat2,2);

    hh = srat1.*A + pi_m .* srat2.*gauss;
    H_ww = - X' *  spdiags( pi_m.*gauss.*hh, 0, nobs, nobs) * X; % Hessian of LLH

    % add Ridge term (Bishop - Pattern recognition and Machine learning - eq 4.143)
    H_ww = H_ww - k_ridge * Ridgemat;

    %Newton-Raphson update on weights
    w = w - H_ww\grad_w';

    if rcond(H_ww)<1e-16
        2;
    end


    A = X * w + FR; % multiply design matrix by weights
    Pm = (1+erf(A/sqrt(2)))/2; % probit function (prob of response 1 without lapse)

    if do_lapse
        %% E-step again 2: estimate probability of lapse and no-lapse for each trial

        % compute likelihood for each trial
        [lhh,pi_all] = likeli(Pm, nlapse, lapse, param);
      
        z_l = bsxfun(@times, pi_all(:,1:2), yy ./ lhh); % posterior probability for lapse given current parameters

        %% M-step 2: update lapse probability only (and compute LLH for new parameters)
        if isempty(param.gating)
            switch nlapse
                case 1 % symmetrical lapse
                    lapse = (alpha(1)-1 + sum(z_l(:)))/(nobs + sum(alpha)-2)/2;
                case 2 % asymmetrical lapse
                    lapse = (sum(z_l,1)' + (alpha(1:2)-1))/(nobs+sum(alpha)-3);
                    if lapse<0
                        warning('ede');
                    end
                    lapse = min(lapse,1); % keep lapses between 0 and 1 ( only appears if alpha<0)
                    lapse = max(lapse,0);
            end
        else
            zz_l = [z_l 1-sum(z_l,2)]; % three columns: post proba for rsep 1, resp2 and model

            if 0
                lapse = mnrfit(param.gating,zz_l); % fit multiclass logistic model
            else
                %% one Newton-Raphson update on multiclass logistic regression
                lapse = mnr_update(G, zz_l,lapse, ngating, nobs, k_ridge, Ridge_gating);

              
            end
        end

    end


    % compute likelihood for each trial
    lhh = likeli(Pm, nlapse, lapse, param);

    %% log-likelihood of observed data
    ylh = yy .* log(lhh);
    ylh(yy==0) = 0;
    LLH = sum(ylh(:));

    % posterior (with prior)
    logpost = LLH - k_ridge * Ridge*w.^2/2;
    if isempty(param.gating)
        switch nlapse*do_lapse
            case 1 % beta prior for symmetrical lapses
                fullapse = [lapse*2; 1-2*lapse];
            case 2 % Dirichlet prior for two lapse
                fullapse = [lapse; 1-sum(lapse)];
        end
        if any(nlapse*do_lapse==[1 2])
            alpha_lapse= (alpha-1).*log(fullapse);
            alpha_lapse(alpha==1) = 0; % making sure to avoid 0*inf
            logpost = logpost + sum( alpha_lapse) - logB;

        end
    end

    % check whether should iterate more
    iter = iter + 1; % update iteration counter;

    % has converged if improvement in LLH is smaller than epsilon
    dLogpost= old_logpost-logpost;
    converged = abs(dLogpost)<param.epsilon;

    if FastJump && logpost<old_logpost
        % if fast jump failed, revert to lapse values before fast jump
        lapse = prev_lapse- dLapse;
    end

    old_logpost = logpost;
    logpost_iter(iter) = logpost;
    dLapse = lapse - prev_lapse; % lapse updates

    still = (iter<param.maxiter) && ~converged && ~any(isnan(w)) && ~any(isnan(lapse(:)));



    if do_lapse && still && iter>2

        cos_successive_updates(iter) = dot(prev_dLapse,dLapse)/(norm(prev_dLapse)*norm(dLapse));
        rho_2 =  norm(dLapse)^2/norm(prev_dLapse)^2;
        rat_successive_updates(iter) =rho_2;
        dLapse_hist(:,iter) = dLapse;

        rho = dLogpost / prev_dLogpost;
        rat_successive_logjoint(iter) = rho;
        consistency = rho/  rho_2; % make sure the rho estimate is similar to that from series of lapse

        % whether we're in quadratic region of full
        % parameter space
        inQuadratic=  cos_successive_updates(iter)>thr_FHU(1) && consistency>thr_FHU(2) && consistency<thr_FHU(3);

        FastJump =  inQuadratic && rho>.8 && rho<1 && (iter>LastFastJump+10);

        % look for convergence point if in quadratic region
        % and would converge slowly
        if FastJump
            rho_0 = min(mean([rho, rho_2]),.99);

            % convergence point of time series
            tmp = lapse;
            lapse = lapse + dLapse/(1-sqrt(rho_0))*sqrt(rho_0);
            if any(lapse<0)
                lapse = tmp;
                FastJump = 0;
            end

            LastFastJump = iter;
        end


    end
    prev_lapse = lapse;
    prev_dLapse = dLapse;
    prev_dLogpost = dLogpost;


end

logpost_iter(iter+1:end) = [];

end



% null structure (for void data)
function S = nullstruct(nvar,npar, with_bias,nlapse)
S.plapse = zeros(2,0);
S.lh = zeros(0,1);
S.resid = nan(0,1);
S.w = nan(1,nvar-with_bias);
S.bias =  nan(1,1*with_bias);
S.lapse = nan(1,nlapse);
S.n = 0;
S.exitflag = 0;
S.niter = 0; % total number of iterations
S.covb = nan(npar);
S.se = nan(1,npar);
S.T = nan(1,npar);
S.p = nan(1,npar);
S.w_se = nan(1, nvar-with_bias);
S.bias_se = nan(1,1*with_bias);
S.lapse_se = nan(1,nlapse);
S.w_p = nan(1, nvar-with_bias);
S.bias_p = nan(1,1*with_bias);
S.lapse_p = nan(1,nlapse);
S.joint = [[S.bias, S.w, S.lapse]; S.se; S.p]';
S.ridge = nan(1,npar);
S.LLH_MAP = 0;
S.LLH = 0;
S.BIC = 0;
S.AIC = 0;
S.AICc = 0;
S.xent = nan; % cross-entropy estimate
S.certainty = nan; % exponential of cross-entropy estimate (certainty of the model)
S.n_init = 0;
S.LLH_iter = [];
S.x_all = nan(npar, 1);
S.LLH_all = [];
S.exitflag_all = [];
S.all_LLH_iter = [];
S.MinFvalCount = nan;

end

%% compute likelihood of each trial for given parameters
function [lhh, pi_all] = likeli(Pm,nlapse, lapse, param)
if isempty(param.gating) % fixed mixing coefficients
    switch nlapse
        case 0
            %  lh = Pm;
            pi_all = [0 0 1]; % prior probability of using model (i.e. no-lapse)
        case 1
            pi_all = [lapse lapse 1-2*lapse]; % prior probability of no-laspe
            %  lh = lapse + pi_*Pm; % likelihood of response 1
        case 2
            pi_all = [lapse(:)' 1-sum(lapse)];
            %  lh = lapse(2) + pi_m*Pm;
    end
    lh = pi_all(2) + pi_all(3)*Pm;
else % gating model
    aa = exp(lapse(1,:) + param.gating*lapse(2:end,:)); % trial x side of activations

    pi_all = aa./(1+sum(aa,2)); % multinomial logit
    pi_m = 1-sum(pi_all,2); % prior probability of no-laspe
    pi_all(:,end+1) = pi_m;
    lh = pi_all(:,2) +  pi_m.*Pm;
end

if any(pi_all(:)<0)
    2;
end

lhh = [1-lh lh]; % probability of each of the responses

end

%% one Newton-Raphson update on multiclass logistic regression
function  lapse = mnr_update(G, zz_l,lapse, ngating, nobs, k_ridge, Ridge_gating)

AA = G *lapse; % multiply design matrix by lapse weights
Y = [exp(AA) ones(nobs,1)]; Y = Y./sum(Y,2); %multilogistic model

% compute gradient over weights
grad = G'*(zz_l(:,1:end-1)-Y(:,1:end-1));

% add Ridge term
grad = grad(:)' -  k_ridge * diag(Ridge_gating)' .* lapse(:)';

% compute Hessian
H = zeros(ngating,ngating,2,2);
for i=1:2
    for j=1:i
        R = Y(:,i).*((i==j)-Y(:,j));
        H(:,:,i,j) = - G' * spdiags( R, 0, nobs, nobs) * G; % Hessian of LLH w.r.t set wieghts for class i and j
    end
    H(:,:,1:i,i) = H(:,:,i,1:i); % symmetric
end
H = permute(H,[1 3 2 4]);
H = reshape(H,ngating*2,ngating*2);

% add Ridge term (Bishop - Pattern recognition and Machine learning - eq 4.143)
H = H - k_ridge * Ridge_gating;

%Newton-Raphson update on weights
lapse = lapse(:) - H\grad';
lapse = reshape(lapse,ngating,2);

end

%% sample from Gamma distribution
function R = f_gamrnd(A,m)
B = 1;
if exist('gamrnd') ==2 % if stats toolbox with gamrnd function
    R = gamrnd(repmat(A,m,1),1,m,length(A));

else % sample using Marsaglia's simple transformation-rejection method
    % see https://stackoverflow.com/questions/41967117/how-to-draw-random-numbers-from-a-gamma-distribution-without-the-statistics-tool

    R = zeros(m,length(A));
    for i=1:length(A)

        R2 = [];
        d = A(i) - 1/3;

        % Marsaglia's simple transformation-rejection:
        while length(R2) < m % generate until we have N values
            x = randn(m,1);
            U = rand(m,1);
            v = (1+x./sqrt(9*d)).^3;
            accept = log(U)<(0.5*x.^2+d-d*v+d*log(v));
            accept = accept & ((d.*v)>0); % missing that sample must be positive
            new = d.*(v(accept)).*B;

            R2 = [R2; new];
        end
        R2(m+1:end) = [];

        R(:,i) = R2;

    end
end
end