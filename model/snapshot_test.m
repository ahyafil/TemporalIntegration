% Test fitting of snapshot model on synthetic data
%
% See also snapshot_fit, snapshot_rnd, extremadetection_test

clear; close all;

% model type
mode ='stochastic';  % 'deterministic' or 'stochastic' versions of snapshot model
span = 1; % span only 1 frame (NOT a parameter at the momdent)
beta = 3; % sensitivity to sample information (for non-deterministic model)

nSample = 7; % number of samples in each stimulus
nTrial = 5000; % number of trials

% properties of generative model of sensory information (drawn from a
% gaussian)
muX = 1; % mean shift due to stimulus category
sigmaX = 3; % standard deviation of sample evidence

is_deterministic = strcmp(mode,'deterministic');
with_lapse = is_deterministic; % use lapses only for deterministic model

pi = drchrnd(ones(1,nSample+2*with_lapse),1); % probability of selecting each frame (including lapses as last two values)
if ~is_deterministic
    %  beta = gamrnd(3, 1, 1, nSample); % sensitivity for attended frame
    beta = beta*ones(1,nSample);
end

%% 1. GENERATE STIMULUS SEQUENCE
categ = sign(rand(nTrial,1)-.5); % stimulus category (-1/1)
X = round(muX*categ + sigmaX*randn(nTrial,nSample)); % design matrix


%% 2. SIMULATE MODEL
if is_deterministic
    [Y,attended] = snapshot_rnd(X,pi);
else
    [Y,attended] = snapshot_rnd(X,pi,beta);
end

%% 2. FIT SNAPSHOT MODEL
param = struct;
param.lapse = with_lapse;
[beta_hat, pi_hat, S]= snapshot_fit( X, Y, mode,span,  param) ;

%% 3. PLOT
figure;

% plot mixture coefficients
if ~is_deterministic
    subplot(121);
end

hold on;
plot(pi,'k'); % true value
errorbar(pi_hat,S.mixt_se', 'r.'); % estimated parameter
xlabel('sample'); ylabel('coefficients');
if with_lapse
xticks(1:nSample+2); xticklabels([string(1:nSample) "lapse0" "lapse1"]);
end
legend({'true','estimated'});

% plot weights
if ~is_deterministic
    subplot(122);

    hold on;
    plot(beta,'k');
    errorbar(beta_hat(end,:)',S.w_se, 'r.');
    xlabel('sample'); ylabel('sensitivity');
end