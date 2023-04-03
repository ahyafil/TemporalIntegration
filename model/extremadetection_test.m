% Test fitting of extrema detection model on synthetic data
%
% See also extremadetection_fit, extremadetection_rnd, snapshot_test
clear; close all;
nTrial = 10000;  % number of trials
nSample = 8; % number of samples
LastSample = 0; % if threshold is not reached during sequence, whether response is based on last sample (otherwise just chance)
FixedThreshold = false; %fixed threshold over parameter, or varies across the sequence

% threshold parameters
T = 4; %  decision threshold
if ~FixedThreshold % varying threshold across sequence
    stdThreshold = 1.5;
    T = T + stdThreshold*randn(1,nSample);
    T = max(T,0);
end
sigma = 1; % std of gaussian noise on threshold

% lapse rate
lapse =  [.04 .02]; % for each side

% properties of generative model of sensory information (drawn from a
% gaussian)
muX = 2; % mean shift due to stimulus category
sigmaX = 3.5; % standard deviation of sample evidence
Xmax = 10; % maximum value of sample

%% 1. GENERATE STIMULUS SEQUENCE
categ = sign(rand(nTrial,1)-.5); % stimulus category (-1/1)
X = round(muX*categ + sigmaX*randn(nTrial,nSample)); % design matrix
X(abs(X)>Xmax) = Xmax * sign( X(abs(X)>Xmax) );

%% 2. SIMULATE EXTREMA DETECTION MODEL
Y = extremadetection_rnd(X,T, sigma, LastSample, lapse);

%% 3. FIT EXTREMA DETECTION MODEL ON SIMULATED DATA
param.LastSample = LastSample;
param.Lapse = true; % include lapses
param.FixedThreshold = FixedThreshold;
[ T_hat, sigma_hat, lapse_hat, S] = extremadetection_fit(X, Y, param);

%% 4. PLOT
figure;

% plot threshold parameters
subplot(131); title('threshold'); hold on;
bar(T); % true parameters
errorbar(T_hat(:), S.T_se(:),'.r'); % inferred parameters
if ~FixedThreshold
    xlabel('sample');
end
legend({'true','estimated'});

% plot noise
subplot(132); title('noise'); hold on;
bar(sigma); % true parameters
errorbar(sigma_hat, S.sigma_se,'.r'); % inferred parameters
xticks([]); ylabel('\sigma');

% plot lapse parameters
subplot(133); title('lapses'); hold on;
bar(lapse); % true parameters
errorbar(lapse_hat', S.lapse_se,'.r'); % inferred parameters
xticks(1:2); xticklabels({'resp 0','resp 1'});



