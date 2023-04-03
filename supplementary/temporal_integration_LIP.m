% Yates et al data: testing temporal integration in LIP responses

clear; close all;


nSample = 7;
EarlySamples = 1:3;
LateSamples = 4:7;
dx = .02;  % bin size (evidence space)
bnd = .8; % boundary of evidence space

nConv = 20; % length of convolution kernel (in bins)
sigma = .1; % width of kernel (in evidence space)

% we focus on neural activity on 500 ms after stimulus offset
WindowAfterStimulusOffset = [0 .5];

nModel = 3;
ModelLabel = {'integration','snapshot', 'extrema detection'};

IntegrationMapColour = [103 169 221; 241 135 34]/255; % gradient between red and blue


%% load spike count data
csvfile = fullfile(TemporalIntegrationRootDirectory, 'data', 'LIP.csv');
T = readtable(csvfile);

% merge columns for stimulus into single variable (one row vector for each
% trial)
T = mergevars(T,"Stimulus_"+(1:nSample),'NewVariableName',"Stimulus");

%% get model parameters from extrema detection model

% directory where model fits are stored
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

md = 'extremadetection';
ExtremaDetectionFile = fullfile(ModelfitsDir,md,sprintf('%s_%s','monkey', md));
ED = load(ExtremaDetectionFile);
T_ED = ED.Tsigma_all(:,1,1); % threshold
sigma_ED = ED.Tsigma_all(:,2,1); % sensory noise

%%

nNeuron = max(T.NeuronIndex);

% edges for
binEdges = -bnd:dx:bnd;
binCenters = [binEdges(1)-dx/2 binEdges+dx/2];

binEdges = [-Inf binEdges Inf]; % bin boundaries


targPref_all = [];

GLM_Weight = zeros(nSample+1, nNeuron);
GLM_WeightSe = zeros(nSample+1, nNeuron);
PreferredDir = zeros(1, nNeuron);
BaseSpikeCount = zeros(1, nNeuron);
nTrial = zeros(1, nNeuron);

for n=1:nNeuron % for each session
    fprintf('neuron #%d/%d,',n, nNeuron);


    Tneuron = T(T.NeuronIndex==n,:);

    nTrial(n) = height(Tneuron);

    SpikeCount = Tneuron.SpikeCount; %spike counts in window
    Stimulus = Tneuron.Stimulus; % stimulus (motion in each pulse)
    Monkey = Tneuron.Monkey{1}; % animal

    % weights of pulses on spike count
    [w, ~, outp] = glmfit(Stimulus, SpikeCount, 'poisson');

    GLM_Weight(:,n) = w;
    GLM_WeightSe(:,n) = outp.se;

    % early and late (weighted) evidence
    EarlyEvidence = Stimulus(:,EarlySamples)*w(1+(EarlySamples));
    LateEvidence = Stimulus(:,LateSamples)*w(1+(LateSamples));
    EarlyLateEvidence = [EarlyEvidence LateEvidence];

    % generate spikes from integration model (Poisson distr with rate
    % according to weighted sum)
    SpkMod = zeros(nTrial(n),nModel);
    PreModIntegration = w(1)+EarlyEvidence+LateEvidence;
    SpkMod(:,1) = poissrnd(exp(PreModIntegration)); % integration model

    % generate spikes from snapshot model
    PreferredDir(n) = sign(sum(w(2:end))); % preferred direction of neuron
    WeightUnsigned = w*PreferredDir(n); % sign with respect to preferred dir
    WeightUnsigned(WeightUnsigned<0) = 0; % if negative, never "attend" that sample
    SumWeights = sum(WeightUnsigned(2:end));
    MixtureCoeffs = WeightUnsigned(2:end)/SumWeights; % mixture coefficients
    AttendedSamples = mnrnd(1,MixtureCoeffs,nTrial(n)); % generate attended sample for each trial
    AttendedSampleVal = sum(AttendedSamples.*Stimulus,2); % attended sample
    PreModSnapshot = w(1) + SumWeights*PreferredDir(n)*AttendedSampleVal;
    SpkMod(:,2) = poissrnd(exp(PreModSnapshot)); % generate spike counts

    % generate spikes from extrema detection
    i_monkey = 1+(Monkey=='n'); % which animal

    % simulate extrema detection with parameters fitted for corresponding monkey
    [X_samp,i_samp] = extrema_detection_sample(Stimulus, T_ED(i_monkey), sigma_ED(i_monkey));
    PreModExtrema = w(1) + 0.5*SumWeights*X_samp*PreferredDir(n);
    SpkMod(:,3) = poissrnd(exp(PreModExtrema));

    % count number of trials in each bin
    nT = grpsum(ones(size(EarlyLateEvidence,1),1), {EarlyLateEvidence(:,1) EarlyLateEvidence(:,2)}, 'dim',1,'histc',{binEdges,binEdges});

    % sum of spike counts in each bin
    SumSpikeCountBin = grpsum(SpikeCount, {EarlyLateEvidence(:,1) EarlyLateEvidence(:,2)}, 'dim',1,'histc',{binEdges,binEdges});

    % sum of simulated spike counts in each bin
    for m=1:nModel
        SumSpikeCountBinModel(:,:,m) = grpsum(SpkMod(:,m), {EarlyLateEvidence(:,1) EarlyLateEvidence(:,2)}, 'dim',1,'histc',{binEdges,binEdges});
    end


    BaseSpikeCount(n) = exp(w(1)); % predicted baseline spike count

    SumSpikeCountBinAll(:,:,n) = SumSpikeCountBin;
    SumSpikeCountBinModelAll(:,:,:,n) = SumSpikeCountBinModel;
    nTrialBin(:,:,n) = nT;

    fprintf('...done\n');
end

%% AVERAGE ACROSS POPULATION
nTrialBinTotal = sum(nTrialBin,3);

% normalize sum of spike by baseline rate for each neuron, then sum across
% neurons
SumSpikeCountBinNorm = sum( SumSpikeCountBinAll./permute(BaseSpikeCount,[1 3 2]) ,3);

% smooth by convolution with gaussian kernel  (one d)
ConvolutionKernel = exp(-(dx*(-nConv:nConv)).^2/2/sigma^2);
ConvolutionKernel2D = ConvolutionKernel'*ConvolutionKernel; % 2d kernel
SumSpikeCountBinNorm = conv2(SumSpikeCountBinNorm,ConvolutionKernel2D,'same');

nTrialBinTotal = conv2(nTrialBinTotal,ConvolutionKernel'*ConvolutionKernel,'same'); % corresponding number of trial
MeanSpikeCountBin = SumSpikeCountBinNorm./nTrialBinTotal;

% repeat the same for simulated spike counts for each model
for m=1:nModel

    % normalize sum of spike by baseline rate for each neuron, then sum across
    % neurons
    SumSpikeCountBinModelNorm = sum(permute(SumSpikeCountBinModelAll(:,:,m,:),[1 2 4 3])./permute(BaseSpikeCount,[1 3 2]),3);
    SumSpikeCountBinModelNorm = conv2(SumSpikeCountBinModelNorm,ConvolutionKernel'*ConvolutionKernel,'same');
    MeanSpikeCountBinModel(:,:,m) = SumSpikeCountBinModelNorm./nTrialBinTotal;
end

MeanSpikeCountBin(isnan(MeanSpikeCountBin)) = 0;
MeanSpikeCountBinModel(isnan(MeanSpikeCountBinModel)) = 0;

% normalize
norm_factor = 2; max([MeanSpikeCountBin(:);MeanSpikeCountBinModel(:)]);
MeanSpikeCountBinNorm = MeanSpikeCountBin/ norm_factor;
MeanSpikeCountBinModelNorm = MeanSpikeCountBinModel/norm_factor;

%% PLOT FIGURE
figure();

% hue intensity
Hue = min(sum(SumSpikeCountBinNorm,3)/20,1);

% experimental spike count panel
subplot(2,2,1);

% convert mean response and number of trials to RGB
for c=1:3
    RGB(:,:,c) = 1- (IntegrationMapColour(1,c)-IntegrationMapColour(2,c)) * MeanSpikeCountBinNorm - IntegrationMapColour(2,c);
end
RGB = 1- RGB.* Hue;

% image plot
image(binCenters, binCenters, permute(RGB,[2 1 3])); axis xy; hold on;

% add contour
isoLines = [.4 .6 1.4 1.8];
contour(binCenters, binCenters, MeanSpikeCountBin', isoLines,'k', 'linewidth',.5); % draw contour lines
contour(binCenters, binCenters, MeanSpikeCountBin',[1 1],'k', 'LineWidth',1); % draw contour lines
xlabel('early evidence'); ylabel('late evidence');
title('LIP neurons');

% simulated spike counts panels
for m=1:nModel
    subplot(2,2,1+m);
    for c=1:3
        RGB_Model(:,:,c) = 1- (IntegrationMapColour(1,c)-IntegrationMapColour(2,c)) * MeanSpikeCountBinModelNorm(:,:,m) - IntegrationMapColour(2,c);
    end
    RGB_Model = 1- RGB_Model.* Hue;
    image(binCenters, binCenters, permute(RGB_Model,[2 1 3])); axis xy; hold on;
    contour(binCenters, binCenters, MeanSpikeCountBinModel(:,:,m)', isoLines,'k', 'linewidth',.5); % draw contour lines
    contour(binCenters, binCenters, MeanSpikeCountBinModel(:,:,m)',[1 1],'k', 'LineWidth',1); % draw contour lines
    xlabel('early evidence'); ylabel('late evidence');
    title(ModelLabel{m});
end

%%%
function [X_samp,i_samp] = extrema_detection_sample(X, mu, sigma)

[nTrial,nSamp] = size(X);

still = true(nTrial,1);
X_samp = zeros(nTrial,1);
i_samp = zeros(nTrial,1);

for s=1:nSamp % for each sample
    Xnoisy = X(:,s) + sigma*randn(nTrial,1); % noise stim

    ReachedThreshold = still & (abs(Xnoisy)>mu); % reaches negative threshold
    i_samp(ReachedThreshold) = s;
    X_samp(ReachedThreshold) = Xnoisy(ReachedThreshold);

    still(ReachedThreshold) = false;
end

% if threshold is not reached
X_samp(still) = Xnoisy(still); % reaches negative threshold
i_samp(still) = s;

end
