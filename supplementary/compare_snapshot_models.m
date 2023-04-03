%% SUPPLEMENTARY FIGURE 2: COMPARE VARIANTS OF THE SNAPSHOT MODEL
%(deterministic vs smooth, fixed vs free lapse, span length) against the temporal integration model

clear; close all;
animal = 'monkey';

nLapseType = 2; % free or fixed lapses
nBehType = 2; % probabilistic or deterministic mapping
nType = nLapseType*nBehType; % four different model types
LapseTypeLabel = {'free','fixed'};
BehTypeLabel = {'deterministic','probabilistic'};
AllTypeLabel = {'free det','fixed det','free proba','fixed proba'};
ColorMod = repelem({'b','k'},nLapseType);    % blue: deterministic; black: probabilistic
LineStyle = repmat({'-','--'},1,nBehType); % full line: free lapses; dashed: fixed lapses

nQuantile = 50; % number of quantiles for psychometric curves

%% LOAD AND PREPROCESS DATA
csvfile = fullfile(TemporalIntegrationRootDirectory, 'data', animal);
T = readtable([csvfile '.csv']);
fprintf('loading behavioral data from %s\n',csvfile);

% preprocess data
[T, nSample, subject_id, nSubject, nSamples_max, nSamples_mean] = preprocess_table(T);

nSpan = nSample-1; % model covering just 1 frame to all but one sample

% all model names
Models = cell(nSpan, nLapseType, nBehType); % span x lapse (free, fixed) x (deterministic/probabilistic)
for sp=1:nSpan
    for lt = 1:nLapseType
        for bt = 1:nBehType
            mdl = ['snapshot' num2str(sp)]; % add span number
if bt ==2 % 1: deterministic (non-probabilistic); 2: probabilistic
                mdl = [mdl '_stochastic'];
            end
            if lt ==1 % 1: free lapses; 2: fixed lapses
                mdl = [mdl '_freelapse'];
            end
            
            Models{sp,lt,bt} = mdl;
        end
    end
end
nModels = numel(Models);

ModelColor = num2cell([0 0 1] + (0:nSpan-1)'/(nSpan+2).*[1 1 0],2); % dark blue for span = 1 to light blue for largest span

% model files
ModelFiles = cell(size(Models));
for m=1:nModels
    ModelFiles{m}  = sprintf('%s/%s_%s',Models{m}, animal,Models{m});
end

% where to load analysis from from
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

% integration model file, used as reference
integration_model = 'integration';
IntegrationFile  = sprintf('%s/%s_%s',integration_model, animal,integration_model);
Sref = load(fullfile(ModelfitsDir,IntegrationFile));

% integration map correlation for integration model
corr2dmap_ref = [Sref.M.score.IntegrationMap.r];

%% COMPUTE ACCURACY, PSYCHOMETRIC CURVES AND KERNELS FOR EXPERIMENTAL AND SIMULATED DATA, AIC AND CORRELATION BETWEEN INTEGRATION MAPS

Accuracy = zeros(1,nSubject);
AccuracyModel = zeros(nSpan,nLapseType,nBehType, nSubject);
AccuracyModelSem = zeros(nSpan,nLapseType,nBehType, nSubject);
StimulusEvidence = cell(1,nSubject);
PK_Data = zeros(nSample+1,nSubject);
PK_DataSe = zeros(nSample+1,nSubject);
PK_Model = zeros(nSample+1,nSpan, nLapseType,nBehType,nSubject);
PK_ModelSe = zeros(nSample+1,nSpan, nLapseType,nBehType,nSubject);
AIC = zeros(nSpan, nLapseType, nBehType,nSubject);
DeltaAIC = zeros(nSpan, nLapseType, nBehType,nSubject);
CorrelationIntegrationMap  = zeros(nSpan, nLapseType, nBehType,nSubject);
for s=1:nSubject
    % trial corresponding to this subject
    mask = T.subject== subject_id(s);
    resp = T.resp(mask);
    target = T.target(mask);
    stimulus = T.stimulus(mask,:);
    nTrial = sum(mask);

    % mean accuracy of animal
    Accuracy(s) = mean(resp == target);

    % animal sensory kernel (logistic regression)
    [PK_Data(:,s),~,LRstats] = glmfit(stimulus,resp, 'binomial');
    PK_DataSe(:,s) = LRstats.se; % standard error

    StimulusEvidence{s} = Sref.M.Predictions(s).rho; % activation

    % compute sensory kernel for each model
    for sp=1:nSpan
        for lt = 1:nLapseType
            for bt = 1:nBehType

                % load model file
                S = load(fullfile(ModelfitsDir,ModelFiles{sp,lt,bt}));

                % probability of response on each trial,according to model
                if isfield(S.S_all(s),'Y')
                    Ymodel = S.S_all(s).Y;
                else
                    Ymodel = S.S_all(s).lh;
                end

                nTrial = length(Ymodel);
                resp_model = Ymodel > rand(nTrial,1); % generate response from model

                % model accuracy
                AccuracyModel(sp,lt,bt,s) = mean(resp_model == target);
                AccuracyModelSem(sp,lt,bt,s) = sqrt( AccuracyModel(sp,lt,bt,s).* (1-AccuracyModel(sp,lt,bt,s)) / nTrial); % s.e.m.

                % model sensory kernel
                [PK_Model(:,sp,lt,bt,s),~,LRstats] = glmfit(stimulus,resp_model, 'binomial'); % logistic regression
                PK_ModelSe(:,sp,lt,bt,s) = LRstats.se;

                % AIC
                AIC(sp,lt,bt,s) = S.S_all(s).AIC;

                % AIC w.r.t integration model
                DeltaAIC(:,:,:,s) = AIC(:,:,:,s) - Sref.M.score.AIC(s);

                % Correlation between integration maps
                CorrelationIntegrationMap(sp,lt,bt,s) = S.S_all(s).IntegrationMap.r;
            end
        end
    end
end


%% SUPP FIGURE 2: ACCURACY, MODEL COMPARISON FOR ALL VARIANTS - PSYCHOMETRIC CURVE AND KERNEL FOR BEST VARIANT
figure;
hold on;

for s=1:nSubject

    %% Panel A: Accuracy
    hh(1,s) = subplot(5,nSubject,s); hold on; title(S.subject_label{s});
    this_AccuracyModel = reshape(AccuracyModel(:,:,:,s),nSpan,nType); % reshape as matrix (rows: span; columns: model type)
    this_AccuracyModelSem = reshape(AccuracyModelSem(:,:,:,s),nSpan,nType);

    % plot
    [~,~,h_acc] = wu(this_AccuracyModel,this_AccuracyModelSem, {1:nSpan AllTypeLabel}, 'color',ColorMod, 'linestyle',LineStyle);

    plot(xlim, Accuracy(s)*[1 1],'k');
    xlabel('span');
    if s==1
        ylabel('accuracy');
        legend(h_acc.mean([1 3]), {'non-probabilistic','probabilistic'});
    else
        legend off;
    end

    %% Panel B: AIC
    hh(2,s) = subplot(5, nSubject,s+2); hold on;
    this_DeltaAIC = reshape(DeltaAIC(:,:,:,s),nSpan,nType);

    % plot
    wu(this_DeltaAIC,[], 'color',ColorMod, 'linestyle',LineStyle);
    xlabel('span');
    if s==1
        ylabel('\Delta AIC');
    end
    set(gca, 'yscale','log');
    yticks([1e2 1e4]);
    ylim([10 1e5]);
    legend off;

    %% Use best model (with smallest AIC: Span 3, free lapses and probabilistic mapping) for psychometric curve and sensory kernel
    spn = 3; % span =3
    lt = 1; % free lapses
    bt = 2; % with sensory noise (probabilistic)

    %% Panel C: Psychometric curve
    hh(3,s) =subplot(5, nSubject,s+4); hold on;

    % load model fit and simulate it
    S = load(fullfile(ModelfitsDir,ModelFiles{sp,lt,bt}));
    if isfield(S.S_all(s),'Y')
        Ymodel = S.S_all(s).Y;
    else
        Ymodel = S.S_all(s).lh;
    end
    % generate response from model
    resp_model = Ymodel > rand(length(Ymodel),1);

    % model psychometric curve
    grpmeann(resp_model, StimulusEvidence{s}, 'quantiles',{nQuantile}, 'plot','ste','errorstyle','fill','color','b');

    % animal psychometric curve
    mask = T.subject== subject_id(s);
    resp = T.resp(mask);
    grpmeann(resp, StimulusEvidence{s}, 'quantiles',{nQuantile}, 'plot','ste','linestyle','none','marker','.', 'errorstyle','bar', 'color','k');
    axis tight;
    ylim([0 1]);
    if s==1, ylabel('p (right)'); end
    if strcmp(animal, 'monkey'), xlim([-5 5]); end
    xlabel('evidence');


    %% Panel D: Sensory kernel
    hh(4,s) =subplot(5, nSubject,s+6); hold on;

    % model kernel (remove intercept)
    wu(PK_Model(2:end,spn,lt,bt,s),PK_ModelSe(2:end,sp,lt,bt,s),'errorstyle','fill','color','b');

    %animal kernel
    wu(PK_Data(2:end,s),PK_DataSe(2:end,s),'linestyle','none','marker','.', 'errorstyle','bar', 'color','k');
    if s==1, ylabel('pulse weight'); end
    xlabel('pulse');
    axis tight;

    %% Panel E: correlation in 2d maps
    hh(5,s) =subplot(5, nSubject,s+8); hold on;
    this_CorrelationIntegrationMap = reshape(CorrelationIntegrationMap(:,:,:,s),nSpan,nType);

    % plot
    wu(this_CorrelationIntegrationMap,[], 'color',ColorMod, 'linestyle',LineStyle);
    plot(xlim, corr2dmap_ref(s)*[1 1],'k');
    xlabel('span');
    if s==1
        text(1.5,  corr2dmap_ref(s), 'integration model','fontangle','italic', 'verticalalignment','top','horizontalalignment','center');
        ylabel({'Correlation','2D map'});
    end
    legend off;
    xlim([0 6]);

end

% add panel letters
panelletter(hh(:,1));
