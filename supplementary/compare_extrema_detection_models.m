%% SUPPLEMENTARY FIGURE 3: COMPARING VARIANTS OF THE EXTREMA DETECTION MODEL
%(deterministic vs smooth, fixed vs free lapse, span length)

clear; close all;
animal = 'monkey';

nLapseType = 4;2; % free or fixed
nDefaultRule = 2; % default rules (if threshold is not reached at the end of sequence): random, or based on last sample

nType = nLapseType*nDefaultRule; % number of variant types ( 2x2)
LapseTypeLabel = {'free','fixed','free','fixed'};
%LapseTypeLabel = {{'free lapse,','fixed threshold'},{'fixed lapse,','fixed threshold'},{'free lapse,','varying threshold'},{'fixed lapse,','varying threshold'}};

%LapseTypeLabel = {'free&fixed','fixed&fixed','free&varying','fixed&varying'};

DefaultRuleLabel = {'random','last sample'};

% number of quantiles for psychometric curves
nQuantile = 50;

% directory where model fits are stored
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');


%% recover model files
models = cell(nLapseType, nDefaultRule);
for lt = 1:nLapseType
    for r = 1:nDefaultRule
        mdl = 'extremadetection';

        if any(lt ==[1 3]) % fixed lapse
            mdl = [mdl '_freelapse'];
        end
        if r ==2 % last sample
            mdl = [mdl '_lastsample'];
        end
        if lt>2
             mdl = [mdl '_varyingthreshold'];
        end   
        models{lt,r} = mdl;
    end
end

nModels = numel(models);
ModelColor = [0 .8 0];  % green

% model files
modelfiles = cell(size(models));
for m=1:nModels
    modelfiles{m}  = sprintf('%s/%s_%s',models{m}, animal,models{m});
end

% load integration model (reference model)
integration_model = 'integration';
IntegrationModelFile  = sprintf('%s/%s_%s',integration_model, animal,integration_model);
Sref = load(fullfile(ModelfitsDir,IntegrationModelFile));
subject_label = Sref.subject_label;

% integration map correlation for integration model
corr2dmap_ref = [Sref.M.score.IntegrationMap.r];

%% Load and preprocess data
csvfile = fullfile(TemporalIntegrationRootDirectory, 'data', animal);
T = readtable([csvfile '.csv']);
fprintf('loading behavioral data from %s\n',csvfile);

% preprocess data
[T, nSample, subject_id, nSubject, nSamples_max, nSamples_mean] = preprocess_table(T);


%% COMPUTE ACCURACY, PSYCHOMETRIC CURVES AND KERNELS FOR EXPERIMENTAL AND SIMULATED DATA, AIC AND CORRELATION BETWEEN INTEGRATION MAPS
Accuracy = zeros(1,nSubject);
AccuracyModel = zeros(nLapseType,nDefaultRule,nSubject);
AccuracyModelSem = zeros(nLapseType,nDefaultRule,nSubject);
StimulusEvidence = cell(1,nSubject);
PK_Data = zeros(nSample+1,nSubject);
PK_DataSe = zeros(nSample+1,nSubject);
PK_Model = zeros(nSample+1,nLapseType,nDefaultRule,nSubject);
PK_ModelSe = zeros(nSample+1,nLapseType,nDefaultRule,nSubject);
AIC = zeros(nLapseType, nDefaultRule,nSubject);
DeltaAIC = zeros(nLapseType, nDefaultRule,nSubject);
CorrelationIntegrationMap  = zeros(nLapseType, nDefaultRule,nSubject);
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

    % stimulus evidence on each trial(weighted evidence according to logistic
    % regression)
    StimulusEvidence{s} = Sref.M.Predictions(s).rho;

    % compute sensory kernel for each model
    for lt = 1:nLapseType
        for r = 1:nDefaultRule

            % load model file
            S = load(fullfile(ModelfitsDir,modelfiles{lt,r}));

            % probability of response on each trial,according to model
           % if isfield(S.S_all(s),'Y')
           %     Ymodel = S.S_all(s).Y;
           % else
                Ymodel = S.S_all(s).pModel;
           % end

            resp_model = Ymodel > rand(length(Ymodel),1); % generate response from model
            AccuracyModel(lt,r,s) = mean(resp_model == target); % model accuracy
            AccuracyModelSem(lt,r,s) = sqrt( AccuracyModel(s) .* (1-AccuracyModel(s)) / nTrial);

            % model sensory kernel
            [PK_Model(:,lt,r,s),~,LRstats] = glmfit(stimulus,resp_model, 'binomial'); % logistic regression
            PK_ModelSe(:,lt,r,s) = LRstats.se;

            % AIC
            AIC(lt,r,s) = S.S_all(s).AIC;

            % correlation between experimental and model integration maps
            CorrelationIntegrationMap(lt,r,s) = S.S_all(s).IntegrationMap.r;
        end
    end

    % AIC w.r.t integration model
    DeltaAIC(:,:,s) = AIC(:,:,s) - Sref.M.score.AIC(s);


end

%% MAKE FIGURE
figure;
for s=1:nSubject

    %% Accuracy panel
    sp(s,1) = subplot(5,nSubject,s); hold on; title(subject_label{s});
    wu(AccuracyModel(:,:,s),AccuracyModelSem(:,:,s), {LapseTypeLabel,DefaultRuleLabel});


    plot(xlim, Accuracy(s)*[1 1],'k');
    if s==1, ylabel('accuracy'); end
    if s==nSubject
        legend off;
    end

VerticalOffset = 0.02;
text(1.5, min(ylim) - VerticalOffset, 'fixed threshold','HorizontalAlignment','center','FontSize',8);
text(3.5, min(ylim) - VerticalOffset, 'varying threshold','HorizontalAlignment','center','FontSize',8);


    %% AIC panel
    sp(s,2) = subplot(5, nSubject,s+2); hold on;
    wu(DeltaAIC(:,:,s),[],{LapseTypeLabel});
    if s==1, ylabel('\Delta AIC'); end
    set(gca, 'yscale','log');
    ylim([1e2 5e3]);
    VerticalOffset = 4;
text(1.5, min(ylim) / VerticalOffset, 'fixed threshold','HorizontalAlignment','center','FontSize',8);
text(3.5, min(ylim) / VerticalOffset, 'varying threshold','HorizontalAlignment','center','FontSize',8);

    legend off;

    %% For sensory kernel and psychometric curve, use the best variant: free lapse and last sample rule
    lt = 3; % free lapses, varying threshold
    r = 2; % default rule: last sample

    %% Sensory kernel
    sp(s,3) = subplot(5, nSubject,s+4); hold on;

    % model kernel (quit bias)
    wu(PK_Model(2:end,lt,r,s),PK_ModelSe(2:end,lt,r,s),'errorstyle','fill','color',ModelColor);

    %animal kernel
    wu(PK_Data(2:end,s),PK_DataSe(2:end,s),'linestyle','none','marker','.', 'errorstyle','bar', 'color','k');
    if s==1, ylabel('sample weight'); end
    xlabel('sample');
    axis tight;


    %% Psychometric curve
    sp(s,4) = subplot(5, nSubject,s+6); hold on;
    S = load(fullfile(ModelfitsDir,modelfiles{lt,r}));

    % probability of response in each trial according to model
   % if isfield(S.S_all(s),'Y')
   %     Ymodel = S.S_all(s).Y;
   % else
        Ymodel = S.S_all(s).pModel;
   % end

    % generate response from model
    resp_model = Ymodel > rand(length(Ymodel),1);

    % model psychometric curve
    grpmeann(resp_model, StimulusEvidence{s}, 'quantiles',{nQuantile}, 'plot','ste','errorstyle','fill','color',ModelColor);

    % animal psychometric curve
    mask = T.subject== subject_id(s);
    resp = T.resp(mask);
    grpmeann(resp, StimulusEvidence{s}, 'quantiles',{nQuantile}, 'plot','ste','linestyle','none','marker','.', 'errorstyle','bar', 'color','k');

    axis tight;
    ylim([0 1]);
    if s==1, ylabel('p (right)'); end
    if strcmp(animal, 'monkey'), xlim([-5 5]); end
    xlabel('evidence');

    %% Integration map correlation
    sp(s,5) = subplot(5, nSubject,s+8); hold on;
    wu(CorrelationIntegrationMap(:,:,s),[],{LapseTypeLabel});
    plot(xlim, corr2dmap_ref(s)*[1 1],'k');
    text(3,  corr2dmap_ref(s), 'integration model','fontangle','italic', 'verticalalignment','top','horizontalalignment','center');
    legend off;
    if s==1, ylabel({'Correlation','2D map'}); end

    VerticalOffset = 0.05;
text(1.5, min(ylim) - VerticalOffset, 'fixed threshold','HorizontalAlignment','center','FontSize',8);
text(3.5, min(ylim) - VerticalOffset, 'varying threshold','HorizontalAlignment','center','FontSize',8);

end

% add panel letters
panelletter(sp(1,:));

