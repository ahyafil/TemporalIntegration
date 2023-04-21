%% Plot Figure 2, 5 or 6
% Displays fits of psychometric curves, sensory kernels, model comparison for integrant and no-integration models

clear; close all;
animal = 'rat'; %'rat','monkey', 'human'

%% % %%
models = {'integration', ... % integration model
    'snapshot1',... % snapshot model with span=1, fixed lapse and deterministic mapping
    'extremadetection'}; % extrema detection with fixed lapses
model_label = {'integration','snapshot','extrema detection'}; % model labels
nModel = length(models);
ModelColor = {[1 0 0];[0 0 .8]; [0 .5 0]}; % red for integration, blue for snapshot, green for extrema detection

% number of quantiles of stimulus evidence for psychometric curves
if ~strcmp(animal, 'human')
    nQuantile = 50;
else
    nQuantile = 10;
end

%% get path for figures and data
ModelSaveDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');


%% 1. LOAD AND PREPROCESS DATA
csvfile = fullfile(TemporalIntegrationRootDirectory, 'data', animal);
T = readtable([csvfile '.csv']);
fprintf('loading behavioral data from %s.csv\n',csvfile);

% preprocess data
[T, nSample, subject_id, nSubject, nSamples_max, nSamples_mean] = preprocess_table(T);

% define subject label
switch animal
    case 'monkey'
        subject_label = "Monkey "+subject_id;
    case 'rat'
        subject_label = subject_id;
    case 'human'
        % select only valid trials
        T = T(T.valid==1,:);
        subject_label = "S" + subject_id;
    otherwise
        error('don''t have data to analyze for this animal');
end

% model files
ModelFiles = cell(1,nModel);
for m=1:nModel
    ModelFiles{m}  = sprintf('%s/%s_%s',models{m}, animal,models{m});
end

% file with integration model, to compute weighted early and late evidence,
% and stimulus evidence
IntegrationModelFile =  ModelFiles{1};
integration_file = load(fullfile(ModelSaveDir,IntegrationModelFile));


%% 2. COMPUTE PSYCHOMETRIC CURVES, KERNELS AND ACCURACY

Accuracy = zeros(1,nSubject);
AccuracySem = zeros(1,nSubject);
if strcmp(animal, 'rat') % for rat we only plot PK for trials with 10 samples (in main figure)
    nSamplePK = 10;
else
    nSamplePK = nSample;
end
PK = nan(1+nSamplePK,nSubject);
PK_se = nan(1+nSamplePK,nSubject);
PK20 = zeros(1+nSample,nSubject);
PK20_se = zeros(1+nSample,nSubject);
StimulusEvidence = cell(1,nSubject);
RespModel = cell(nModel, nSubject);
AccuracyModel = zeros(nModel, nSubject);
AccuracyModelSem = zeros(nModel, nSubject);
PK_Model = zeros(1+nSamplePK,nModel, nSubject);
PK_ModelSem = zeros(1+nSamplePK,nModel, nSubject);
PK20_Model = zeros(1+nSample,nModel, nSubject);
PK20_ModelSem = zeros(1+nSample,nModel, nSubject);


for s=1:nSubject
    % retrieve weighted evidence from integration model
    StimulusEvidence{s} = integration_file.M.Predictions(s).rho;

    nS = nSamples_max(s); %max number of samples for this subject

    mask = T.subject==subject_id(s); % trials for this subject
    nTrial = sum(mask); % number of trials for this subject
    stimulus = T.stimulus(mask,1:nS); % stimulus matrix for this subject
    target = T.target(mask); % target response for this subject
    resp = T.resp(mask); % response for this subject

    % 2.1: mean accuracy of subject
    Accuracy(s) = mean(resp == target);

    % 2.2 psychophysical kernels for participant
    if ~strcmp(animal, 'rat')
        % logistic regression
        [PK(1:nS+1,s),~,LRstats] = glmfit(stimulus(:,1:nS),resp, 'binomial');
        PK_se(1:nS+1,s) =  LRstats.se; % standard error for weights
    else
        % for rats we only take trials with 10 samples (500 ms)
        nS_rat = 10;
        this_trial = sum(~isnan(stimulus),2)' == nS_rat; % select trials matching number of samples

        % logistic regression
        if any(this_trial)
            [PK(1:nS_rat+1,s),~,LRstats] = glmfit(stimulus(this_trial,1:nS_rat),resp(this_trial), 'binomial');
            PK_se(1:nS_rat+1,s) = LRstats.se;
        end

        % compute kernel for 20-sample stimuli aside
        if nS==20
            this_trial = sum(~isnan(stimulus),2)' == nS;  % select trials matching number of samples

            [PK20(1:nS+1,s),~,LRstats] = glmfit(stimulus(this_trial,1:nS),resp(this_trial), 'binomial'); % logistic regression
            PK20_se(1:nS+1,s) = LRstats.se;
        end
    end

    %% 2.3 psychophysical kernel for each model
    for m=1:nModel
        % load model fit
        S = load(fullfile(ModelSaveDir,ModelFiles{m}));

        % p(right) on each trial according to model
        if isfield(S, 'M') % GUM (integration model)
            pModel = S.M.Predictions(s).Expected;
        elseif isfield(S.S_all(s),'Y')
            pModel = S.S_all(s).Y;
        else
            pModel = S.S_all(s).pModel;
        end

        % draw response from model
        RespModel{m,s} = pModel > rand(nTrial,1);

        % model accuracy
        AccuracyModel(m,s) = mean(RespModel{m,s} == target);
        AccuracyModelSem(m,s) = sqrt( AccuracyModel(s) .* (1-AccuracyModel(s)) / nTrial); % standard error of the mean

        % model psychophysical kernel
        if ~strcmp(animal, 'rat')
            [PK_Model(1:nS+1,m,s),~,LRstats] = glmfit(stimulus(:,1:nS),RespModel{m,s}, 'binomial'); % logistic regression
            PK_ModelSem(1:nS+1,m,s) = LRstats.se;

        else % compute kernel separately for 10 & 20 sample stim
            nS_rat = 10;
            this_trial = sum(~isnan(stimulus),2)' == nS_rat;  % select trials matching number of samples

            if any(this_trial)
                [PK_Model(1:nS_rat+1,m,s),~,LRstats] = glmfit(stimulus(this_trial,1:nS_rat),RespModel{m,s}(this_trial), 'binomial'); % logistic regression
                PK_ModelSem(1:nS_rat+1,m,s) = LRstats.se;
            end

            % compute kernel for 20-sample stimuli aside
            if nS==20
                this_trial = sum(~isnan(stimulus),2)' == nS;  % select trials matching number of samples

                [PK20_Model(1:nS+1,m,s),~,LRstats] = glmfit(stimulus(this_trial,1:nS),RespModel{m,s}(this_trial), 'binomial'); % logistic regression
                PK20_ModelSem(1:nS+1,m,s) = LRstats.se;
            end
        end
    end
end

%% 3. COMPUTE BOUNDARY PERFORMANCE FOR SNAPSHOT AND EXTREMA DETECTION MODELS
maxperf_snapshot_sample = zeros(nSubject, nSample);
maxperf_snapshot = zeros(1,nSubject);
opt_sample = zeros(1,nSubject);
maxperf_extrema  = zeros(1,nSubject);
opt_threshold  = zeros(1,nSubject);
maxperf_extrema_threshold = cell(1,nSubject);

for s=1:nSubject

    mask = T.subject==subject_id(s); % trials for this subject
    stimulus = T.stimulus(mask,:); % stimulus matrix for this subject
    target = T.target(mask); % target response for this subject

    %% compute boundary performance for snapshot model

    % for each trial, whether each frame is aligned with target response
    aligned_snapshot = (stimulus>0 & target ==1) +  (stimulus<0 & target ==0) + (stimulus==0)/2;

    % maximum performance overall per frame, if attending this frame only
    maxperf_snapshot_sample(s,:) = mean(aligned_snapshot);

    % boundary performance corresponds to max over samples (give optimal
    % sample as well)
    [maxperf_snapshot(s), opt_sample(s)] = max(maxperf_snapshot_sample(s,:));

    %% compute boundary performance for extrema detection model

    % list of all possible thresholds
    if strcmpi(animal,'monkey')
        Threshold = 1:max(abs(stimulus(:))); % all possible threshold (integers)
    else
        Threshold = 0:.01:1; % continuous values
    end

    nThr = length(Threshold);

    % response on each trial for given threshold value
    PP = zeros(size(stimulus,1), nThr);
    for f=1:nSamples_max(s) % for each sample
        still = ~PP; % trials/models without a response yet
        reach = abs(stimulus(:,f)) >= Threshold; % whether reaches threshold value (trial x threshold boolean matrix)
        newReached = reach & still; % new responses: when reaches threshold for first time
        newResponse = repmat(sign(stimulus(:,f)), 1, nThr); % response depends on the sign for that sample (-1/+1)
        PP(newReached) = newResponse(newReached); % define responses for trial reaching response
    end

    % for each trial, whether extrema detection model is aligned with target
    % response (for each threshold level)
    aligned_extrema = (PP>0 & target==1) +  (PP<0 & target==0) + (PP==0)/2;

    % average over trials, given the overall performance for all threshold
    maxperf_extrema_threshold{s} = mean(aligned_extrema);

    % take max over possible thresholds
    [maxperf_extrema(s), opt_threshold(s)] = max(maxperf_extrema_threshold{s});

end

% maximum perforance for all models (exclude integration because it would
% be 1)
maxperf_allmodels = [nan(1,nSubject); maxperf_snapshot; maxperf_extrema]; % model x subject matrix

% number of subjects for which boundary performance for extrema detection
% and snapshot model are below subject accuracy
nSubjectBoundSnapshotBelow = sum(maxperf_snapshot<Accuracy);
nSubjectBoundExtremaBelow = sum(maxperf_extrema<Accuracy);

fprintf('In %d/%d %s subjects, boundary performance for extrema detection model is below subject accuracy\n',nSubjectBoundExtremaBelow,nSubject, animal);
fprintf('In %d/%d %s subjects, boundary performance for snapshot model is below subject accuracy',nSubjectBoundSnapshotBelow,nSubject, animal);


%% 4. MODEL COMPARISON
AIC = zeros(nModel,nSubject);

for m=1:nModel
    % load model fit
    S = load(fullfile(ModelSaveDir,ModelFiles{m}));

    if isfield(S, 'M') % gum
        AIC(m,:) = S.M.score.AIC;
    else
        for s=1:nSubject
            AIC(m,s) = S.S_all(s).AIC;
        end

    end
end

% Difference with AIC of integration model
DiffAIC = AIC-AIC(1,:);


%% 5. INTEGRATION MAPS: CORRELATION BETWEEN DATA AND MODEL

CorrelationIntegrationMap = zeros(nSubject,length(models)); % subject x model matrix
for m=1:length(models) % loop through models
    % load integration file
    S = load(fullfile(ModelSaveDir,ModelFiles{m}));

    for s = 1:nSubject
        if isfield(S, 'M')
            IM = S.M.score.IntegrationMap(s);
            IntegrationMap(:,:,s) = IM.IntegrationMapData.IntegrationMap;
            IM_bnd = IM.bnd;
            IM_dx = IM.dx;
        else
            IM = S.S_all(s).IntegrationMap;
        end
        CorrelationIntegrationMap(s,m) = IM.r; % correlation score
        CorrelationIntegrationMapBoostrap(:,s,m) = IM.rBoostrap;% correlation of boostraps
    end
end
IntegrationMapAverage = mean(IntegrationMap,3); % average over subjects
IMAverage = IM.IntegrationMapData;
IMAverage.IntegrationMap = IntegrationMapAverage;

pIntegrationMap = nan(nSubject,2);
for s=1:nSubject
    [~,pIntegrationMap(s,1)] = ttest2(CorrelationIntegrationMapBoostrap(:,s,1),CorrelationIntegrationMapBoostrap(:,s,2));
    [~,pIntegrationMap(s,2)] = ttest2(CorrelationIntegrationMapBoostrap(:,s,1),CorrelationIntegrationMapBoostrap(:,s,3));
    fprintf('ttest integration map integration vs snapshot models (boostrap, %s):%f\n',subject_label{s},pIntegrationMap(s,1));
    fprintf('ttest integration map integration vs extrema detection models (boostrap, %s):%f\n',subject_label{s},pIntegrationMap(s,2));
end
fprintf('ttest integration map integration vs snapshot models (boostraps):%d/%d significant at p<0.001\n',sum(pIntegrationMap(:,1)<0.001),nSubject);
fprintf('ttest integration map integration vs extrema detection models  (boostraps):%d/%d significant at p<0.001\n',sum(pIntegrationMap(:,2)<0.001),nSubject);


%% 6. ANALYSIS OF DISAGREE TRIALS (for monkeys this is done on a separate script)
if ~strcmp(animal, 'monkey')
    disagree_file = sprintf('disagree/%s_disagree.mat',animal);
    i_model = [1 3]; % for integration and extrema detection models

    mlabel = [{'subjects'}, model_label(i_model)];
    pTotalEvidenceChoiceOnDisagree = zeros(nSubject,3);
    for m=1:3
        % load analysis of disagree files
        if m==1 % for experimental data
            S = load(fullfile(ModelSaveDir,disagree_file));
        else % for simulations
            S = load(fullfile(ModelSaveDir,models{i_model(m-1)}, disagree_file));
        end

        % proportion of choice agligned with total evidence on disagree
        % trials (for each subject)
        pTotalEvidenceChoiceOnDisagree(:,m) = S.pTotalEvidenceChoiceOnDisagree;
    end
end

%% COMPOSE MAIN FIGURE (FIGURE 2, 4 or 5)
if strcmp(animal, 'monkey')
    %% FIGURE 2 - MONKEYS

    figure;
    set(gcf, 'name', 'Figure 2');
    nCol = 4; % number of subplot columns

    %% psychometric curves and sensory kernel for each model
    for m=1:nModel

        for s=1:2 % for each animal

            ss = 3-s; % reverse order (monkey N above monkey P)

            %% Panel B: sensory kernel
            sp_kern(m,s)  =  subplot(4,nCol,m+1+nCol*(s-1)); hold on;
            subplot_shift(gca, -(m+1)*.02,0);

            [~,~,h_mod] = wu(PK_Model(2:end,m,ss),PK_ModelSem(2:end,m,ss),...
                'errorstyle','fill','color',ModelColor{m}); % model kernel (without intercept)
            [~,~,h_data] = wu(PK(2:end,ss),PK_ModelSem(2:end,ss),'linestyle','none',...
                'marker','.', 'errorstyle','bar', 'color','k'); %animal kernel
            if m==1
                ylabel(subject_label{ss});
            else
                yticklabels({});
            end
            if m==1 && s==1
                legend([h_mod.mean h_data.mean], {'model','data'}, 'location','northeast');
                legend boxoff
            end
            if s==1
                xticklabels([]);
                title(model_label{m});
            else
                xlabel('sample position');
            end
            axis tight;

            %% Panel D: psychometric curves
            sp_pc(m,s) = subplot(4,nCol,(1+s)*nCol+m+1); hold on;
            subplot_shift(gca, -(m+1)*.02,0);

            % model psychometric curve (computes proportion of rightward
            % responses for each bin of stimulus evidence)
            grpmeann(RespModel{m,ss}, StimulusEvidence{ss}, 'quantiles',{nQuantile},...
                'plot','ste','errorstyle','fill','color',ModelColor{m});

            resp = T.resp(T.subject==subject_id(ss));
            grpmeann(resp, StimulusEvidence{ss}, 'quantiles',{nQuantile},...
                'plot','ste','linestyle','none','marker','.', 'errorstyle','bar', 'color','k'); % animal psychometric curve

            plot(xlim, [.5 .5], 'color',.7*[1 1 1]);
            plot([0 0], ylim,  'color',.7*[1 1 1]);
            xlim([-4 4]); ylim([0 1]);
            if strcmp(animal, 'monkey'), xlim([-5 5]); end
            if m==1
                ylabel(subject_label{ss}); %'p(rightward)');
            else
                yticklabels({});
            end
            if s==1
                xticklabels([]);
            else
                xlabel('stimulus evidence');
            end

        end
    end
    subplot_shift(sp_kern(:),.12,0);
    subplot_shift(sp_pc(:),.12,0);
    subplot_shift(sp_kern(:,2),0,.03);
    subplot_shift(sp_pc(:,2),0,.03);

    for s=1:2, sameaxis(sp_kern(:,s)); end

    %% Panel C: accuracy
    sp_ac = subplot(2,nCol,1+nCol); hold on;
    subplot_resize(gca, [1.3 1],'left');

    AccuracyAll = cat(1,Accuracy,AccuracyModel)'; % join model and data
    AccuracyAllSem = cat(1,AccuracySem,AccuracyModelSem)';

    % plot
    [~,~,h_acc] = wu(AccuracyAll,AccuracyAllSem, 'bar',{subject_label [{'animal'} model_label]});
    bndperf_color = {[],[.5 .5 1], [.5 1 .5] };

    % add bars for boundary performance for extrema detection and snaphost
    % model
    for m=2:3
        plot((1:2)+h_acc.mean(m+1).XOffset+.08*[-1;1], [1;1]*maxperf_allmodels(m,:), 'color',bndperf_color{m},'linewidth',2);  % boundary perf
    end

    ylim([.68 .92])
    ylabel('accuracy');
    legend off;

    %% Panel A: model comparison
    sp_mc = subplot(2,nCol,1); hold on;
    subplot_resize(gca, [1.3 1],'left');
    set(gca, 'yscale','log');
    eps = 1; % to plot a 0 on log scale

    % plot
    [~,~,h_aic] = wu(eps + DiffAIC',[],'bar','basevalue',.8*eps,{subject_label model_label},'color',ModelColor);

    plot(xlim,eps*[1 1], 'color',.7*[1 1 1]);
    set(gca, 'ytick',[eps 1e2 1e4],'yticklabel',{'0','100','10.000'});
    ylabel('AIC difference');
    ylim([.8*eps 2e4]);

    for m=1:length(model_label)
        htext = text(1+h_aic.mean(m).XOffset, 1.2*eps, model_label{m});
        set(htext, 'Rotation',90,'Color',[1 1 1]*(m>1)); % black for first, white for last two
    end

    legend off;

    %% add panel letters
    [~,ax] = panelletter([sp_mc sp_kern(1) sp_ac sp_pc(1)]);

    axes(ax);
    text(.36, .67, 'sample weight','rotation',90);
    text(.36, .24, 'p(rightward)','rotation',90);

else
    %% FIGURE FOR HUMANS AND RATS (FIGURE 5 & 6)
    figure;
    set(gcf, 'name',"Figure "+54+strcmp(animal, 'rat')));
    nRow = 3; % number of panel rows
    nCol = 4; % number of panel columns


    %% psychometric curves and sensory kernel for each model
    for m=1:nModel

        % Panel A: psychometric curve
        sp(1,m) = subplot(nRow,nCol,m); hold on;
        subplot_shift(gca, -(m)*.02,0);

        % bins for stimulus evidence
        if strcmp(animal,'human')
            StimulusEvidenceBins = -3:3;
        else
            StimulusEvidenceBins = -3:.1:3;
        end
        BinCenters = StimulusEvidenceBins(1:end-1) + diff(StimulusEvidenceBins)/2;

        % compute mean response for model and data in each bin of stimulus
        % evidenec
        for s=1:nSubject
            PC_model(s,:) = grpmeann(RespModel{m,s}, StimulusEvidence{s}, 'histc',{StimulusEvidenceBins}); % model psychometric curve

            mask = T.subject==subject_id(s); % trials for this subject
            PC_data(s,:) = grpmeann(T.resp(mask), StimulusEvidence{s}, 'histc',{StimulusEvidenceBins}); % animal psychometric curve
        end

        % plot
        wu(BinCenters, PC_model,'ste','errorstyle','fill','color',ModelColor{m}); % model psychometric curve
        wu(BinCenters, PC_data,'ste','linestyle','none','marker','.', 'errorstyle','bar', 'color','k'); % animal psychometric curve
        axis tight;
        ylim([0 1]);
        if m==1
            ylabel('p(rightward)');
        else
            yticklabels({});
        end
        set(gca, 'xtick',[-2 0 2]);
        title(model_label{m});
        xlabel('stimulus evidence');

        % Panel C: sensory kernel
        sp(3,m) =  subplot(nRow,nCol,m+4); hold on;
        subplot_shift(gca, -(m)*.02,0);

        wu(permute(PK_Model(2:end,m,:),[3 1 2]),'ste','errorstyle','fill','color',ModelColor{m}); % model kernel (quit bias)
        wu(PK(2:end,:)','ste','linestyle','none','marker','.', 'errorstyle','bar', 'color','k'); %animal kernel
        if m==1
            ylabel('sample weight');
        else
            yticklabels({});
        end
        title(model_label{m});
        xlabel('sample');
        axis tight;

    end
    sameaxis(sp(3,:));

    %% Panel B: Accuracy
    sp(2,1) = subplot(nRow,nCol,4); hold on;
    symbol = '*xov^sphd><.';
    dx = 1e-3*(-1:1); % shift just enough to avoid overlap of lines
    for s=1:nSubject
        % accuracy of fitted models
        for m=1:3
            plot(Accuracy(s) + dx(m), AccuracyModel(m,s), symbol(s),'color', ModelColor{m},'markersize',8);
        end
    end
    xlim(ylim);
    plot(xlim,ylim,'k');
    ylabel('model accuracy');
    xlabel('animal accuracy');


    %% Panel D: model comparison
    sp(4,1) =subplot(nRow,nCol,8); hold on;
    eps = strcmp(animal,'rat'); % to plot for log scale

    for s=1:nSubject
        % accuracy of fitted models
        for m=1:3
            plot(m, DiffAIC(m,s)+eps, symbol(s),'color', ModelColor{m},'markersize',8);
        end
    end
    ylabel('AIC difference');
    xlim([0.3 3.5]);
    xticks(1:3); xticklabels({'integration', 'snapshot','extrema'});
    set(gca, 'xticklabelrotation',45);
    if strcmp(animal,'rat')
        set(gca, 'yscale','log','ytick',[1 1e2 1e4], 'yticklabel',{'0','100','10^4'});

    end

    %% Panel E: analysis of disagree trials
    sp(5,1) =subplot(nRow,nCol-1,2*nCol-1); hold on;
    DisagreeColors = {.4*[1 1 1],ModelColor{1} ModelColor{2}}; % color for experimental data, integration model, extrema detection

    % plot as bar plot
    wu(pTotalEvidenceChoiceOnDisagree,  { mlabel}, 'bar','color',DisagreeColors,'basevalue',0.5);
    plot(xlim,.5*[1 1], 'color',.7*[1 1 1]);
    ylabel('p. majority-aligned choices');
    set(gca, 'xticklabelrotation',60);
    if strcmp(animal,'rat'), ylim([.5 .63]); end

    %% Panels F: early-late integration map
    sp(6,1) =subplot(nRow,nCol-1,2*nCol);
    IntegrationMapColour = [ 241 135 34; 103 169 221]/255; % gradient between red and blue


    isoLevels = [.15 .3 .5 .7 .85];
    plot_integration_map(IMAverage, IntegrationMapColour, isoLevels, {'early evidence','late evidence'});

    subplot_shift(gca, -.02);

    subplot(nRow,13,nRow*13-4);
    cmap = IntegrationMapColour(1,:)+ linspace(0,1,64)'.*diff(IntegrationMapColour,1); % colormap
    colormap(gca, cmap);
    cb = colorbar;
    title(cb, 'p(right)')
    %cb.Position(1) = .75;
    axis off;

    %%  Panel G: correlation of integration maps
    sp(7,1) =subplot(nRow,nCol,nRow*(nCol));
    model_label = {'integration','snapshot','extrema'};

    % bar plot
    [~,~,hbar] = wu(CorrelationIntegrationMap, 'bar','color',ModelColor);
    ylabel({'model-data','map correlation'});
    legend off;


    if strcmp(animal, 'human')
        % add labels for animal
        y_text = .89+.05*strcmp(animal,'rat');
        for m=1:length(model_label)
            htext = text(m-.1, y_text, model_label{m});
            set(htext, 'Rotation',90,'Color',[1 1 1]*(m<3)); % white for first two, black for last
        end
    end
    if strcmp(animal,'rat')
        ylim([.94 1]);
    else
        ylim([.85 1]);
    end

    %%
    width = 16;
    height = 16;

    % Position plot on the screen for drawing
    set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

    % Position plot on the paper for printing
    set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
        'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

    panelletter(sp(:,1));

    %% SUPP FIGURE: MAXIMUM PERFORMANCE OF NON-INTEGRATION MODELS VS SUBJECT PERF
    subfigure(['Supp Fig - Boundary performance' animal]);
    symbol = '*xov^sphd><.';
    for m=2:3 % for both non-integration models
        subplot(1,2,m-1); hold on;
        for s=1:nSubject % plot one symbol per subject
            % accuracy of fitted models
            plot(Accuracy(s), maxperf_allmodels(m,s), symbol(s),'color', ModelColor{m},'markersize',8);
        end
        all_data = [Accuracy(:); maxperf_allmodels(:)];
        xylim = [min(all_data)-.01 max(all_data)+.01];
        xlim(xylim);
        ylim(xylim);

        plot(xlim,ylim,'k');
        if m==2, ylabel('maximum model accuracy'); end
        xlabel('subject accuracy');
        title(model_label{m});
    end

    %% SUPP FIGURE: KERNEL FOR 20-SAMPLE STIMULI IN RATS
    if strcmp(animal, 'rat')
        subfigure('kernel long sequence rats');
        for m=1:3 % for each model
            subplot(1,3,m);
            % plot average over animals
            wu(permute(PK20_Model(2:end,m,:),[3 1 2]),'ste','errorstyle','fill','color',ModelColor{m}); % model kernel (no intercept)
            wu(PK20(2:end,:)','ste','linestyle','none','marker','.', 'errorstyle','bar', 'color','k'); %animal kernel
            title(model_label{m});
            xlabel('sample');
            if m==1, ylabel('sample weight'); end
        end
        sameaxis;
    end

end
