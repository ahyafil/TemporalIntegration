%% illustration of information integration analysis
clear; close all;

animal = 'monkey';

IntegrationMapColor = [241 135 34; 103 169 221]/255; % integration maps: gradient between orange (left) light blue (right)
PsychometricColor = [0 0 0; 1 0 0]; % color of conditional psychometric curves (graded from black to red)

% levels for contour lines
isoLevels = [.15 .3 .5 .7 .85];

EvidenceLabels = {'early evidence', 'late evidence'};

% values of late evidence we condition on for psychometric curves
LateEvidenceValue = -2:2;

%% LOAD MODEL FITS
models = {'integration', ... % integration model
    'snapshot1',... % snapshot model
    'extremadetection'}; % extrema detection
model_label = {'integration','snapshot',{'extrema', 'detection'}};

% number of models
nModel = length(models);

% directory where model fits are stored
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

% load
i = 2;
anymodel_file =  fullfile(ModelfitsDir,models{i},sprintf('%s_%s',animal, models{i}));

S = load(anymodel_file);
LateEvidenceBinEdges = S.S_all(1).IntegrationMap.LateEvidenceValue; % corresponding bins of late evidence
LateEvidenceBinCenters = (LateEvidenceBinEdges(1:end-1)+LateEvidenceBinEdges(2:end))/2;

%% load biases, lapses for simulated data
s = 1; % Monkey X


%% MAKE FIGURE
figure;
nCol = nModel+1; % number of panel columns

biases = cell(1,nModel);
biases_se = cell(1,nModel);
lapses = cell(1,nModel);
lapses_se = cell(1,nModel);
for m=1:nModel
    % load model file
    model_file =  fullfile(ModelfitsDir,models{m},sprintf('%s_%s',animal, models{m}));
    Stmp = load(model_file);

    % retrieve integration map structure
    if isfield(Stmp,'M') % gum
        IM = Stmp.M.score.IntegrationMap(s);
    else
        IM = Stmp.S_all(s).IntegrationMap;
    end

    % plot integration map
    subplot(4, nCol, m);
    IM.IntegrationMapModel = plot_integration_map(IM.IntegrationMapModel, IntegrationMapColor, isoLevels, EvidenceLabels);
    if m>1
        ylabel('');    
        yticklabels({});
    end
    ttl = title(model_label{m});
    ttl.Position = ttl.Position + [0 .6 0];

    %axes(h_data{2});

    % plot conditional psychometric curve
    subplot(4, nCol, m+nCol);
    DataPointsCutoff = 20;
    plot_conditional_psychometric_curve(IM.IntegrationMapModel, LateEvidenceValue, DataPointsCutoff, PsychometricColor);

    title(''); legend off;
    if m>1
        ylabel('');
    end
    ylim([0 1]);

    % fits of conditional psychometric curves
    PCfit = IM.ConditionalPsychometricCurvesModel; % parameter fits for psychometric curves conditioned on late evidence
    biases{m} = PCfit(1,:,1); % biases fit to monkey 1 data
    biases_se{m} = PCfit(1,:,2); % biases fit to monkey 1 data (standard error)
    lapses{m} = PCfit(3:4,:,1); % lapses fit to monkey 1 data
    lapses_se{m} = PCfit(3:4,:,2); % lapses fit to monkey 1 data

    % plot biases
    sb_beta(m) = subplot(4,nCol,2*nCol+m);
    wu(LateEvidenceBinCenters, biases{m},biases_se{m}, 'linewidth',2, 'color','k', 'errorstyle','fill');
    xticklabels({}); box off;
    if  m==1, ylabel('bias \beta'); end

    % plot lapses
    sb_pi(m) = subplot(4,nCol,3*nCol+m);
    wu(LateEvidenceBinCenters, lapses{m}',lapses_se{m}', {{},{'left','right'}}, 'linewidth',2, 'color',IntegrationMapColor, 'errorstyle','fill');
    ylim([-.2 0.7]);
    xlabel('late evidence');
    if  m==1
        ylabel('lapse \pi');
    else
        legend off;
    end
    box off;

end

%% add colorbar
subplot(4,nCol,nCol);
cmap = IntegrationMapColor(1,:)+ linspace(0,1,64)'.*diff(IntegrationMapColor,1); % colormap

colormap(cmap);
cb = colorbar;
title(cb, 'p(right)')
cb.Position(1) = .75;
axis off;

% adjust limits of axes
sameaxis(sb_beta);
sameaxis(sb_pi);
subplot_shift(sb_pi,0,.04);
