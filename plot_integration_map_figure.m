%% Plot integration maps and related analysis in monkey choices (Figure 4 and Supp Figure 5)
clear; close all;

models = {'integration', ... % logistic regression
    'snapshot1',... % snapshot model
    'extremadetection'}; % extrema detection

animal = 'monkey';

IntegrationMapColor = [241 135 34; 103 169 221]/255; % integration maps: gradient between orange (left) light blue (right)
PsychometricColor = [0 0 0; 1 0 0]; % color of conditional psychometric curves (graded from black to red)
ModelColor = {[1 0 0];[0 0 .8]; [0 .5 0]}; % red for integration, blue for snapshot, green for extrema detection

% max value of evidence (define boundaries of integration map)
MaxEvidence = 2.5;

% levels for contour lines
isoLevels = [.15 .3 .5 .7 .85];

% values of late evidence we condition on for psychometric curves
LateEvidenceValue = -2:2;

%% 1. DEFINE VANILLA INTEGRATION / NO-INTEGRATION MODEl
SensitivityIntegration = 1; % sensitivity to stimulus in integration model
SensitivityNoIntegration = 1.5; % sensitivity to stimulus in no integration model
toy_model_label = {'integration','no integration'};

% logistic function (mapping one-d stim to probability of response)
logistic = @(x) 1./(1+exp(-SensitivityIntegration*x));

% integration model: based on sum of two evidences
model_fun{1} = @(x,y) logistic(x+y);
model_eqn{1} = {'$p($right$) = \sigma(E_t+L_t)$'};

% no-integration model: based on either first or second evidence only, with
% equal probability (mixture model)
model_fun{2} = @(x,y) (logistic(SensitivityNoIntegration*x)+logistic(SensitivityNoIntegration*y))/2;
model_eqn{2} = {'$p($right$) = 0.5 \sigma(E_t) +0.5 \sigma(L_t)$'};

%% 2. RETRIEVE RESULTS OF INTEGRATION MAPS FOR EXPERIMENTAL DATA

% where to load analysis from from
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

% load any model, just to define bins of evidence
i = 2;
anymodel_file =  fullfile(ModelfitsDir,models{i},sprintf('%s_%s',animal, models{i}));

S = load(anymodel_file);
LateEvidenceBinEdges = S.S_all(1).IntegrationMap.LateEvidenceValue;  % corresponding bins of late evidence
LateEvidenceBinCenters = (LateEvidenceBinEdges(1:end-1)+LateEvidenceBinEdges(2:end))/2;

nAnimal = length(S.subject_label);
EvidenceLabels = {'early evidence', 'late evidence'};

% correlation rho between data and model 2d map
corr_2dmap = zeros(nAnimal,length(models));
for m=1:length(models)
    integration_file =  fullfile(ModelfitsDir,models{m},sprintf('%s_%s',animal, models{m}));
    Stmp = load(integration_file);
    for s = 1:nAnimal
        % retrieve integration map structure
        if isfield(Stmp,'M') % gum
            IM(s) = Stmp.M.score.IntegrationMap(s);
        else
            IM(s) = Stmp.S_all(s).IntegrationMap;
        end

        corr_2dmap(s,m) = IM(s).r;

    end
end


nPC = length(LateEvidenceValue); % number of psychometric curves

hh = -MaxEvidence:.01:MaxEvidence; % span of early/late evidence

% all possible combination of early evidence (X) and late evidence(Y) for 2D plot
[X,Y] = meshgrid(hh, hh);

% all possible combination of early evidence (X) and late evidence(Y) for psychometric curves
[x_pc,y_pc] = meshgrid(hh, LateEvidenceValue);

PsychometricCurveColor = PsychometricColor(1,:) + linspace(0,1,nPC)'*diff(PsychometricColor); % graded colours for psychometric curves
PsychometricCurveColormap = PsychometricColor(1,:) + linspace(0,1,32)'*diff(PsychometricColor); % graded colours for psychometric curves


%% create figure
figure;

for m=1:2 % for each model

    % compute probability of right choice for each value of (early evidence, late evidence)
    M = model_fun{m}(X,Y);

    %% 3. PANEL A: INTEGRATION MAP FOR MODELS
    sb(m,1) = subplot(4,4,m); hold on;

    % map probability of right onto RGB matrix
    C = zeros(length(hh), length(hh),3);
    for c=1:3
        C(:,:,c) = 1- (IntegrationMapColor(2,c)-IntegrationMapColor(1,c)) * M - IntegrationMapColor(1,c);
    end

    % plot integration map for vanilla model
    image(hh,hh,permute(1-C,[2 1 3])); axis xy; hold on;

    % add contour lines
    contour(hh, hh, M', isoLevels,'k', 'linewidth',.5); % draw contour lines
    contour(hh, hh, M',isoLevels(3)*[1 1],'k', 'LineWidth',1); % thicker contour lines

    % add arrows for values of late evidence in conditional psychometric
    % curves
    for p=1:nPC
        ha = arrow([MaxEvidence+.2 LateEvidenceValue(p)],  [MaxEvidence LateEvidenceValue(p)]);
        arrow(ha, 'length',8,'width',2);
        set(ha, 'facecolor',PsychometricCurveColor(p,:),'edgecolor',[1 1 1]);
    end

    % axis labels and limits, title
    xlim([-MaxEvidence MaxEvidence+.2]);
    xlabel('early evidence $E_t$', 'interpreter','latex');
    if m==1, ylabel('late evidence $L_t$', 'interpreter','latex');
    else
        yticklabels({});
    end
    title(toy_model_label{m});

    %% 4. DRAW PANEL B: CONDITIONAL PSYCHOMETRIC CURVES FOR MODELS
    sb(m,2) = subplot(4,4,m+4);

    % compute probability of right choice for each value of (early evidence, late evidence)
    PC = model_fun{m}(x_pc,y_pc);

    % draw conditional psychometric curves
    wu(hh, PC',[], 'color',PsychometricCurveColor, {{},num2strcell(LateEvidenceValue)});
    xlim([-MaxEvidence MaxEvidence+.2]);
    xlabel('early evidence $E_t$', 'interpreter','latex');
    if m==1 % integration model
        ylabel('p(right)');

        % add graded horizontal arrow (integration model)
        ha(1) = arrow([0 .5],  [1.5 .5]);
        ha(2) =  arrow([0 .5],  [-1.5 .5]);

        set(ha(1), 'FaceColor',PsychometricColor(1,:), 'edgecolor',PsychometricColor(1,:));
        set(ha(2), 'FaceColor',PsychometricColor(2,:), 'edgecolor',PsychometricColor(2,:));
        arrow(ha, 'length',8,'width',2,'tipangle',30);
        xx = linspace(-1.1,1.1,20);

        % add colormap
        colormap(gca, flipud(PsychometricCurveColormap));
        surface([xx;xx],.5*ones(2,20),zeros(2,20),[xx;xx],...
            'facecol','no','edgecol','interp', 'linew',3);
    else % no integration model
        %  if m==2
        % vertical graded arrow (no-integration model)
        ha(2) =   arrow([0 .5],  [0 .95]);
        ha(1) =   arrow([0 .5],  [0 .05]);
       % arrow(ha, 'length',8,'width',2);

        set(ha(1), 'FaceColor',PsychometricColor(1,:), 'edgecolor',PsychometricColor(1,:));
        set(ha(2), 'FaceColor',PsychometricColor(2,:), 'edgecolor',PsychometricColor(2,:));
        arrow(ha, 'length',8,'tipangle',30,'width',2);

        yy = linspace(.9,.1,20);

        colormap(gca, PsychometricCurveColormap);
        surface(zeros(2,20),[yy;yy], zeros(2,20),[yy;yy],...
            'facecol','no','edgecol','interp', 'linew',3);

        % end
        yticklabels({});
    end
    ylim([0 1]);
end

% legend of conditional psychometric curves
leg = findobj(gcf, 'type','legend');
set(leg(1), 'visible','off');
title(leg(2), 'late evidence');
leg(2).Position(1) = .68;

%% 5. PLOT INTEGRATION MAP AND CONDITIONAL PSYCHOMETRIC CURVES FOR EXPERIMENTAL DATA
s = 1; % select monkey N

% plot integration map
subplot(4,4,3);
IM(s).IntegrationMapData = plot_integration_map(IM(s).IntegrationMapData, IntegrationMapColor, isoLevels, EvidenceLabels);

% set axis labels and title
%axes(h_data{1});
box off;
ylabel('');
yticklabels({});
ttl = title('monkey N');
ttl.Position = ttl.Position + [0 .6 0];
xlabel('early evidence $E_t$', 'interpreter','latex');


% add arrows for values of late evidence in conditional psychometric
% curves
for p=1:nPC
    ha = arrow([MaxEvidence+.2 LateEvidenceValue(p)],  [MaxEvidence LateEvidenceValue(p)]);
    arrow(ha,'length',8,'width',2);
    set(ha, 'facecolor',PsychometricCurveColor(p,:),'edgecolor',[1 1 1]);
end

% plot conditional psychometric curve
subplot(4,4,7);
DataPointsCutoff = 20;
plot_conditional_psychometric_curve(IM(s).IntegrationMapData, LateEvidenceValue, DataPointsCutoff, PsychometricColor);

% set axis labels and title
title(''); ylabel(''); legend off;
ylim([0 1]);
yticklabels({});
xlabel('early evidence $E_t$', 'interpreter','latex');


%% add colorbar for integration map
subplot(4,4,4);
cmap = IntegrationMapColor(1,:)+ linspace(0,1,64)'.*diff(IntegrationMapColor,1); % colormap
colormap(gca, cmap);
cb = colorbar;
title(cb, 'p(right)')
cb.Position(1) = .75;
axis off;

%% 6. PANELS D & E: BIASES AND LAPSES OF CONDITIONAL PSYCHOMETRIC CURVES

biases = cell(1,4); % 2 models, then 2 monkeys
biases_se = cell(1,4);
lapses = cell(1,4);
lapses_se = cell(1,4);

% integration model
nLateEvidenceValue = length(LateEvidenceBinCenters);
biases{1} = SensitivityIntegration*LateEvidenceBinCenters; % according to integration model, bias corresponds to late evidence
lapses{1} = [zeros(1,nLateEvidenceValue); .02*ones(1,nLateEvidenceValue)]; % no lapses ( just adding epsilon for visualization)

% snapshot model
biases{2} = zeros(1,nLateEvidenceValue); % according to no-integration model, no bias
lapses{2} = .5*[1-logistic(SensitivityNoIntegration*LateEvidenceBinCenters) ; ...
    logistic(SensitivityNoIntegration*LateEvidenceBinCenters)]; %  lapses are set by late evidence level

% experimental data
for s=1:nAnimal
    PCfit = IM(s).ConditionalPsychometricCurvesData; % parameter fits for psychometric curves conditioned on late evidence
    biases{s+2} = PCfit(1,:,1); % biases fit to monkey 1 data
    biases_se{s+2} = PCfit(1,:,2); % biases fit to monkey 1 data (standard error)
    lapses{s+2} = PCfit(3:4,:,1); % lapses fit to monkey 1 data
    lapses_se{s+2} = PCfit(3:4,:,2); % lapses fit to monkey 1 data
end

for m=1:3 % for both models and experimental data
    % bias panel
    sb(m,4) = subplot(4,4,8+m);
    if m==3 % experimental data
        wu(LateEvidenceBinCenters', biases{m}',biases_se{m}', 'linewidth',2, 'color','k', 'errorstyle','fill');
    else % model
        plot(LateEvidenceBinCenters, biases{m}, 'k', 'linewidth',2);
    end
    xticklabels({}); box off;
    if  m==1
        ylabel('bias \beta');
    else
        yticklabels({});
    end

    % lapse panel
    sb(m,5) = subplot(4,4,12+m);
    if m==3 % experimental data
        wu(LateEvidenceBinCenters, lapses{m}',lapses_se{m}', 'linewidth',2, 'color',IntegrationMapColor, 'errorstyle','fill');
        ylim([-.2 0.5]);
    else % model
        wu(LateEvidenceBinCenters, lapses{m}', [], 'color',IntegrationMapColor, 'linewidth',2);
    end
    xlabel('early evidence $E_t$', 'interpreter','latex');
    if  m==1
        legend({'left','right'});
        ylabel('lapse \pi');
    else
        legend off;
        yticklabels({});
    end
    box off;

end

% use same xlim and ylim for panels across models
sameaxis(sb(:,4));
sameaxis(sb(:,5));

%% 5. ILLUSTRATE FIT OF CONDITIONAL PSYCHOMETRIC CURVE
sb(1,3) = subplot(4,4,12); hold on;

% parameters
LeftLapse = .1;
RightLapse = .2;
Bias = .3;

% define psychometric function with lapses
PsychometricFunction = @(x)  RightLapse+(1-LeftLapse-RightLapse)*logistic(x);

% plot psychometric function
plot(hh,PsychometricFunction(-Bias+hh),'color',mean(PsychometricCurveColor),'linewidth',2);
plot(xlim,.5*[1 1], 'color',.7*[1 1 1]);
plot([0 0],[0 1], 'color',.7*[1 1 1]);

% add dashed line for each parameter
plot(Bias*[1 1], [0 PsychometricFunction(Bias)], 'k--');
plot(hh(1)+[0 1], RightLapse*[1 1], 'k--');
text(hh(1)+1, RightLapse/2, '\pi_R');
plot(hh(end)+[0 -1], (1-LeftLapse)*[1 1], 'k--');
text(hh(end)-1, 1-LeftLapse/2, '\pi_L', 'horizontalalignment','right');
set(gca, 'xtick',Bias, 'xticklabel','\beta');
xlabel('early evidence $E_t$', 'interpreter','latex');
axis tight; ylim([0 1]);

%% 6. ADD CORRELATION OF INTEGRATION MAPS BETWEEN MODEL AND DATA
sb(1,6) = subplot(4,4,16);
model_label = {'integration','snapshot','extrema'};
subject_label = cellfun(@(x) x(end), S.subject_label,'unif',0);

% bar plot
[~,~,hbar] = wu(corr_2dmap, [], 'bar',{subject_label},'color',ModelColor);
ylabel({'model-data','map correlation'});
ylim([.78 1]);
legend off;
subplot_resize([sb(:,4);sb(:,5)],[1.1 1],'right'); % change size of subplot

% add model label
for m=1:length(model_label)
    htext = text(2+hbar.mean(m).XOffset, .8, model_label{m});
    set(htext, 'Rotation',90,'Color',[1 1 1]);
end

%%
width = 20;
height = 20;

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [0 2 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [2 2 width height]);

panelletter(sb(1,:));


%% 7. SUPPLEMENTARY FIGURE: OTHER MONKEY
figure;
s = 2;% second monkey
m = s+2; % corresponding identifier

% plot integration map
h_data(1) = subplot(221);
IM(s).IntegrationMapData = plot_integration_map(IM(s).IntegrationMapData, IntegrationMapColor, isoLevels, EvidenceLabels);

%h_data = copyfigtoaxes( data_analysis_fig,gcf, [3 1], 2, 2, [1 3]);
%axes(h_data{1});
title('');
xlabel('early evidence $E_t$', 'interpreter','latex');
ylabel('early evidence $L_t$', 'interpreter','latex');

% add arrows for conditional psychometric curves
for p=1:nPC
    ha = arrow([MaxEvidence+.2 LateEvidenceValue(p)],  [MaxEvidence LateEvidenceValue(p)]);
    arrow(ha, 'length',8, 'width',2);
    set(ha, 'facecolor',PsychometricCurveColor(p,:),'edgecolor',[1 1 1]);
end

% plot conditional psychometric curves
h_data(2) = subplot(223);
plot_conditional_psychometric_curve(IM(s).IntegrationMapData, LateEvidenceValue, DataPointsCutoff, PsychometricColor);
legend off;

%axes(h_data{2});
title('');
ylim([0 1]);
ylabel('prop. rightwards');


% bias panel
h_data(3) = subplot(222);
wu(LateEvidenceBinCenters, biases{m},biases_se{m}, 'linewidth',2, 'color','k', 'errorstyle','fill');
xlabel('late evidence');
ylabel('bias \beta');

% lapses panel
subplot(224);
wu(LateEvidenceBinCenters, lapses{m}',lapses_se{m}', {{},{'left','right'}},'linewidth',2, 'color',IntegrationMapColor, 'errorstyle','fill');
ylim([-.2 0.5]);
xlabel('late evidence'); ylabel('lapse \pi');

% panel letter
panelletter(h_data,'ABC');