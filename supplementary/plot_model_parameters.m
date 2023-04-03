% SUPPLEMENTARY FIGURE 1: PLOT MODEL PARAMETERS FOR MONKEYS

clear; close all;
animal = 'monkey';
models = { 'integration',... % integration model
    'snapshot1',... % snapshot model
    'extremadetection'}; % extrema detection

nModel = length(models);
AnimalColor = {[.5 1 .5], [.8 .2 .7]}; % color for each animal

% where to load analysis from from
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

% model files
modelfiles = cell(1,nModel);
for m=1:nModel
    modelfiles{m}  = sprintf('%s/%s_%s',models{m}, animal,models{m});
end


%% create figure
figure;
nRow = 1; % one panel row
nCol = 3; % three panel columns
set(gcf, 'name','Supp Fig 1 - Model parameters');


%% Panel A: session-wise parameter in integration model
sp(1) = subplot(nRow,nCol,1); hold on;
title('Integration model');

% load model fit
S = load(fullfile(ModelfitsDir,modelfiles{1}));
subject_label = S.subject_label;
nSubject = length(subject_label);

for s=1:nSubject
    SessionWeight = S.M.regressor(1).Weights(2); % weight structure for session (GUM)
    wu(SessionWeight.PosteriorMean(s,:)', SessionWeight.PosteriorStd(s,:)', 'color',  AnimalColor{s});
end
legend off;
xlabel('session');
ylabel('gain');

%% Panel B: mixture coefficients in snapshot model
sp(2) = subplot(nRow,nCol,2);
title('Snapshot model');

% load model fit
S = load(fullfile(ModelfitsDir,modelfiles{2}));
wu(S.pi_all(1:end-2,:,1),S.pi_all(1:end-2,:,2), 'color',AnimalColor,{{},subject_label});
xlabel('sample');
ylabel({'mixture coefficient'});
ylim([0 .4]);

%% Panel C: extrema detection parameters
sp(3) = subplot(nRow,nCol,3); hold on;
title('Extrema detection');

% load model fit
S = load(fullfile(ModelfitsDir,modelfiles{3}));
wu(S.Tsigma_all(:,:,1)',S.Tsigma_all(:,:,2)','bar', {{'threshold T','std \sigma'},{}},'color',AnimalColor);
legend off;

%% Set paper size
width = 15;
height = 6;

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);

% add panel letter
panelletter(sp);