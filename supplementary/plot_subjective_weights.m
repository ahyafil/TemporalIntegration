%% PLOT SUPP FIGURE X
% plot mapping from stimulus dimension to perceptual weight for monkeys, integration model and extrema detection

clear; close all;

nModel = 2; % integration model and extrema detection
nCol = nModel+1; % first two cols for models, then experimental data

mapping_model = 'subjectiveweights';
animal = 'monkey';
nSubject = 2;

models = {'integration', ... % integration model
    'extremadetection' };   % extrema detection with fluctuating threshold
model_label = {'integration','extrema detection'};

ModelColor = {[1 0 0];[0 .5 0]; 'k'}; % red for integration, green for extrema detection, black for experimental data

% directory where model fits are stored
ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

figure;


%% 1. PLOT SUBJECTIVE WEIGHTS FOR SIMULATED DATA
for m=1:nModel % for all models
    % load file with fitted model with subjective weights
    modelfile_subjective  = sprintf('%s/%s/%s_%s',models{m}, mapping_model, animal,mapping_model);
    Smap = load(fullfile(ModelfitsDir,modelfile_subjective));

    % fitted regressor
    Rmap = Smap.M.regressor(1);

    if m==2 % for extrema detection model, retrieve threshold parameter
        % load model file
        modelfile  = sprintf('%s/%s_%s',models{m}, animal, models{m});
        S = load(fullfile(ModelfitsDir,modelfile));

        % fitted value of detection threshold
        T = S.Tsigma_all(:,1,1);
    end

    for s=1:nSubject % loop through subjects (2)
        % create subplot
        sb2(s,m) = subplot(nSubject, nCol,m+(s-1)*nCol);
        hold on;

        % plot weights for each value of pulse
        W = Rmap.Weights(2); % weight structure (GUM)
        wu(W.PosteriorMean(s,:)',W.PosteriorStd(s,:)', 'color',ModelColor{m});
        if m==2 % extrema detection
            % plot vertical line at threshold parameter
            plot(T(s)*[1 1], ylim, 'k--');
        end
        if s==1
            title(model_label{m});
        else
            xlabel('sample evidence');
        end

        ylim([0 2]);
        if m==1
            ylabel({'impact on','decision'});
        end
        xlim([0 13]);
    end

end

%% 2. PLOT SUBJECTIVE WEIGHTS FOR EXPERIMENTAL DATA

% load analysis of subjective weights on experimental data
modelfile_subjective  = sprintf('%s/%s_%s',mapping_model, animal,mapping_model);
Smap = load(fullfile(ModelfitsDir,modelfile_subjective));
Rmap = Smap.M.regressor(1); % fitted regressor
for s=1:nSubject

    sb2(s,m+1) = subplot(nSubject, nCol,nModel+1+(s-1)*nCol); hold on;
    if s==1
        title('Experimental data');
    else
        xlabel('sample evidence');
    end

    % plot weights for each value of pulse
    W = Rmap.Weights(2); % weight structure (GUM)
    wu(W.PosteriorMean(s,:)',W.PosteriorStd(s,:)', 'color',ModelColor{3});

    % axis limits
    xlim([0 13]);
    sameaxis(sb2(s,:));
    if s==1
        ylim(sb2(s,:),[0 2]);
    else
        ylim(sb2(s,:),[0 1.5]);
    end
end

% add panel letter
panelletter(sb2(1,:),'ABC');
