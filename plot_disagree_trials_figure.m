%% Figure 4:
% plot analysis of disagree trials for monkeys (comparing experimental data
% with integration and extrema detection models)

clear; close all;
animal = 'monkey';


models = {'integration', ... % integration model
    'extremadetection',...  % extrema detection with fluctuating threshold (fixed lapses)
    'extremadetection_lastsample'}; % extrema detection with last sample (fixed lapses)

model_label = {'integration','extrema detection', 'ED last sample'};

ModelColor = {[1 0 0];[0 0 .8]; [0 .5 0]}; % red for integration, blue for snapshot, green for extrema detection

ModelfitsDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');


%% Load analysis files
disagree_file = sprintf('disagree/%s_disagree.mat',animal);
mlabel = [{'animal'}, model_label(1:2)];

%load analysis of experimental responses in disagree trials
S = load(fullfile(ModelfitsDir,disagree_file));

nSubject = S.nSubject;
subject_label = S.subject_label;

pTotalEvidenceChoiceOnDisagree = zeros(nSubject,3);
pTotalEvidenceChoiceOnDisagreeCI = nan(nSubject,3,2);

for m=1:3
    if m==1
        for s=1:S.nSubject
            fprintf('Number of disagree trials for %s: %d\n',S.subject_label{s}, S.nDisagreeTrials(s));
        end
    else
        %load analysis of simulated responses in disagree trials
        S = load(fullfile(ModelfitsDir,models{m-1}, disagree_file));

        % 95% confidence intervals from boostraps
        pTotalEvidenceChoiceOnDisagreeCI(:,m,:) = permute(S.pTotalEvidenceChoiceOnDisagreeCI,[2 3 1]);
        
    end
    pTotalEvidenceChoiceOnDisagree(:,m) = S.pTotalEvidenceChoiceOnDisagree; % animal x side x weighting x model
end

%% PLOT FIGURE FOR DISAGREE TRIALS (FIGURE 3)
figure;
nRow = 1;
nCol = 2;

for s=1:nSubject
    sb(s) = subplot(nRow,nCol,s); hold on;

    % bar plot
    LowerBar = pTotalEvidenceChoiceOnDisagree(s,:) - pTotalEvidenceChoiceOnDisagreeCI(s,:,1);
     UpperBar = pTotalEvidenceChoiceOnDisagreeCI(s,:,2) - pTotalEvidenceChoiceOnDisagree(s,:);
    wu([],pTotalEvidenceChoiceOnDisagree(s,:),LowerBar,UpperBar,...
        'bar','color',{.4*[1 1 1],ModelColor{1} ModelColor{3}},'basevalue',0.5, {{},mlabel});

    title(subject_label{s});
    if s==1
        legend off;
    end
    xticks([]);

    if s==1 %ceil(s/nColSubplot)==nRowSubplot
        ylabel('p. majority-aligned choices');
    else
        set(gca, 'yticklabel',{});
    end
    plot(xlim, .5*[1 1],'color',.7*[1 1 1]);
end
sameaxis(sb);
sb = sb(1);

%%
width = 10;
height = 4;

% Position plot on the screen for drawing
set(gcf, 'Units', 'centimeters', 'Position', [2 4 width height]);

% Position plot on the paper for printing
set(gcf, 'PaperUnits', 'centimeters', 'PaperPositionMode', 'manual',...
    'PaperSize', [width height], 'PaperPosition', [0 0 width height]);