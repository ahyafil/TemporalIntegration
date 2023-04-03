function S = plot_integration_map(S, Colours, isoLevels, EvidenceLabels)
% plot_integration_map(S)
% plots an integration map, i.e. a 2D map for the 
%
% plot_integration_map(S, Colours) to specify the 2 extreme colours to use
% to create the colormap. Colous must be a 2-by-3 matrix of RGB values.
%
% plot_integration_map(S, Colours, isoLevels) to specify the levels at
% which to draw contour lines. Default: [.15 .3 .5 .7 .85]
%
% plot_integration_map(S, Colours, isoLevels, EvidenceLabels) to provide
% labels for 1st and 2nd sources of information. Default: {'early evidence', 'late evidence'}
%
if nargin<2 || isempty(Colours)
    Colours = [103 169 221; 241 135 34]/255; % gradient between red and blue
else
    assert(ismatrix(Colours) && size(Colours,1)==2 && size(Colours,2)==3);
end
if nargin<4 
    EvidenceLabels =  {'early evidence', 'late evidence'};
end
if nargin<3 || isempty(isoLevels)
    isoLevels = [.15 .3 .5 .7 .85]; % level for isolines
end

M = S.IntegrationMap;
M(isnan(M)) = 0;

%convert mean response and number of trials to RGB
C = zeros(size(M,1),size(M,2),3);
for c=1:3 % for R,G,B
    C(:,:,c) = 1- (Colours(2,c)-Colours(1,c)) * M - Colours(1,c);
end

% hue intensity (clearer if less datapoints)
int = min(S.nDatapoints/20,1);
C = 1- C.* int;

%plot image
image(S.binCenters, S.binCenters, permute(C,[2 1 3]));

axis xy; hold on;

%draw contour lines
midValue = round( (length(isoLevels)+1) /2); % mid-value (usually .5) with thicker value
ThickIsoLine = isoLevels(midValue);
isoLevels(midValue) = [];
contour(S.binCenters, S.binCenters, M', isoLevels,'k', 'linewidth',.5);
contour(S.binCenters, S.binCenters, M',ThickIsoLine*[1 1],'k', 'LineWidth',1); % thicker line

xlabel(EvidenceLabels{1}); ylabel(EvidenceLabels{2});

S.EvidenceLabels = EvidenceLabels;

end