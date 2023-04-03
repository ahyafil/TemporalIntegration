function plot_conditional_psychometric_curve(S, yValue, DataPointsCutoff, Colour)
% plot_conditional_psychometric_curve(S, yValue)
% plot_conditional_psychometric_curve(S, yValue, DataPointsCutoff, Colour)
%
% plots conditional psychometric curves (integration map analysis)

if nargin<3
    DataPointsCutoff = 20;
end

if nargin<4
    Colour = [0 0 0; 1 0 0]; % color of psychometric curves (graded from black to red)
end

if ~isfield( S, 'EvidenceLabels')
    S.EvidenceLabels =  {'early evidence', 'late evidence'};
end

% do not plot when less than a certain number of data points
EnoughDataPoints = S.nDatapoints>= DataPointsCutoff;
EnoughDataPoints = double(EnoughDataPoints);
EnoughDataPoints(EnoughDataPoints==0) = nan;

% vector of values on x-axis
xBinEdges = -S.bnd:S.dx:S.bnd;
xBinCenter = xBinEdges(1:end-1) + S.dx/2;

% number of conditional psychometric curves
nPC = length(yValue);

% find indices corresponding to each value of yValue
yIndex = zeros(1,nPC);
for i=1:nPC
    yIndex(i) = find(abs(S.binCenters - yValue(i))<1e-3);
end

% graded colours for psychometric curves
ColormapPC = Colour(1,:) + linspace(0,1,nPC)'*diff(Colour);

Evidence2Label = string(S.EvidenceLabels{2})+yValue;

% all conditional psychometric curves (removing points with non enough data points)
PC =  EnoughDataPoints.* S.IntegrationMap;

% plot
hold on;
for i=1:nPC
    plot(xBinCenter, PC(:,yIndex(i)), 'color',ColormapPC(i,:),'LineWidth',2);
end
legend(Evidence2Label);
xlabel(S.EvidenceLabels{1});
ylabel('p(choice)');

end