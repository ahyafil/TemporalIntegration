function S = IntegrationMap(X, Y, bnd, dx, sigma)
% S = IntegrationMap(X, Y, [,bnd] [,dx][,sigma]) computes an integration
% map
% X: ntrial x 2 matrix (evidence from source 1 in column 1, evidence from source 2 in column 2)
% Y: corresponding binary vector of choices (can also be probability for
% choice 1 as defined from a model)
% dx :binsize
% bnd: bounds of analysis (i.e. [Xmin Xmax])
% sigma: width of gaussian kernel for smoothing

if nargin<3
    bnd = [min(X(:)) max(X(:))];
end
if nargin<4
    dx = diff(bnd)/10; % by default use 10 bins
end
if nargin<5
    sigma = 5*dx; % width of kernel (in bins)
end

nConv = sigma*4; % length of convolution kernel
sigma = sigma/dx;% convert to bins
nConvBins = sigma*4; % length of convolution kernel (in bins, long enough to go back to 0)

if nargin<7
    EvidenceLabels =  {'early evidence', 'late evidence'};
end

% Define bin boundaries
binBoundaries = -(bnd+nConv):dx:(bnd+nConv); % we extend the boundaries with nConv bins prior to applying the gaussian kernel, we will trim later on
binBoundariesInf = [-Inf binBoundaries Inf];

% center of bins
binCenters = [binBoundaries(1)-dx/2 binBoundaries+dx/2];

%% compute number of trials and sum responses in each bin
nBin = length(binCenters);
nDatapoints = zeros(nBin);
M = zeros(nBin);

% compute corresponding bin for source 2
bin2 = zeros(size(X,1),1);
for j=1:nBin
    mask = X(:,2)>= binBoundariesInf(j) & X(:,2)< binBoundariesInf(j+1);
    bin2(mask) = j;
end

for i=1:nBin
    mask1 = X(:,1)>= binBoundariesInf(i) & X(:,1)< binBoundariesInf(i+1); % data points in this bin for evidence 1
    for j=1:nBin
        mask = mask1 & bin2==j; % data points in this 2D bin

        % number of data points in this bin
        nDatapoints(i,j) = sum(mask);

        % sum of responses for data points in this bin
        M(i,j) = sum(Y(mask));
    end
end

%% smooth by convolution with gaussian kernel

% define gaussian kernel
filter_1d = exp(-(-nConvBins:nConvBins).^2/2/sigma^2);

% 2D filter
filter_2d = filter_1d'*filter_1d;

% apply 2D convolution
M = conv2(M,filter_2d,'same');
nDatapoints = conv2(nDatapoints,filter_2d,'same');

%% trim data outside of bounds
inBoundaries = (binCenters>-bnd) & (binCenters<bnd);
M = M(inBoundaries,inBoundaries);
nDatapoints = nDatapoints(inBoundaries,inBoundaries);
binCenters = binCenters(inBoundaries);

%% compute mean value (as sum divided by number of data points)
M = M./nDatapoints;

%% add everythin in structure
S = struct;
S.IntegrationMap = M;
S.nDatapoints = nDatapoints;
S.dx = dx;
S.bnd = bnd;
S.binCenters = binCenters;
S.nConv = nConvBins;
S.sigma = sigma;
S.EvidenceLabels = EvidenceLabels;

end