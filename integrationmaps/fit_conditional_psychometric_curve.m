function params = fit_conditional_psychometric_curve(X,Y, x2Value)
% params = fit_conditional_psychometric_curve(X,Y, x2Value)
% fits conditional psychometric curves. X is a n-by-2 matrix of stimulus
% evidence, Y is a n-by-1 vector of binary responses. x2Value is a vector
% of V values of X2. fit_conditional_psychometric_curve fits psychometric
% curves with lapses of Y vs X(:,1), for different values of X(:,2).
%
% params is of size 4 by nValue by 3 array
% dimension 1: bias / sensitivity / lapse 0 /lapse 1
% dimension 2: value of X2 we condition on
% dimension 3: ML param / s.e.m. / p-value
%
% see also plot_conditional_psychometric_curve

PLoptions = struct('ninit',1,'alpha',[1 1 3],'maxiter',5000);

nPC = length(x2Value)-1;
params = zeros(4,3,nPC);

for i=1:nPC
    this_trial = X(:,2)>=x2Value(i) & X(:,2)<x2Value(i+1);

    % fit simple probit regression with lapses
    S = probitlapse(X(this_trial,1),Y(this_trial),0,2,[],PLoptions);
    params(:,:,i) = S.joint;
end

% (bias / lapses ) x late evidence x (ML param / s.e.m.)
params = permute(params, [1 3 2]);
end