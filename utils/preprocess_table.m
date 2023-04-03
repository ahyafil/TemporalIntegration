function [T, nSample, subject_id, nSubject, nSamples_max, nSamples_mean] = preprocess_table(T)
% [T, nSample, subject_id, nSubject] = preprocess_table(T)
% preprocesses table with monkey behavioral data for temporal integration project

% use string for subject ID (for monkeys)
if iscell(T.subject)
    T.subject  = string(T.subject);
end

% list of subjects
subject_id = unique(T.subject);

% number of subjects
nSubject = length(subject_id);

if nSubject==2 % invert order for monkeys
subject_id = flipud(subject_id);
end

% identify how many samples in stimuli (i.e. how many columns named
% 'stimulus_i')
nSample = 1;
while ismember("stimulus_"+nSample,T.Properties.VariableNames)
    nSample = nSample+1;
end
nSample = nSample-1;


% merge columns for stimulus into single variable (one row vector for each
% trial)
T = mergevars(T,"stimulus_"+(1:nSample),'NewVariableName',"stimulus");

%% detect maximum and mean number of subjects for each subject
nSamples_max = zeros(1,nSubject); % maximum number of samples for each subject
nSamples_mean = zeros(1,nSubject); % mean number of samples for each subject

for s=1:nSubject
    this_stimulus = T.stimulus(T.subject==subject_id(s),:); % stimulus data for this subject
    nSamples_trial = sum(~isnan(this_stimulus),2)'; % number of frames for all trials
    nSamples_max(s) = max(nSamples_trial);
    nSamples_mean(s) = mean(nSamples_trial);
end

end