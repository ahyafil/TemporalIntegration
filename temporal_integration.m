%% fit different for temporal integration vs non-temporal integration

clear; close all;

% select animal: 'rat','monkey', 'human'
animal = 'monkey';

UseSimulations = false; % true to analyze simulated data, false to analyze experimental data
SimulatedModel = 'extremadetection';  %  'extremadetection' or 'integration': name of model to simulate and analyze (for disagree trials and subjective weights)

% define variant of extrema detection and snapshot model
LapseParameters = [.01 .01]; % lapse parameters for extrema detection and snapshot model (if using fixed lapses)
UseFixedLapseSnapshot = true; %  (default: true) whether lapses are fixed parameters for snapshot model
snapshot_model_mode = 'smooth'; % 'deterministic' (default) or 'smooth': whether in snapshot model always picks the side of the attended sample
singleSnapshot = false; % true to fit only snapshot model with span of 1, false to fit models with all possible spans

UseFixedLapseExtremaDetection = true; % (default: true) whether lapses are fixed parameters for extrema detection model 
LastSample = true; % (default: false) in extrema detection, whether to use last sample information if threshold has not been reached during sequence
UseFixedThresholdExtremaDetection = false; % (default: true) whether the threshold is fixed throughout the sample sequence or is allowed to vary sample to sample


% Define which analysis to perform
anals.integration = 0; ~UseSimulations; % integration model
anals.snapshot = 0;~UseSimulations; % snapshot model
anals.extremadetection = ~UseSimulations; % extrema detection
anals.subjectiveweights = 0; strcmp(animal, 'monkey'); % capturing subjective weights, mapping from stimulus space to evidence space
anals.disagree = 0; % analysis of disagree trials

Nbtstrp = 100; % number of bootstraps for integration maps
nRepPosterior = 2000; % number of simulations of integration/extrema detection model from weights sampled for posterior (for error bars of disagree trials analysis)

dxIntegrationMap = .2; % bin size for integration maps
sigmaIntegrationMap = .5; % smoothing parameter for integration maps


%numWorkers = 10; % how many workers to use (parallel processing)
%mono_tool('parcluster','',numWorkers); % define parfor clusters

% a regressor modulating stimulus sensitivity is included in monkey and rat
% data
SessionModulation = ~strcmp(animal, 'human');

% defines the directory where results of analysis should be saved
ModelSaveDir = fullfile(TemporalIntegrationRootDirectory,'modelfits');

%% LOAD DATA
csvfile = fullfile(TemporalIntegrationRootDirectory, 'data', animal);
T = readtable([csvfile '.csv']);
fprintf('loading behavioral data from %s.csv\n', csvfile);

%% PREPROCESS DATA
[T, nSample, subject_id, nSubject, nSamples_max, nSamples_mean] = preprocess_table(T);

% define subject label
switch animal
    case 'monkey'
        subject_label = "Monkey "+subject_id;
    case 'rat'
        subject_label = subject_id;
    case 'human'
        % select only valid trials
        T = T(T.valid==1,:);
        subject_label = "S" + subject_id;
    otherwise
        error('don''t have data to analyze for this animal');
end


%% USE SIMULATED DATA IF REQUIRED
resp_posterior= cell(1,nSubject);
if UseSimulations
    fprintf('using simulated data from %s model\n',SimulatedModel);

    SimulationFile =[SimulatedModel '/%s_' SimulatedModel]; % mat file with fitted model to simulate

    SimulationFile = fullfile(ModelSaveDir,SimulationFile);

    mfile = load(strrep(SimulationFile,'%s',animal)); % load model


    for s=1:nSubject % loop through subjects

        if any(strcmp(SimulatedModel,{'integration'})) % GUM (integration model)

            pModel = mfile.M.Predictions(s).Expected; % probability of response according to model

            % sampled responses from posterior over weights
            if strcmp(animal,'monkey')
            resp_posterior{s} = mfile.M.Predictions(s).sample_from_posterior;
            end

        elseif any(strcmp(SimulatedModel,{'extremadetection'})) % extrema detection
            pModel = mfile.S_all(s).pModel; % probability of response according to model

            if strcmp(animal,'monkey')
            resp_posterior{s} = mfile.S_all(s).resp_posterior;
            end

        else % snapshot
            pModel = mfile.S_all(s).Y; % probability of response according to model

        end

        nTrial = length(pModel); % number of trials for this subject

        % generate random response from model (using prob of response from each trial)
        simulated_resp = pModel > rand(nTrial,1);

        % trial corresponding to this subject
        mask = T.subject== subject_id(s);

        % replace in table
        T.resp(mask) = simulated_resp;
    end

    % move to appropriate directory
    if ~isdir([SimulatedModel '_simul'])
        mkdir([SimulatedModel '_simul']);
    end
    cd([SimulatedModel '_simul']);

    ModelSaveDir = [ModelSaveDir '/' SimulatedModel];
    if ~isdir(ModelSaveDir)
        mkdir(ModelSaveDir);
    end
end


%% file for integration model, serving as basis for early/late integration analysis
IntegrationModelFile =  fullfile(ModelSaveDir,'integration',[animal '_integration.mat']);

%% define cut-off between early and late part of stimulus
cutoff_earlylate = floor(nSamples_mean/2);

switch animal % boundaries of integration map (depends on distribution of stimuli in each paradigm)
    case 'monkey'
        boundaries_IntegrationMap = 2.5;
    case 'human'
        boundaries_IntegrationMap = 1.5;
    case 'rat'
        boundaries_IntegrationMap = 2;
end

% parameters for Generalized Unrestriced Models (GUMs)
gum_options = struct;
gum_options.observations = 'binomial';

gum_fit_options = struct;
gum_fit_options.initialpoints = 3; % number of initial points (just to speed up, they all converge to the same points)
gum_fit_options.verbose = 'on'; 'full';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% START FITTING MODELS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5



%% 1. FIT TEMPORAL INTEGRATION MODEL
if anals.integration
    an = struct('name','integration','short','integration');

    % create directory if needed to store results
    if ~isdir(fullfile(ModelSaveDir,an.short)), mkdir(fullfile(ModelSaveDir,an.short)); end

    subdir = [an.short '_' animal '/'];  % name of subdirectory for figures
    fprintf('%s:\n', upper(an.name));

    % loop through all subjects
    for s=1:nSubject
        fprintf('%s, %s %d/%d:\n',an.name, animal, s,nSubject);

        Tsub = T(T.subject == subject_id(s),:); % table for this subject

        % simple stimulus regressor
        R = regressor(Tsub.stimulus(:,1:nSamples_max(s)), 'linear','variance',1, 'label','stimulus');

        % add session-dependent regressor
        if SessionModulation

            Sess = regressor(Tsub.session, 'categorical','label','block','variance',1);

            Sess.HP(1).LB = .2; % lower bound on hyperparameter for block

            % build model as multiplication of sample evidence and session regressors
            R = R * Sess;
        end

        % create GUM model
        M(s) = gum(R, Tsub.resp, gum_options); %,'label','linear regression');

        options = gum_fit_options;

        % infer model (compute posterior over weights)
        M(s) = M(s).infer(options);

        % compute probability of rightward choices for each trial according
        % to estimated weights
        M(s) = M(s).ExpectedValue;

        if strcmp(animal,'monkey')
        M(s) = M(s).Sample_Observations_From_Posterior(nRepPosterior);
        end
        M(s) = M(s).clear_data;
    end

    % concatenate models over subjects
    M = M.concatenate_over_models(true);

    %% FOR INTEGRATION MAP: COMPUTE EARLY AND LATE PART OF EVIDENCE

    for s=1:nSubject
        Tsub = T(T.subject == subject_id(s),:); % table for this subject

        id_early = 1:cutoff_earlylate(s); % all samples considered early
        id_late = cutoff_earlylate(s)+1:nSamples_max(s); % all samples considered late

        sensory_weights = M.regressor(1).Weights(1).PosteriorMean; % MAP weights for sensory regressor
        early_weights = sensory_weights(s,id_early); % weights for early samples
        late_weights = sensory_weights(s,id_late); % weights for late samples
        if SessionModulation % if we use modulation of weights by session, apply that modulatino
            session_weights = M.regressor(1).Weights(2).PosteriorMean(s,:); %session-wise MAP weights
            session_weights = session_weights( Tsub.session)'; % corresponding modulation for each trial
        else
            session_weights = 1;
        end

        % compute early and evidence as sum of evidence for corresponding
        % samples weighted by corresponding weights
        M.Predictions(s).early =  nansum(early_weights.*Tsub.stimulus(:,id_early),2) .* session_weights; % weighted evidence frmo first 3 samples
        M.Predictions(s).late =  nansum(late_weights.*Tsub.stimulus(:,id_late),2) .* session_weights; % weighted evidence frmo last 4 samples
    end

    %% SAVE MODEL
    modelfile = fullfile(ModelSaveDir, an.short,[animal '_' an.short]);
    save(modelfile, 'M','animal','nSubject','subject_label');
    fprintf('saved results for %s regression for all %ss in %s\n',an.name, animal, modelfile);

    %% APPLY INTEGRATION MAP ANALYSIS
    if ~UseSimulations
        M = run_integration_map(IntegrationModelFile,T,M, subject_id, boundaries_IntegrationMap, Nbtstrp, dxIntegrationMap,sigmaIntegrationMap);

        %update file
        save(modelfile, 'M','animal','nSubject','subject_label','subdir');
        fprintf('updated result file\n');
    end
end


%% 2. FIT SNAPSHOT MODEL
if anals.snapshot
    % whether uses deterministic variant (i.e. no noise in threshold)
    is_deterministic = strcmp(snapshot_model_mode,'deterministic');

    if singleSnapshot
        MaxSpan = 1; % fit only with span=1
    else
        MaxSpan = min(nSample-1,6); % fit with all possible: 1, 2, ... (max 6)
    end

    for span = 1:MaxSpan %  loop through all values of span (i.e.number of samples in each snapshot)
        fprintf('snapshot of %d samples:\n',span);

        %% define name of analysis
        an.name = sprintf('snapshot%d', span);
        an.short = an.name;


        if ~is_deterministic
            an.short = [an.short '_stochastic'];
        end
        if ~UseFixedLapseSnapshot
            an.name = [an.name ' with free lapses'];
            an.short = [an.short '_freelapse'];
        end

        if ~isdir(fullfile(ModelSaveDir,an.short)), mkdir(fullfile(ModelSaveDir,an.short)); end
        fprintf('%s:\n', upper(an.name));
        subdir = [an.short '_' animal '/'];  % name of subdirectory for figures


        % if modulation of sensory information, recover session-wise gain
        % from integration model fit
        if SessionModulation
            F = load(IntegrationModelFile);
        end

        % parameter for model fitting
        param = struct;
        if UseFixedLapseSnapshot
            param.FixedLapse = LapseParameters;
        end
        U_all = cell(1+~is_deterministic,nSubject,2);
        clear S_all

        for s=1:nSubject % loop through subjects
            fprintf('%s, %s %d/%d:\n',an.name,animal,s,nSubject);
            mask = T.subject==subject_id(s); % trials for this subject

            if strcmp(animal,'rat') && ~is_deterministic
                param.maxiter = 1000;
            end

            stimulus = T.stimulus(mask,1:nSamples_max(s));
            stimulus(isnan(stimulus)) = 0; % if no sample because sequence is shorter, then no info

            if SessionModulation
                % session gain for each trial
                SessionGain = F.M.regressor(1).Weights(2).PosteriorMean(s,T.session(mask))'; % session gain for each trial

                % multiply evidence for each sample by session gain
                stimulus = stimulus .*  SessionGain;
            end

            % fit snapshot model
            [U, pi_snap, S]= snapshot_fit( stimulus, T.resp(mask),  snapshot_model_mode, span, param);

            % store fitted parameters in cell array
            U_all{1,s,1} = pi_snap'; % mixture components
            U_all{1,s,2} = S.mixt_se; % standard error for mixture components
            if ~is_deterministic
                U_all{2,s,1} = U; % sensitivity parameter
                U_all{2,s,2} = S.w_se; % standard error for sensitivity parameter
            end
            S_all(s) = S;
        end


        U_all = cellfun(@(x) [x nan(1,nSample+2-length(x))], U_all, 'unif',0);
        if ~is_deterministic
            k_temp_all = catcell(1,permute(U_all(2,:,:),[4 1 2 3]),2,true); % temporal kernels (sample x mixture x animal x mean/ste)
        else
            k_temp_all = [];
        end
        pi_all =  permute(catcell(1,U_all(1,:,:),2,true),[2 1 3]); % mixture coefficients

        % save data across animals
        modelfile = fullfile(ModelSaveDir, an.short,[animal '_' an.short]);
        save(modelfile, '*_all','subdir', 'animal','subject_label');
        fprintf('saved results for %s regression for all %ss in %s\n',an.name, animal, modelfile);

        %% compute integration map according to model
        if ~UseSimulations
            FitConditionalPC = (span==1) && is_deterministic && UseFixedLapseSnapshot; %whether we fit conditional PCs
            S_all = run_integration_map(IntegrationModelFile,T,S_all, subject_id, boundaries_IntegrationMap, Nbtstrp, dxIntegrationMap, sigmaIntegrationMap,FitConditionalPC);
            save(modelfile, '*_all','subdir', 'animal','subject_label');
        end
    end
end


%% 3. FIT EXTREMA DETECTION MODEL
if anals.extremadetection
    fprintf('extrema detection model:\n');

    % define model name
    an = struct('name', 'extrema detection','short', 'extremadetection');
    if ~UseFixedLapseExtremaDetection
        an.name = [an.name ' with free lapses'];
        an.short = [an.short '_freelapse'];
    end
    if LastSample
        an.name = [an.name ', last sample rule'];
        an.short = [an.short '_lastsample'];
    end
    if ~UseFixedThresholdExtremaDetection
 an.name = [an.name ', varying threshold'];
        an.short = [an.short '_varyingthreshold'];

    end
    fprintf('%s:\n', upper(an.name));
    if ~isdir(fullfile(ModelSaveDir,an.short)), mkdir(fullfile(ModelSaveDir,an.short)); end
    subdir = [an.short '_' animal '/'];  % name of subdirectory for figures

    % parameters for model fitting
    param = struct;
    if UseFixedLapseExtremaDetection
        param.FixedLapse = LapseParameters;
    end
    param.LastSample = LastSample;
    param.FixedThreshold = UseFixedThresholdExtremaDetection;

    U_all = cell(2,nSubject,2);
    clear S_all

    % if modulation of sensory information, recover session-wise gain
    % from integration model fit
    if SessionModulation
        F = load(IntegrationModelFile);
    end

    %% loop through subjects
    for s=1:nSubject
        fprintf('%s, %s %d/%d:\n',an.name,animal, s,nSubject);

        % mask for trials for this subject
        mask = T.subject == subject_id(s);

        stimulus = T.stimulus(mask,1:nSamples_max(s)); % sample evidence for this subject (trial x sample)
        if SessionModulation
            SessionGain = F.M.regressor(1).Weights(2).PosteriorMean(s,T.session(mask))'; % session gain for each trial
            stimulus = stimulus .*  SessionGain;
        end

        % FIT EXTREMA DETECTION MODEL
        [Threshold, sigma, lapse, S] = extremadetection_fit(stimulus, T.resp(mask), param);

        U_all{1,s,1} = [Threshold sigma];
        U_all{1,s,2} = [S.T_se' S.sigma_se];
        U_all{2,s,1} = lapse;
        U_all{2,s,2} = S.lapse_se';

        % generate simulations with parameters from the approximated
        % posterior
        if strcmp(animal,'monkey') && UseFixedLapseExtremaDetection && ~LastSample && UseFixedThresholdExtremaDetection
        this_resp_posterior = zeros(nRepPosterior,sum(mask));
        for r=1:nRepPosterior
            % sample set of parameters (mu, sigma, lapse) from Laplace
            % approximate posterior
            if ~UseFixedLapseExtremaDetection % free lapse parameters
                MAP_w = [Threshold,sigma,lapse]; % map parameters
                ParametersSampled = mvnrnd(MAP_w, S.PosteriorCovariance); % sample from multivariate gaussian
                LapseSampled = min(max(ParametersSampled(3:4),0),1);
            else % fixed lapse parameters
                MAP_w = [Threshold,sigma]; % map parameters
                ParametersSampled = mvnrnd(MAP_w, S.PosteriorCovariance);
                LapseSampled = LapseParameters;
            end

            % simulate extrema detection model with these parameters
            this_resp_posterior(r,:) = extremadetection_rnd(stimulus,ParametersSampled(1), ParametersSampled(2), LastSample, LapseSampled);
        end
        S.resp_posterior = this_resp_posterior;
        end

        S_all(s) = S;
    end

    Tsigma_all =  catcell(1,permute(U_all(1,:,:),[2 1 3]),1,true); % animal x treshold/noise x mean/se
    lapse_all = catcell(1,permute(U_all(2,:,:),[2 1 3]),1,true); % animal x  resp x mean/se

    % save data across sessions
    modelfile = fullfile(ModelSaveDir, an.short,[animal '_' an.short]);
    save(modelfile, '*_all','subdir');
    fprintf('saved results for %s regression for all %ss in %s\n',an.name, animal, modelfile);

    %% APPLY INTEGRATION MAPS
    if ~UseSimulations
        FitConditionalPC = UseFixedLapseExtremaDetection && ~LastSample && UseFixedThresholdExtremaDetection; % whether we fit conditional psychometric curves
        S_all = run_integration_map(IntegrationModelFile,T,S_all, subject_id, boundaries_IntegrationMap, Nbtstrp, dxIntegrationMap,sigmaIntegrationMap, FitConditionalPC);
        save(modelfile, '*_all','subdir');
    end
end

%% 4. SUBJECTIVE WEIGHTS ANALYSIS (testing linearity of pulses)
if anals.subjectiveweights
    an = struct('name','subjective weights','short','subjectiveweights');
    fprintf('%s:\n', upper(an.name));

    if ~isdir(fullfile(ModelSaveDir,an.short)), mkdir(fullfile(ModelSaveDir,an.short)); end
    subdir = [an.short '_' animal '/'];  % name of subdirectory for figures

    %% loop through subjects
    for s=1:nSubject
        fprintf('%s, %s %d/%d:\n',an.name,animal,s,nSubject);

        Tsub = T(T.subject == subject_id(s),:); % table for this subject


        % build regressor for subjective weights
        R = mapping_regressor(Tsub.stimulus, animal);

        if  SessionModulation
            % build session categorical regressor
            Sess = regressor(Tsub.session, 'categorical','label','session');

            % this regressor multiplies the stimulus regressor
            R = R * Sess;
        end

        % build GUM
        M(s) = gum(R, Tsub.resp, gum_options);

        % infer GUM parameters
        options = gum_fit_options;
        M(s) = M(s).infer(options);

        M(s) = M(s).ExpectedValue;
        M(s) = M(s).clear_data;
    end

    M = M.concatenate_over_models(true);

    %% save data
    modelfile = fullfile(ModelSaveDir, an.short,[animal '_' an.short]);
    save(modelfile, 'M','nSubject','subject_label','subdir');
    fprintf('saved results for %s regression for all %ss in %s\n',an.name, animal, modelfile);
end


%% 5. ANALYSIS OF DISAGREE TRIALS (WHETHER CHOICES DEPEND MORE ON STRONG SENSORY SAMPLES OR MANY WEAK)
if anals.disagree
    an = struct('name','disagree trials','short','disagree');
    fprintf('%s:\n', upper(an.name));

    if ~isdir(fullfile(ModelSaveDir,an.short)), mkdir(fullfile(ModelSaveDir,an.short)); end

    pTotalEvidenceChoiceOnDisagree = zeros(1,nSubject);
    pTotalEvidenceChoiceOnDisagreeSem = zeros(1,nSubject);

    nDisagreeTrials = zeros(1,nSubject);

    propchoice_samp = [];

    %% loop through subjects
    for s=1:nSubject

        % mask for trials for this subject
        mask = T.subject == subject_id(s);

        % stimulus matrix for this subject
        stimulus = T.stimulus(mask,:);

        % total evidence over all trial samples
        TotalEvidence = nansum(stimulus,2);

        % choice according to total Evidence
        ChoiceTotalEvidence = sign(TotalEvidence);

        % largest evidence reached on each trial
        LargestEvidence = max(abs(stimulus),[],2);

        % choice according to largest evidence
        LargestEvidence0 = any(stimulus==-LargestEvidence,2); % if max strength reached by neg stim
        LargestEvidence1 = any(stimulus==LargestEvidence,2); % if max strength reached by pos stim
        ChoiceLargestEvidence = LargestEvidence1 - LargestEvidence0; % choice of stronger: -1 if neg-stim, +1 if pos-stim, 0 if both

        % disagree trials
        disagree_trials = (ChoiceTotalEvidence==-ChoiceLargestEvidence) & (ChoiceLargestEvidence>0); %  (first 1/-1, then -1/1)

        % number of disagree trials
        nDisagreeTrials(s) = sum(disagree_trials);

        % proportion of choice agligned with total evidence on disagree trials
        resp = 2*T.resp(mask)-1; % resp as -1/1
        ActualChoiceTotalEvidence = (resp == ChoiceTotalEvidence); % whether subject response is aligned with total evidence
        pTotalEvidenceChoiceOnDisagree(s) = mean(ActualChoiceTotalEvidence(disagree_trials));

        % when running analysis on simulated data, repeat analysis
        % for other simulations from posterior
        if ~isempty(resp_posterior{s})
            % compute p(choice) on disagree trials for each simulation
            ActualChoiceTotalEvidence = 2*resp_posterior{s}-1==ChoiceTotalEvidence';
            pTotalEvidenceChoiceOnDisagreeAll = mean(ActualChoiceTotalEvidence(:,disagree_trials),2);

            % extract 95% confidence interval
            pTotalEvidenceChoiceOnDisagreeCI(:,s) = quantile(pTotalEvidenceChoiceOnDisagreeAll, [.025 .975]);

            % mean value
            pTotalEvidenceChoiceOnDisagree(s) = mean(pTotalEvidenceChoiceOnDisagreeAll);
        end
    end

    %% save
    modelfile = fullfile(ModelSaveDir, an.short,[animal '_' an.short]);
    save(modelfile,'nSubject','nDisagreeTrials','pTotalEvidenceChoiceOnDisagree*','subject_label');
    fprintf('saved results for %s for all %ss in %s\n',an.name, animal, modelfile);
end


%% Build regressor for stimulus mapping
function R =  mapping_regressor(stim, animal)
if strcmp(animal, 'monkey')
    vartype = 'categorical'; % pulses take discrete values
    binning = [];
    scale = 1:max(abs(stim));
    tau = [];
else
    vartype = 'continuous'; % stim take continuous values ( cos(angle) or auditory relative power)
    binning = .01;
    scale = [];
    tau = .1;
end
doOneHotEncoding = 1;

R = regressor(abs(stim), vartype,... % we fit a symmetric function so f(x) = sign(x)*f(|x|)
    'OneHotEncoding',doOneHotEncoding,'scale',scale,'constraint','fm',...
    'tau',tau, 'binning',binning,'label',{'temporal kernel','stimulus mapping'});

sgn = sign(stim);
sgn(isnan(stim)) = 0;
R = R * sgn;
end

%% RUN INTEGRATION MAP ANALYSIS
function M = run_integration_map(IntegrationModelFile,T,M, subject_id, bnd, nBootstrap, dx,sigma, FitConditionalPC)
if nargin<9
    FitConditionalPC=true;
end

nSubject = length(subject_id);

if ~isfile(IntegrationModelFile)
    warning('could not compute early and late evidence, missing file %s',  IntegrationModelFile);
    return;
end

% load file for integration model analysis where early and late
% evidence is stored
F = load(IntegrationModelFile);

%% loop through all subjects
for s=1:nSubject

    % use model data to define two columns arrays with early
    % and late evidence
    if isfield(F, 'M') % if GUM (i.e. integration model)
        earlylate =  [F.M.Predictions(s).early F.M.Predictions(s).late];
    else
        earlylate = [F.S_all(s).early F.S_all(s).late];
    end

    %probability of rightward response according to model for each
    %trial
    if isa(M, 'gum')
        Ymodel = M.Predictions(s).Expected;
    elseif isfield(M, 'Y')
        Ymodel = M(s).Y;
    else
        Ymodel = M(s).pModel;
    end

    % experimental responses for this subject

    % trial corresponding to this subject
    mask = T.subject== subject_id(s);
    Y = T.resp(mask);

    %% run integration map for experimental data and simulations, compute correlation
    S = integration_map_model_vs_data(earlylate,Y,Ymodel, dx, sigma, bnd, nBootstrap);

    %% FIT CONDITIONAL PSYCHOMETRIC CURVES
   if FitConditionalPC
    LateValue = -bnd:.5:bnd;

    S.ConditionalPsychometricCurvesData = fit_conditional_psychometric_curve(earlylate,Y, LateValue);
    S.ConditionalPsychometricCurvesModel = fit_conditional_psychometric_curve(earlylate,Ymodel, LateValue);
    S.LateEvidenceValue = LateValue;
   end

    %% ADD STRUCTURE TO OUTPUT
    if isa(M, 'gum')
        M.score.IntegrationMap(s) = S;
    else
        M(s).IntegrationMap = S;
    end
end
end

