function varargout = mono_tool(tool, varargin)
% tool subfunctions for analysis of Yates et al monkey data
%'datadir' % data directory
%'stimdir' % data directory
%'figdir' % main figure directory
% 'getsessions' % get filenames for all behavioral sessions
% 'loadsessions' % load behavioral data as a structure
% 'session_subplot' % create subplot for each sessionpar
% 'parcluster' % define parfor clusters
% 'build_summation_designmatrix' build design matrix for summation of pulses
% 'spatialadaptationmatrix' build design matrix for spatial adaptation
%  'initialize_component_from' initialize component (hyperparametes and weights from other fitted model)
%
%PLOT
% 'plot_session_bias'  plot session bias/sensitivity
% 'plot_temporal' plot temporal weights
%'plot_summation'  plot spatial summation
% 'plot_spatial' plot spatial weights
% 'plot_spatial_distance_angle'  plot spatial weights as function of distance and angle
% 'plot_history' plot history effect
% 'plot_temporal_adaptation' plot temporal adaptation
%  'plot_spatial_adaptation'   plot spatial adaptation
% 'plotmixture' plot parameters from mixture model
% 'plot_snapshot' plot parameters from snapshot model
% 'plot_responsibilities' plot responsabilities from mixture model
% 'savefig' save figures created by analysis
% 'earlylateintegration' to test for integration vs snapshot model, compute proportion of choices as a function of early and late evidence
% 'compute_earlylateintegration'  :  compute and plot map plot response vs early and late evidence



switch tool
    
    
    case 'datadir' % data directory
        list_of_dirs = {'C:\Users\Alex\Documents\Datasets\mtlipglm','/Users/alexhyafil/Documents/Databases/mtlipglm/',...
           'C:\Users\alexa\Documents\databases\Yates','D:\databases\mtlipglm\'};
        a=1;
        nogood = 1;
        % identify which of these directories actually exist on this
        % computer and have write access (will actually stop with the first one encountered form
        % the list)
        while nogood
            assert(a<=length(list_of_dirs), 'no directory for monkey data found on your computer, please update ''list_of_dirs'' on mono_tool.m line 36');
            nogood =  ~isdir(list_of_dirs{a});
            if nogood
                a = a+1;
            else % file exist let's see if it writable
                fileName = fullfile(list_of_dirs{a}, 'tmpFile.txt');
                [fid,errmsg] = fopen(fileName, 'w');
                if ~isempty(errmsg)&&strcmp(errmsg,'Permission denied') % no writing permission
                    nogood = 1;
                    a = a+1;
                else % we found our directory!
                    fclose(fid);
                    delete(fileName);
                end
            end
            
        end
        if a>length(list_of_dirs)
            error('could not identify an existing home directory for monkey data'),
        end
        dirpath =  list_of_dirs{a}; % got it!
        
        
        %         if ispc
        %             dirpath = 'D:\mtlipglm\'; %laptop
        %             if ~isdir(dirpath)
        %                 dirpath =  'C:\Users\U109469\Documents/Datasets/Yates/';
        %             end
        %         elseif ismac
        %             dirpath = '/Users/alexhyafil/Documents/Databases/mtlipglm/';
        %         else % cluster
        %             dirpath = '/tigress/ahyafil/mtlipglm/';
        %         end
        if nargin>1 && strcmp(varargin{1},'purebehav')
            dirpath = fullfile(dirroot,'full_behavioral_dataset');
        elseif nargin==1 || ~strcmp(varargin{1},'root')
            dirpath = fullfile(dirpath,'mtlipglm_data_full');
        end
        varargout = {dirpath};
        
    case 'stimdir' % data directory
        dirpath = mono_tool('datadir','root');
        if nargin>1 && strcmp(varargin{1},'purebehav') % directory for purely behavioral sessions
            dirpath = fullfile(dirpath,'full_behavioral_dataset');
        else
            dirpath = fullfile(dirpath,'mtlipglm_data_full','stim');
        end
        varargout = {dirpath};
        
    case 'figdir' % main figure directory
        dirpath = [dropboxroot 'mono/figures/'];
        varargout = {dirpath};
        
    case 'getsessions' % get filenames for all behavioral sessions
        % all session files
        purebehavdataset = nargin>1 && strcmp(varargin{1},'purebehav'); %full dataset
        fulldataset = nargin>1 && strcmp(varargin{1},'full'); %full dataset
        if fulldataset
            [allfilename1, monkey1, dates1, session_nb1] = mono_tool('getsessions'); % sessions with electrophysiology
            [allfilename2, monkey2, dates2, session_nb2] = mono_tool('getsessions','purebehav'); % session without electrophysiology
            allfilename = [allfilename1 allfilename2];
            monkey = [monkey1 monkey2];
            dates = [dates1 dates2];
            session_nb = [session_nb1 session_nb2];
            
        else
            
            thisdir = mono_tool('stimdir',varargin{:});
            allfile = dir(fullfile(thisdir,'*.mat'));
            allfilename = {allfile.name};
            assert(~isempty(allfilename), 'could not find any file ... check if directory is set up correctly');
            
            % select files matching any of the patterns
            if purebehavdataset
                expr = '([a-z]{3})(\d{8})-\d{4}.mat';
            else
                expr = '([\w])(\d*)([a-z]?)_stim.mat';
            end
            
            reg = regexp(allfilename, expr,'tokens');
            if any(cellfun(@isempty, reg))
                error('filename not matching regular expression');
            end
            reg = [reg{:}];
            reg = cat(1,reg{:});
            
            if purebehavdataset
                monkey = reg(:,1)'; % monkey letter (n/p)
                sess_incl = strcmp(monkey, 'nan') | strcmp(monkey, 'pat'); % only take sessions from monkeys
                monkey = cellfun(@(x) x(1), monkey, 'uniformoutput',0); % just retain first letter
                dates = datenum(reg(:,2)','yyyymmdd')'; % convert date string to number
                session_nb = zeros(1,length(allfilename)); % good get it from time byt too lazy
                
                % set only 'flat' experiment
                for f=1:length(allfilename)
                    variableInfo = who('-file',fullfile(thisdir,allfilename{f}));
                    if ismember('temporalWeighting',variableInfo) % if variable 'temporalWeighting' exists
                        load(fullfile(thisdir,allfilename{f}),'temporalWeighting'); % load variable from dataset
                        sess_incl(f) = sess_incl(f) && strcmp(temporalWeighting, 'flat');
                    end
                end
                
                % select only required sessions
                allfilename = allfilename(sess_incl);
                monkey = monkey(sess_incl);
                dates = dates(sess_incl);
                session_nb = session_nb(sess_incl);
                
            else % ephys sessions
                monkey = reg(:,1)'; % monkey letter (n/p)
                dates = datenum(reg(:,2)','yyyymmdd')'; % convert date string to number
                reg(cellfun(@isempty,reg(:,3)'),3) = {'a'}; % if no letter identifying session within day, means only one in that day
                session_nb = cellfun(@(x) double(x)-96, reg(:,3)'); % convert letter a/b/c to session letter
            end
            
            % add directory to name
            for f=1:length(allfilename)
                allfilename{f} = fullfile(thisdir, allfilename{f});
            end
        end
        
        % exclude some sessions (double count between recording sessions
        % and behavioral sessions; session 'p20140213' has no seedsmatch)
        excl_sess = {'p20140213_stim','nan20150305-1308','nan20150305-1602','nan20150306-1452','nan20150309-1351',...
            'nan20150316-1558','nan20150324-1228','nan20150325-1406','nan20150326-1309','nan20150331-1425',...
            'nan20150401-1343','nan20150407-1212','nan20150408-1312','nan20150408-1507',...
            'nan20150416-1130','nan20150501-1323','nan20150505-1129','nan20150609-1458','pat20140303-1545',...
            'pat20140304-1407','pat20140305-1454','pat20140306-1341','pat20140307-1409','pat20140310-1453'};
        for f=1:length(allfilename)
            [~,sessname{f}] = fileparts(allfilename{f}); % extract filename
        end
        [~,excl_ind] = ismember(excl_sess,sessname); % find index for these sessions in list of sessions
        excl_ind = excl_ind(excl_ind>0); % in
        allfilename(excl_ind) = [];
        monkey(excl_ind) = [];
        dates(excl_ind) = [];
        session_nb(excl_ind) = [];
        
        % sort by monkey, date and session
        score = 10000*cellfun(@(x) strcmp(x,'p'), monkey) + dates + session_nb/1000;
        [~,order] = sort(score);
        allfilename = allfilename(order);
        monkey = monkey(order);
        dates = dates(order);
        session_nb = session_nb(order);
        
        varargout = {allfilename, monkey, dates, session_nb};
        
    case 'loadsessions' % load behavioral data as a structure
        [allfilename,monkey,dates,session_nb] = mono_tool('getsessions',varargin{:});
        nsessions = length(allfilename); % number of sessions
        
        fprintf('loading data from %d sessions ...',nsessions);
        
        %load and append data from each session
        S = [];
        for i=1:nsessions
            this_S = load(allfilename{i}); % load data from session
            %   this_S = load(fullfile(mono_tool('stimdir', varargin{:}),allfilename{i})); % load data from session
            if isfield(this_S, 'eyepos')
                this_S = rmfield(this_S, 'eyepos');
            end
            this_S.monkey = monkey{i};
            this_S.dates = dates(i);
            this_S.session_nb = session_nb(i);
            if ~isfield(this_S, 'pmat'), this_S.pmat = []; end
            if ~isfield(this_S, 'seedsmatch'), this_S.seedsmatch = []; end
            
            % if GaborXY data missing, compute again (some sessions for
            % 'nan')
            if all(isnan(this_S.gaborXY)) % replace sessions with nan gaborXY
                if ~strcmp(this_S.monkey,'n')
                    error ('gaborXY data missing, do not know which parameters to regenerate');
                end
                dv.pa.center = this_S.motionApertureCenterXY;
                dv.pa.v1rf = 0.13;
                dv.pa.mtrf = 0.8; % values from Leor email, 03/12/18
                dv.pa.theta = this_S.theta;
                dv = makeGaborPositionsFixed(dv);
                this_S.gaborXY = dv.pa.pos';
            end
            
            S = structcat(S, this_S);
            
        end
        
        fprintf('done\n');
        
        varargout = {S,nsessions};
        
    case 'session_subplot' % create subplot for each session
        sess = varargin{1}; % session index
        exname = varargin{2}; % structure
        switch sess
            case 0
                s_sess = 48;
                titl = 'average all';
            case -1
                s_sess = 47;
                titl = 'average P';
            case -2
                s_sess = 24;
                titl = 'average N';
            otherwise
                s_sess = sess+1*(sess>23); % for monkey N, switch to new line
                titl = exname{sess};
        end
        subplot(6,8,s_sess);
        title(titl);
        
    case   'parcluster' % define parfor clusters
        if strcmp(getenv('SLURM_CLUSTER_NAME'),'della')
            pc = parcluster('local');
            if nargin>1
                script = varargin{1};
            else
                script = '';
            end
            disp(fullfile('/home/',getenv('USER'),'junk',[script getenv('SLURM_ARRAY_TASK_ID')]));
            pc.JobStorageLocation = fullfile('/home/',getenv('USER'),'junk',[script getenv('SLURM_ARRAY_TASK_ID')]);
            
            if nargin>2
                NumWorkers = varargin{2};
            else
                NumWorkers = pc.NumWorkers; % default number of workers
            end
            parpool(pc,NumWorkers);
        end
        
        %% build design matrix for summation of pulses
    case 'build_summation_designmatrix'
        warning('should be outdated');
        stim = varargin{1};
        Xsummation = cell(1,3);
        Xsummation{1} = sign(stim); % sign of pulse
        ngab = abs(stim); % (unsigned) number of Gabor
        ngab_uq = setdiff(unique(ngab(:)),0); % unique values used (excluding zero, that has no associated weight)
        ngab_id = zeros(size(ngab)); % replace number of Gabor by its index
        for gg =1:length(ngab_uq)
            ngab_id(ngab==ngab_uq(gg)) = gg;
        end
        Xsummation{3} = ngab_id;
        
        varargout = {Xsummation, ngab_uq};
        
        %% build design matrix for spatial adaptation
    case 'spatialadaptationmatrix'
        S= varargin{1}; % data structure for all sessions
        s = varargin{2}; % which monkey
        dist_boundaries = varargin{3}; % bins for distance between Gabors
        nsession_monkey = varargin{4};
        nGabor_monkey = varargin{5}; % number of Gabor per monkeys
        across_list = varargin{6};
        ntrial_monkey = varargin{7};
        nPulses = S(1).nPulses;
        ndist_bins = length(dist_boundaries)-1;
        
        this_sess = find(across_list{s});
        nincltrial_session = cellfun(@sum, {S(this_sess).incltrial}); % number of good trials per session
        nGabor_cum =  cumsum([0 [S(this_sess).nGabors]]); % cumulated number of Gabors per session
        nincltrial_cum = cumsum([0 nincltrial_session]); % cumulated values
        XSA = zeros(ntrial_monkey(s),nPulses, nGabor_monkey(s),1+ndist_bins); % pre-allocate
        
        for i=1:nsession_monkey(s)
            % gabor_id = 1:nGabor_monkey(s);
            gg_idx = nGabor_cum(i)+1:nGabor_cum(i+1); %indices for Gabors in this session
            tt_idx = nincltrial_cum(i)+1:nincltrial_cum(i+1); % indices for trials in this session
            this_S = S(this_sess(i));
            incl_trial = this_S.incltrial; % trials to include in analysis
            this_X = this_S.pulses(incl_trial,:,:); % pulses for this session
            
            % design matrix for spatial adaptation
            Gabdist = dist(this_S.gaborXY'); % matrix of distance between each Gabors for this session
            
            % compute corresponding bin
            GDB = zeros(this_S.nGabors);
            for b=1:ndist_bins
                GDB = GDB + (Gabdist>dist_boundaries(b));
            end
            
            this_XA = zeros(nincltrial_session(i),nPulses, this_S.nGabors,1+ndist_bins);
            this_XA(:,:,:,1) = this_X; % first component: plain 3d design matrix
            for g=1:this_S.nGabors
                coactive = this_X(:,:,g) .* this_X; % coactive Gabors
                this_gdb = permute(GDB(g,:),[1 3 2]); % distance bins from this gabor
                for b=1:ndist_bins
                    this_XA(:,:,g,b+1) = this_X(:,:,g) .* sum(coactive .* (this_gdb==b),3);
                end
            end
            XSA(tt_idx,:, gg_idx,:) = this_XA; % into big matrix
        end
        
        varargout = {XSA};
        
        %% initialize component (hyperparametes and weights from other fitted model)
    case 'initialize_component_from'
        param = varargin{1}; % param model structure
        HP_ini  = varargin{2}; % initial values of hyper parameters
        c_from =   varargin{3}; % which component to take from
        c_to = varargin{4}; % which component to update
        initfromfile = varargin{5}; % which file to update from
        if exist(initfromfile,'file')  || exist([initfromfile '.mat'],'file')
            IF = load(initfromfile);
            if length(HP_ini)<max(c_to) || length(IF.HP) < max(c_from)
                varargout = {param,HP_ini};
                fprintf('Incorrect number of components from model file, initializing from default values (file:''%s'')\n', initfromfile);
                return;
            end
            siz = cellfun(@length, HP_ini(c_to));
            siz2 = cellfun(@length, IF.HP(c_from));
            excl = siz ~=siz2 & (siz>0);
            if any(excl)
                fprintf('Number of %d component hyperparameters from model file does not correspond, initializing from default values (file:''%s'')\n', sum(excl), initfromfile);
                c_to(excl) = [];
                c_from(excl) = [];
            end
            if isfield(param,'U')
                siz = cellfun(@length, param.U(c_to));
                siz2 = cellfun(@length, IF.U(c_from));
                excl = (siz ~=siz2) & (siz>0); % if number of weights do not correspond
                if any(excl)
                    fprintf('Number of %d component weights from model file does not correspond, initializing from default values (file:''%s'')\n', sum(excl), initfromfile);
                    c_to(excl) = [];
                    c_from(excl) = [];
                end
            end
            param.U(c_to) = IF.U(c_from); % initial value of weights
            HP_ini(c_to)  = IF.HP(c_from); % initial value of hyperparameters
        else
            fprintf('Initializing model file not found, initializing from default values (file:''%s'')\n',initfromfile);
        end
        varargout = {param,HP_ini};
        
        %% %%% PLOTTING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% plot session bias/sensitivity
    case 'plot_session_bias'
        k_sess_all = varargin{1}; % session weight
        subdir = varargin{2}; % analysis subdirectory
        subject_label = varargin{3};
        if nargin>4
            titl = varargin{4};
        else
            titl = 'session bias';
        end
        nSubject = length(k_sess_all);
        subfigure([subdir titl]);
        for s=1:nSubject
            subplot(1,nSubject,s); hold on; title(subject_label{s});
            plot([1 size(k_sess_all{s},2)],[0 0], 'color',.8*[1 1 1]);
            wu(k_sess_all{s,1}', k_sess_all{s,2}','color','b','errorstyle','fill');
            xlabel('session');  axis tight;
            if s==1, ylabel('weight'); end
        end
        sameaxis;
        
        %% plot temporal weights
    case 'plot_temporal'
        k_temp_all = varargin{1}; % weights
        subdir = varargin{2}; % analysis subdirectory
        subject_label = varargin{3};
        
        nSubject = length(subject_label);
        subfigure([subdir 'temporal kernel']);
        nFrame = size(k_temp_all,2);
        for s=1:nSubject
            subplot(1,nSubject,s); hold on; title(subject_label{s});
            plot([1 nFrame],[0 0], 'color',.8*[1 1 1]);
            wu(k_temp_all(:,s,1), k_temp_all(:,s,2),'color','b');
            xlabel('pulse \it t'); ylabel('weight \it w_t'); axis tight;
        end
        
        %% plot spatial summation
    case 'plot_summation'
        f_summation_all = varargin{1}; % weights
        subdir = varargin{2}; % analysis subdirectory
        subject_label = varargin{3};
        ngab_uq = varargin{4};
        nSubject = length(subject_label);
        lincoeff = zeros(1,nSubject);
        subfigure([subdir 'evidence summation']);
        for s=1:nSubject
            subplot(1,nSubject,s); hold on;
            wu(ngab_uq{s}, f_summation_all{s,1}',f_summation_all{s,2}','color','b','errorstyle','fill');
            lincoeff(s) = f_summation_all{s,1}(1) / ngab_uq{s}(1); % linear coefficient computed from first point
            plot(ngab_uq{s}, lincoeff(s)*ngab_uq{s}, 'color',.7*[1 1 1]);
            xlabel('#pulse'); ylabel('evidence');
        end
        varargout = {lincoeff};
        
        %% plot spatial weights
    case 'plot_spatial'
        allgaborXY = varargin{1}; % cell with position of all Gabors
        k_spa_all = varargin{2}; % weights
        subdir = varargin{3}; % analysis subdirectory
        subject_label = varargin{4};
        if nargin>5
            titl = varargin{5};
        else titl = 'spatial kernel';
        end
        
        nSubject = length(allgaborXY);
        subfigure([subdir titl]);
        clim = [min(cat(2,k_spa_all{:,1})) max(cat(2,k_spa_all{:,1}))];
        for s=1:nSubject
            subplot(1,nSubject,s); hold on; title(subject_label{s});
            plot([0 0], [min(allgaborXY{s}(:,2)) max(allgaborXY{s}(:,2))],'color',.7*[1 1 1]);
            plot([min(allgaborXY{s}(:,1)) max(allgaborXY{s}(:,1))],[0 0], 'color',.7*[1 1 1]);
            colorplot(allgaborXY{s}(:,1), allgaborXY{s}(:,2), k_spa_all{s,1}, 'cmap',jet, 'clim',clim);
            axis tight;
        end
        
        %% plot spatial weights as function of distance and angle
    case 'plot_spatial_distance_angle'
        allgaborXY = varargin{1}; % cell with position of all Gabors
        k_spa_all = varargin{2}; % weights
        subdir = varargin{3}; % analysis subdirectory
        subject_label = varargin{4};
        if nargin>5
            titl = varargin{5};
        else titl = 'spatial distance angle';
        end
        
        nSubject = length(allgaborXY);
        subfigure([subdir titl]);
        for s=1:nSubject
            allgaborZ = allgaborXY{s}(:,1)+1i*allgaborXY{s}(:,2); % convert to complex representation
            subplot(2,nSubject,s); hold on; title(['distance ' subject_label{s}]);
            grpmeann(k_spa_all{s}, abs(allgaborZ), 'quantiles',20,'plot','ste');
            subplot(2,nSubject,s+nSubject); title(['angle ' subject_label{s}]);
            grpmeann(k_spa_all{s}, angle(allgaborZ), 'histc',pi*(-1:.25:1),'plot','ste');
            set(gca, 'xtick',pi*(-1:.5:1), 'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});
        end
        
        %% plot history effect
    case 'plot_history'
        k_history_all = varargin{1}; % weights for previous trials
        subdir = varargin{2}; % analysis subdirectory
        subject_label = varargin{3};
        
        nSubject = length(subject_label);
        subfigure([subdir 'history']);
        for s=1:nSubject
            subplot(1,nSubject,s); hold on; title(subject_label{s});
            plot([1 7],[0 0], 'color',.8*[1 1 1]);
            wu(squeeze(k_history_all(s,:,:,1)), squeeze(k_history_all(s,:,:,2)),{{},{'correct','error'}});
            xlabel('trial lag'); ylabel('weight'); axis tight;
        end
        sameaxis;
        
        %% plot temporal adaptation
    case 'plot_temporal_adaptation'
        k_adapt_all = varargin{1};
        subdir = varargin{2};
        subject_label = varargin{2};
        nSubject = length(subject_label);
        maxlag = size(k_adapt_all, 1)-1; % number of temporal lags
        subfigure([subdir 'temporal adaptation']);
        for s=1:nSubject
            subplot(1,nSubject,s); hold on; title(subject_label{s});
            plot([1 7],[0 0], 'color',.8*[1 1 1]);
            wu(squeeze(k_adapt_all(s,:,:,1)), squeeze(k_adapt_all(s,:,:,2)),...
                'line',{{},{'same','opposite'}});
            xlabel('pulse'); ylabel('weight'); axis tight; xlim([.5 maxlag+.5]);
        end
        sameaxis;
        
        %% plot spatial adaptation
    case 'plot_spatial_adaptation'
        k_adapt_all = varargin{1}; % weights for previous trials
        dist_boundaries = varargin{2}; % boundaries between distance bins.
        subdir = varargin{3}; % analysis subdirectory
        subject_label = varargin{4};
        nSubject = length(subject_label);
        
        % plot adaptation effect
        subfigure([subdir 'spatial adaptation']);
        for s=1:nSubject
            subplot(1,nSubject,s); hold on; title(subject_label{s});
            plot([1 length(dist_boundaries)-1],[0 0], 'color',.8*[1 1 1]);
            wu(dist_boundaries(1:end-1), k_adapt_all(s,:,1), k_adapt_all(s,:,2));
            xlabel('distance'); ylabel('modulation'); axis tight;
        end
        sameaxis;
        
        %% plot parameters from mixture model
    case 'plotmixture'
        k_temp_all = varargin{1};
        pi_all = varargin{2};
        subject_label = varargin{3};
        with_components = (nargin>4); % each components (not plotted for snapshot model)
        if with_components
            components = varargin{4};
        end
        nSubject = length(subject_label);
        for s=1:nSubject
            subplot(2+with_components,nSubject,s); hold on; title(subject_label{s});
            plot([1 7],[0 0], 'color',.8*[1 1 1]);
            wu(k_temp_all(:,:,s,1)', k_temp_all(:,:,s,2)',{{},{'bias','sensitivity'}});
            xlabel('component'); ylabel('weights');
            
            subplot(2+with_components,nSubject,s+2);
            wu(pi_all(:,s,1),pi_all(:,s,2));
            ylabel('mixture coeffs');
            
            if with_components
                subplot(3,nSubject,s+4); hold on; title(subject_label{s});
                plot([1 7],[0 0], 'color',.8*[1 1 1]);
                wu(components(:,:), []);
                ylabel('components');
            end
        end
        
        %% plot parameters from snapshot model
    case 'plot_snapshot'
        k_temp_all = varargin{1};
        pi_all = varargin{2};
        subdir = varargin{3};
        animal_label = varargin{4};
        nSubject = length(animal_label);
        is_deterministic = isempty(k_temp_all);
        subfigure([subdir 'logistic mixture']);
        nLine = 1+~is_deterministic;
        for s=1:nSubject
            
            
            subplot(nLine,nSubject,s);
            wu(pi_all(:,s,1),pi_all(:,s,2));
            if s==1, ylabel('mixture coeffs \pi_t'); end
            xlabel('pulse');
            
            if ~is_deterministic
                subplot(nLine,nSubject,s+nSubject); hold on; title(animal_label{s});
                plot([1 7],[0 0], 'color',.8*[1 1 1]);
               % wu(k_temp_all(:,:,s,1)', k_temp_all(:,:,s,2)',{{},{'\it \beta_t','\it w_t'}});
               wu(k_temp_all(:,:,s,1)', k_temp_all(:,:,s,2)');
                if s==1, ylabel('weights');
                else legend off;
                end
            end
        end
        
        %% plot responsabilities from mixture model
    case 'plot_responsibilities'
        S_all = varargin{1};
        session_id = varargin{2};
        subdir = varargin{3};
        animal_label = varargin{4};
        
        nSubject = length(S_all);
        subfigure([subdir 'responsibilities']);
        for s=1:nSubject
            subplot(1,nSubject,s); title(animal_label{s});
            grpmeann(S_all(s).responsibility,session_id{s}', 'dim',1,'plot','ste','color','jet');
            xlabel('session');
            if s==1
                ylabel('responsibilities');
            else
                legend off;
            end
        end
        
        
        %% save figures created by analysis
    case 'savefig'
        if nargin>1
            savefigdir = varargin{1}; % directory for saving provided as argument
            if ~isdir(fileparts(savefigdir)) % create subdirectory if does not exist already
                mkdir(fileparts(savefigdir));
                fprintf('created directory ''%s''\n',fileparts(savefigdir));
            end
        else
            savefigdir = pwd;
        end
        allfig = get(0,'child');
        fprintf('saving %d figures to %s:\n', length(allfig), savefigdir);
        for f=1:length(allfig)
            ff = allfig(f);
            
            ffname = fullfile(savefigdir,get(ff,'filename'));
            if ~isdir(fileparts(ffname)) % create sub-subdirectory if does not exist already
                mkdir(fileparts(ffname));
                disp(['created ''' fileparts(ffname) ''' directory']);
            end
            saveas(ff,ffname,'fig'); % save as FIG
            fprintf([get(ff,'filename') '\n']);
        end
        fprintf('done\n');
        
        %% to test for integration vs snapshot model, compute proportion of choices as a function of early and late evidence
    case 'earlylateintegration'
        
        XXX = varargin{1}; % design matrix (or model file with already early and late)
        YYY = varargin{2};
        S_all = varargin{3}; % model structure
        subject_label = varargin{4};
        subdir = varargin{5};
        bnd = varargin{6}; % max value
        if length(varargin)>6
            Nbtstrp = varargin{7};
        else
            Nbtstrp = 0;
        end
        nSubject = length(subject_label);
        
        withmodelfile = ischar(XXX);
        if withmodelfile
            if ~isfile(XXX)
                warning('could not compute early and late evidence, missing file %s',  XXX);
                varargout = {S_all};
                return;
            end
            F = load(XXX);
            dx = .2; % bin size
            % bnd = 2.5; % max value
        else
            dx = 2;
            bnd = 30;
        end
        
        subfigure([subdir 'earlylateintegration']);
        for s=1:nSubject
            if withmodelfile
                if isfield(F, 'M') % GUM
                  earlylate =  [F.M.Predictions(s).early F.M.Predictions(s).late];
                else
                earlylate = [F.S_all(s).early F.S_all(s).late];
                end
            else % simply sum first three samples and last 4
                stim = sum(XXX{s},3); % net overall number of pulses
                earlylate = [sum(stim(:,1:3),2) sum(stim(:,4:7),2)]; % sum evidence from first 3 samples and last 4 samples
            end
            
            subplot(4,nSubject,s);
            if isa(S_all, 'gum')
                Ymodel = S_all.Predictions(s).Expected;
            elseif isfield(S_all, 'Y')
                Ymodel = S_all(s).Y;
            else
                Ymodel = S_all(s).lh;
            end
            [EL.mean, EL.se, ndatapoints] = mono_tool('compute_earlylateintegration', earlylate, Ymodel, dx, bnd,1);
            title(['model ' subject_label{s}]);
            
            subplot(4,nSubject,s+nSubject);
            [EL.mean_data, ELdata_se] = mono_tool('compute_earlylateintegration', earlylate, YYY{s},dx,bnd,1);
            title(['data ' subject_label{s}]);
            
            doplot = (ndatapoints>20);
            
            %% compute correlation and SSE w.r.t data
            EL.r = corr(EL.mean(doplot), EL.mean_data(doplot)); % Pearson's correlation coeff (only including bins with sufficient datapoints)
            EL.sse = sum(sum((EL.mean-EL.mean_data).^2)); % sum of square difference
            
            EL.dx = dx;
            EL.bnd = bnd;
            
            %% add boostrapped value
            if Nbtstrp>0
               fprintf('computed boostrap values for correlation with experimental data...');
               EL.r_btstrp = zeros(1,Nbtstrp);
               EL.sse_btstrp = zeros(1,Nbtstrp);
               nTrial = length(Ymodel); % number of trials
               for b=1:Nbtstrp
                    rd = randi(nTrial,1,nTrial); % random draw of trials (with replacement)
                    [mean_bs, ~, ndatapoints_bs] = mono_tool('compute_earlylateintegration', earlylate(rd,:), Ymodel(rd), dx, bnd,0);            
                    mean_data_bs = mono_tool('compute_earlylateintegration', earlylate(rd,:), YYY{s}(rd),dx,bnd,0);
            
                 doplot_bs = (ndatapoints_bs>20);
            
            % compute correlation and SSE w.r.t data
            EL.r_btstrp(b) = corr(mean_bs(doplot),mean_data_bs(doplot)); % Pearson's correlation coeff (only including bins with sufficient datapoints)
            EL.sse_btstrp(b) = sum(sum((mean_bs-mean_data_bs).^2)); % sum of square difference
                                                       
                fprintf('*');
               end
               fprintf('done\n'); 
            end
            
            %% plot psychometric curve for early evidence conditioned on late evidence
            colz_pc = [0 0 0; 1 0 0]; % color of psychometric curves (graded from black to red)

            doplot = double(doplot);
            doplot(doplot==0) = nan; % do not plot when less than 20 data points
            doplot = doplot(2:end-1,:);
            
            xval = -bnd:dx:bnd;
            xval = xval(1:end-1) + dx/2;
            x2_val_max = min(2,bnd-dx/2);
            x2val = - x2_val_max:x2_val_max; % values of late evidence we condition on
            ival = find(xval==-x2_val_max,1)+ (1:1/dx:2*x2_val_max/dx+1);
            
            
            nPC = length(ival);
            colz_pc_all = colz_pc(1,:) + linspace(0,1,nPC)'*diff(colz_pc); % graded colours for psychometric curves 

            % ival = 4:1/dx:2*bnd/dx;
            lateevidence_label = num2strcell('late ev =%d',x2val);
            subplot(4,nSubject,s+2*nSubject);
            
            wu(xval, doplot(:,ival).*EL.mean(2:end-1,ival),EL.se(2:end-1,ival),'color',colz_pc_all,{{},lateevidence_label});
            if s>1, legend off; end
            xlabel('early evidence'); ylabel('pright');
            title(['model ' subject_label{s}]); axis tight;
            
            subplot(4,nSubject,s+3*nSubject);
            wu(xval,doplot(:,ival).*EL.mean_data(2:end-1,ival),ELdata_se(2:end-1,ival),'color',colz_pc_all,{{},lateevidence_label});
            legend off;
            xlabel('early evidence'); ylabel('pright');        title(['data ' subject_label{s}]); axis tight;
            
            %% find bias and lapse as function of evidence
            lateval = -bnd:.5:bnd;
            PLoptions = struct('ninit',1,'alpha',[1 1 3],'maxiter',30000);
            for i=1:length(lateval)-1
                this_trial = earlylate(:,2)>=lateval(i) & earlylate(:,2)<lateval(i+1);
                % offset = [0 mean(lateval([i i+1])];
                PLdata = probitlapse(earlylate(this_trial,1),YYY{s}(this_trial),0,2,[],PLoptions);
                EL.joint_data(:,:,i) = PLdata.joint;
                PLmodel = probitlapse(earlylate(this_trial,1),Ymodel(this_trial),0,2,[],PLoptions);
                EL.joint_model(:,:,i) = PLmodel.joint;
            end
            EL.xlateevidence = lateval;
            
            if isa(S_all, 'gum')
                S_all.score.earlylate(s) = EL;
            else
            S_all(s).earlylate = EL;
            end
        end
        varargout = {S_all};
        
        %% compute and plot map plot response vs early and late evidence
    case 'compute_earlylateintegration'
        X = varargin{1}; % ntrial x 2 matrix (early evidence in first col, late evidence in second)
        Y = varargin{2}; % choices
        dx = varargin{3}; % binsize
        bnd = varargin{4}; % bound
        do_plot = varargin{5}; % whether to plot or not
        
        colz = [103 169 221; 241 135 34]/255; % gradient between red and blue
        
        h1 = -bnd:dx:bnd;
        hh = [-Inf h1 Inf]; % bin boundaries
        xx = [h1(1)-dx/2 h1+dx/2];
        %  [M,M_se] = grpmeann(Y, {X(:,1) X(:,2)}, 'dim',1,'histc',{hh,hh},'ste');
        
        % count number of trials in each bin
        S = grpsum(ones(size(X,1),1), {X(:,1) X(:,2)}, 'dim',1,'histc',{hh,hh});
        
        % sum of rightwards responses in each bin
        M = grpsum(Y, {X(:,1) X(:,2)}, 'dim',1,'histc',{hh,hh});
        
        % smooth by convolution with gaussian kernel
        nConv = 10; % length of convolution kernel (in bins)
        sigma = .5; % width of kernel (in bins)
        %conv_vec = cos(pi/2*linspace(-1,1,nConv+1));
        conv_vec = exp(-linspace(-1,1,nConv+1).^2/2/sigma^2);
        M = conv2(M,conv_vec'*conv_vec,'same');
        S = conv2(S,conv_vec'*conv_vec,'same');
        M = M./S; % compute mean (sum/number of trials)
        varargout = {M,nan(size(M)),S}; % output mean and standard error
        
        
        %% plot if required
        if do_plot
        int = min(S/20,1); % hue intensity
        
        M(isnan(M)) = 0;
        % convert mean response and number of trials to RGB
        %C = cat(3, M.*int,int,(1-M).*int);
        for c=1:3
            C(:,:,c) = 1- (colz(1,c)-colz(2,c)) * M - colz(2,c);
        end
        C = 1- C.* int;
        % image(1-permute(C,[2 1 3])); axis xy; hold on;
        image(xx, xx, permute(C,[2 1 3])); axis xy; hold on;
        contour(xx, xx, M', [.15 .3 .7 .85],'k', 'linewidth',.5); % draw contour lines
        contour(xx, xx, M',[.5 .5],'k', 'LineWidth',1); % draw contour lines
        xlabel('early evidence'); ylabel('late evidence');
      %  set(gca, 'xtick',[1.5 length(hh)/2 length(hh)-1.5], 'xticklabel', num2strcell([hh(2) 0 hh(end-1)]));
      %  set(gca, 'ytick',[1.5 length(hh)/2 length(hh)-1.5], 'yticklabel', num2strcell([hh(2) 0 hh(end-1)]));
        end
        
        %% wrong analysis
    otherwise
        error('unknown option %s', tool);
end

