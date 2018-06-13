% function rivalry_Batch
% Template batch script for SPM5 for the fixed effects analysis
% includes realignment parameters
% modified abartels 2008.10.
clear all;
do_contrasts_only = 0; % ONLY CALCULATE CONTRASTS
stats_type = 1;
globalregress_flag = 1;
globalnorm_flag = 0;

spm_read_shared_parameters
%	1: 6 conditions: percepts A and B and switches (pooled), separately for rivalry and for replay.
%	2: 4 conditions: percepts A and B, separately for rivalry and for replay.
%	3: 2 conditions: switches (pooled), separately for rivalry and for replay.


% ------------------------------------------------------------------------
% -------------- SUBJECT SPECIFIC INFO -----------------------------------------

% main_dir = '/Users/soyoung/Atten_Corr_BOLD/';
% task_dir = {
%     '10_Vis_Short/', '11_Motion_OnOff/', '12_Vis_Sess/', ...
% %     '/15_Aud_Local/', '/16_Aud_ROI/', '/17_Aud_Sess/'
%     };
% ------------------ END SUBJECT SPECIFIC INFO ------------------
n_subjects = length(subj);

% ----------------- Analysis params ---------------------------------
high_pass_sec = 256; % sec high-pass cut-off
n_dummies = 0;

% Define numbers of sessions/conditions
cond_n = 2; %conditions

% Results dir: create the names for the directories of the 'fixed effects analysis'
dir_analysis = 'sraf_reg';

file_img_prefix = 'sraf'; % identifier for img files to be analyzed.



img_diff_flag = 0;	% include frame-wise pixel difference as regressor.
realign_flag = 1;	% include motion realign params as regressors


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%SPM DESIGN SPECIFICATION & ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
% definition of experimental conditions (i.e. regressors in design matrix)
%==========================================================================
%% Conditions:
%==========================================================================

% define condition names, i.e. see above.
cond_names ={'attention', 'resting'};
% cname{1}.one = 'attention vs resting';
% cname{1}.two = 'resting vs attention';
% cname{2}.one = 'motion vs static';
% cname{2}.two = 'static vs motion';
% cname{3}.one = 'attention vs resting';
% cname{3}.two = 'resting vs attention';
% cond_names ={'ON', 'OFF'};
%==========================================================================
% load subject-specific onsets times


for thesubject = 3%1:n_subjects% [1:n_subjects],%[4 5 7 8]
    sn=thesubject; %(Subject Number)
    for task = [ 1 ]
        n_sess = length(subj(sn).func_dir{task}); % length(subj(sn).sess_dir_nums);
        subj(sn).dir_ana_parent = [main_dir, task_dir{task}, subj_dir{sn}];
        n_scans = subj(sn).n_scans{task};%repmat(subj(sn).nscans{task}, 1, size(subj(sn).func_dir{task},2));
        n_ses = size(subj(sn).func_dir{task}, 2);
        
        clear SPM	% clear old stuff
        spm fmri	% open spm
        
        % ----------- assign subject specific
        % take the normal defaults
        global defaults
        
        if ~do_contrasts_only,	% do stats + contrasts
            
            
            % make the directory
            mkdir([subj(sn).dir_ana_parent, '/Stat/'],dir_analysis);
            
            
            
            % create cell array cond_ons in which all conditions are stored
            % onset time: 0 sec is the middle of first scan, i.e. TR/2.
            cond_ons    = cell(n_sess,cond_n);
            cond_dur    = cell(n_sess,cond_n);
            
            % ---------- EDIT HERE CONDITOIN ONSESTS ------ Session-specific conditions:
            for ses = 1:n_ses%{sn}
                psychname = subj(sn).fname{task};
                if task == 3, loadname = [psychname, num2str(ses)]; else, loadname = psychname; end
                load(loadname, 'attention_switch_t');%,num2str(ses)], 'attention_switch_t');
                cond_ons{ses,1} = attention_switch_t(attention_switch_t(:,1) ==1,2)'; % attention onset
                cond_ons{ses,2} = attention_switch_t(attention_switch_t(:,1) ==0,2)'; % resting onset
                cond_ons{ses,3} = [cond_ons{ses,1} cond_ons{ses,2}]; % onset of both
                
                cond_rep = size(cond_ons{ses,1},2);
                cond_duration = round(mean(diff(attention_switch_t(:,2))));
                cond_dur{ses,1} = [repmat(cond_duration, 1, cond_rep)]; % attention dur
                cond_dur{ses,2} = [repmat(cond_duration, 1, cond_rep)];% round(mean(diff(attention_switch_t(:,2))));%[repmat(120, 1, 2)]; % resting dur
                cond_dur{ses,3} = [cond_dur{ses,1} cond_dur{ses,2}]; % both dur
            end
            
            
            %==========================================================================
            % number of scans and session, e.g. [128 128 128] for 3 sessions and TR
            SPM.nscan          = n_scans;%{sn}; % n_scans*ones(1,n_sess);
            SPM.xY.RT          = subj(sn).TR; % sec AUSRECHNEN MIT MATFILE! mean(diff(timestamps ????))   hist(diff(..))
            
            %==========================================================================
            % basis functions and timing parameters (#mod#)
            SPM.xBF.name       = 'hrf'; %'hrf (with time derivative)';
            SPM.xBF.order      = 1;% 1: hrf only; 2: with time derivative; 3: dispersion der.
            SPM.xBF.length     = 32;% length in seconds
            SPM.xBF.T          = 16;% number of time bins per scan
            SPM.xBF.T0         = 1; % first slice/timebin (if slicetiming)
            SPM.xBF.UNITS      = 'secs';           % OPTIONS: 'scans'|'secs' for onsets
            SPM.xBF.Volterra   = 1;                % OPTIONS: 1|2 = order of convolution
            
            %==========================================================================
            % Trial specification: onsets, duration (in sec) and parameters for mod.
            
            for sess=1:n_sess,
                for cond=1:cond_n
                    if length(cond_ons{sess,cond}) > 0
                        SPM.Sess(sess).U(cond).name      = {cond_names{cond}};
                        SPM.Sess(sess).U(cond).ons       = cond_ons{sess,cond};
                        SPM.Sess(sess).U(cond).dur       = cond_dur{sess,cond};
                        SPM.Sess(sess).U(cond).P(1).name = 'none';%if name then define
                    end
                end
            end
            
            %==========================================================================
            % User-specified regressors
            
            
            % Realignment Parameters
            if realign_flag,
                % names of the user-specified regressors
                rnames = {'trans_X','trans_Y','trans_Z','rot_X','rot_Y','rot_Z'};
                %loop through sessions and assign Realignment Parameters
                for sess = 1:n_sess
                    thedir = sprintf('%s%s%1.0f',subj(sn).dir_ana_parent, '/NII/', subj(sn).func_dir{task}(sess) );
                    % 							thedir = sprintf('%s%1.0f',subj(sn).dir_data, subj(sn).func_dir{task}(sess) );
                    my_file = spm_select('List', thedir, '^rp.*.txt')
                    
                    SPM.Sess(sess).C.C =[];
                    SPM.Sess(sess).C.name = {};
                    
                    %load rp's from external file: VariableName: MovePara ()
                    tmp = load ([thedir ,filesep, my_file]);
                    tmp = tmp(1:n_scans, :); % remove trailing entries from realignment txt
                    
                    SPM.Sess(sess).C.C = [SPM.Sess(sess).C.C, tmp] ;
                    SPM.Sess(sess).C.name = { SPM.Sess(sess).C.name{:}, rnames{:} };
                end
            end;
            
            %==========================================================================
            
            
            % global normalization: OPTINS:'Scaling'|'None'
            %---------------------------------------------------------------------------
            if globalnorm_flag,
                SPM.xGX.iGXcalc    = 'Scaling';
            else
                SPM.xGX.iGXcalc    = 'None';
            end
            
            % low frequency confound: high-pass cutoff (secs) [Inf = no filtering]
            %---------------------------------------------------------------------------
            for sess=1:n_sess
                SPM.xX.K(sess).HParam = high_pass_sec; % high pass cut off: 2 x max distance between cond
            end
            
            % intrinsic autocorrelations: OPTIONS: 'none'|'AR(1) + w'
            %-----------------------------------------------------------------------
            SPM.xVi.form       = 'AR(1)';
            
            
            % specify data: matrix of filenames
            %===========================================================================
            all_files = [];% for concatenation of data over number of sessions
            
            for sess=1:n_sess
                thedir = sprintf('%s%s%1.0f',subj(sn).dir_ana_parent, '/NII/', subj(sn).func_dir{task}(sess) );
                % my_files = spm_select('List', thedir, '^swuaf.*.img');
                my_files = spm_select('List', thedir, ['^' file_img_prefix '.*.nii']);
                % add the full path to the functional imaging data
                my_dir_files = [repmat([thedir filesep],size(my_files,1),1) my_files];
                all_files = strvcat(all_files, my_dir_files);
                
                % global regress
                if globalregress_flag,
                    rnames = {'global_mean'};
                    for i=1:size(my_dir_files,1)
                        a = spm_read_vols(spm_vol(my_dir_files(i,:)));
                        bold1d = a(:);
                        if i==1,
                            bold2d = zeros([size(bold1d,1) size(my_dir_files,1)]);
                        end
                        bold2d(:,i) = bold1d;
                    end
                    globalmean = mean(bold2d,1);
                    
                    SPM.Sess(sess).C.C = [SPM.Sess(sess).C.C, globalmean'];
                    SPM.Sess(sess).C.name = { SPM.Sess(sess).C.name{:}, rnames{:} };
                end
            end
            
            SPM.xY.P = all_files;
            
            
            % Go to the correct directory before configurating the design matrix
            %======================================================================
            cd([subj(sn).dir_ana_parent,'/Stat/',dir_analysis])
            
            % Configure design matrix
            %===========================================================================
            %SPM = spm_fmri_spm_ui(SPM);
            SPM = spm_fmri_spm_ui(SPM);
            
            % Estimate parameters
            %===========================================================================
            SPM = spm_spm(SPM);
        end; % 	if ~do_contrasts_only
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                              CONTRASTS                              %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if do_contrasts_only,
            load(fullfile(subj(sn).dir_ana_parent, '/Stat/',dir_analysis,'SPM.mat'));
        end;
        
        % Debug: If you re-do contrasts, first remove xCon field:
        % SPM = rmfield(SPM,'xCon');
        
        %==========================================================================
        % SINGLE CONDITIONS POOLED OVER SESSIONS
        %==========================================================================
        if realign_flag,
            rps = zeros(1,6); % realign coord x,y,z,r,t,p
        end
        
        regc = [1 -1];
        if globalregress_flag,
            regrm = [rps 0];
        else
            regrm = [rps];
        end
        c                   = [repmat([regc regrm], 1,n_sess), zeros(1,n_sess)];
        cname               = 'attention vs resting';
        SPM.xCon            = spm_FcUtil('Set',cname,'T','c',c',SPM.xX.xKXs);
        
        c                   = -[repmat([regc regrm], 1,n_sess), zeros(1,n_sess)];
        cname               = 'resting vs attention';
        SPM.xCon(end + 1)   = spm_FcUtil('Set',cname,'T','c',c',SPM.xX.xKXs);
        
        % F tests
        tmp = eye( cond_n );          % F contrast
        tmp2 = repmat(regrm, cond_n,1);
        tmp = repmat( [tmp,tmp2] ,1,n_sess); % all sessions
        tmp3 = zeros( cond_n,n_sess); % session zeros
        tmp = [tmp tmp3];
        c = tmp;
        cname              = 'F all';
        SPM.xCon(end + 1)   = spm_FcUtil('Set',cname,'F','c',c',SPM.xX.xKXs);
        
        
%         cd /Users/soyoung/Atten_Corr_BOLD/00_BatchFiles/
        
        
        %==========================================================================
        %==========================================================================
        %==========================================================================
        %==========================================================================
        
        
        % ---------------------------------------------------------------------
        % and evaluate
        % ---------------------------------------------------------------------
        
        spm_contrasts(SPM);
        
    end
    %      spm_contrasts(SPM, 1:length(SPM.xCon(thec)));
end; % for thesubject

%%spm_contrasts(SPM,15); kann mit Index des Kontrastes aufgerufen werden
