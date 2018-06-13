% batch for SPM5
% Batch script for preprocessing of single or multiple subjects
% realignment, unwarpping, coregistration,
% segmentation, normalisation, and smoothing
% Normalisation done after segmentation of the structural image
% Projet Dyslexia (Cathy Price)....
% Mohamed Seghier, FIL, 08.06.2006
% Extended new version, 08.03.2007
% Extended abartels 2008.10
% session :2

clear all
spm_defaults
global defaults
defaults.modality = 'FMRI' ;
spm_read_shared_parameters;


% ====== DO WHAT: ========
timeslice_flag		= 1; % Do timeslicing. Reference slice = first acquired slice.
tslice_sequ		= 'interleaved'; % 'ascending'|'descending'|'interleaved' % interleaved: does odd/even thing, see below.

realign_unwarp_flag = 1; % 1: realign only; 0: unwarp and realign.
realign_reslice_flag = 1; % 0: estimate only.
%     does not work when unwarping.
%     MAKE SURE THERE IS A MEAN IMG IF YOU PLAN TO DO coreg_struct2mean_flag!
% 1: estimate and reslice also. (makes mean img).
% realign_unwarp_flag = 2; % 2: unwarp and realign

coreg_struct2mean_flag = 1;
normalization = 0;
% normalize_oldfashioned = 0; % 1: normalize using source= meanFunc, target= templates/EPI.nii
%    same for structural (source = struct, target=templates/T1.nii)
% 0: new way: using struct segmentation.
% 		normalize_struct =1; % only appicable if normalize_oldfashioned == 1.
% 0: do not normalize struct.
% 1: normalize struct to struct_template.
%		- TIMESLICE OR NOT
%		- .nii or .img
fprefix_orig = 'f^*'; %
fsuffix = '.*nii';
fsuffix_struct = '.*nii'; % nii
% fsuffix = '.*img';

%%%%%%%%%%%%%    Main INFO %%%%%%%%%%%%%%%%%%%%%%
% main_dir = '/Users/soyoung/Atten_Corr_BOLD/';
main_dir = strrep(pwd, '00_BatchFile/New_spm', '');
% task_dir = {'10_Vis_Short', '11_Motion_OnOff', '12_Vis_Sess'};

% smooth_fwhm = [10 10 10];

%%%%%%%%%%%%%%%  END SUBJECT SPECIFIC INFO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsub = length(subj) ;
for sub=2:nsub
    for task = 2%1:3,
        fprefix = fprefix_orig;
        
        Nsessions = size(subj(sub).func_dir{task},2);%Nses{sub}; % number of sessions (by default = 4)
        folSubject{sub} = [main_dir, task_dir{task}, subj_dir{sub}, '/NII/'];%subj(sub).dir_data;
        folStruct{sub}  = [folSubject{sub}, 'struct_mean/'];
        
        clear P dirs folData ;
        disp('###############################################################')
        disp(['preprocessing of subject number ', num2str(sub), ' .........'])
        disp('###############################################################')
        [files,dirs]=spm_select('List',folSubject{sub},'') ; % to get functional directory list
        for i=1:Nsessions
            folData{i} = fullfile(folSubject{sub}, deblank(dirs(i+2,:)));%+2,:))) ;
            disp(folData{i}) ;
        end
        nses = length(folData) ;
        disp('') ;
        
        struct_image = spm_select('List',folStruct{sub},['^ms^*' fsuffix_struct]) ;
        if isempty(struct_image), struct_image = spm_select('List',folStruct{sub},['^s^*' fsuffix_struct]) ; end
        
        if 1,
            if timeslice_flag,
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %    % Timeslice
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                TR = subj(sub).TR; % specific to subject
                nslices = subj(sub).nslices; % specific to subject
                jobs{1}.temporal{1}.st.tr = TR; % sec
                jobs{1}.temporal{1}.st.ta = TR - TR/nslices;
                jobs{1}.temporal{1}.st.nslices = nslices;
                
                switch lower(tslice_sequ)
                    case 'ascending',
                        jobs{1}.temporal{1}.st.so = 1:nslices;
                        jobs{1}.temporal{1}.st.refslice = 1;
                    case 'descending',
                        jobs{1}.temporal{1}.st.so = nslices:-1:1;
                        jobs{1}.temporal{1}.st.refslice = 1; % not sure - does this refer to acquired slice?
                    case 'interleaved',
                        % REFERENCE SLICE: first in time:
                        if mod(nslices,2) == 0; % EVEN n slices:
                            jobs{1}.temporal{1}.st.so = [ [2:2:nslices],[1:2:nslices-1] ];
                            jobs{1}.temporal{1}.st.refslice = 2; % just to be sure, could be 1 or 2
                        else % ODD n slices:
                            jobs{1}.temporal{1}.st.so = [ [1:2:nslices],[2:2:nslices-1] ];
                            jobs{1}.temporal{1}.st.refslice = 1;
                        end;
                    otherwise
                        disp(' -------- GIVE CORRECT TSLICE_SEQU STRING!! -------- ');
                        break; return;
                end;
                
                for ses=1:nses
                    tmp1 = [] ;
                    tmp2 = [] ;
                    tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                    if isempty(tmp1{1}) % try lower case names for files
                        tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                    end
                    for ttt=1:size(tmp1{1},1)
                        tmp2{1}(ttt,:) = fullfile(folData{ses}, tmp1{1}(ttt,:));
                    end
                    jobs{1}.temporal{1}.st.scans(ses)  = {tmp2{1}} ;
                end
                fprefix = ['a' fprefix];
                
                spm_jobman('run' , jobs) ;
                clear jobs ;
                
            end; % timeslice_flag
            
            if realign_unwarp_flag == 2,
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %    % Realignment (Coregister) + Unwarping (Reslice)
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                jobs{1}.spatial{1}.realignunwarp.eoptions.quality = defaults.realign.estimate.quality ;
                jobs{1}.spatial{1}.realignunwarp.eoptions.fwhm = 5 ;
                jobs{1}.spatial{1}.realignunwarp.eoptions.rtm = 1 ;
                jobs{1}.spatial{1}.realignunwarp.eoptions.einterp = defaults.realign.estimate.interp ;
                jobs{1}.spatial{1}.realignunwarp.uwroptions.uwwhich = [2 1] ;
                jobs{1}.spatial{1}.realignunwarp.uwroptions.rinterp = defaults.realign.write.interp ;
                jobs{1}.spatial{1}.realignunwarp.uwroptions.wrap = [0 0 0] ;
                jobs{1}.spatial{1}.realignunwarp.uwroptions.mask = defaults.realign.write.mask ;
                % load images
                for ses=1:nses
                    tmp1 = [] ;
                    tmp2 = [] ;
                    %       tmp1{1} = spm_select('List',folData{ses},'^af^*.*img');
                    tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                    if isempty(tmp1{1}) % try lower case names for files
                        %           tmp1{1} = spm_select('List',folData{ses},'^af^*.*img');
                        tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                    end
                    for ttt=1:size(tmp1{1},1)
                        tmp2{1}(ttt,:) = fullfile(folData{ses}, tmp1{1}(ttt,:));
                    end
                    jobs{1}.spatial{1}.realignunwarp.data(ses).scans  = {tmp2{1}} ;
                    jobs{1}.spatial{1}.realignunwarp.data(ses).pmscan = {''} ;
                end
                disp(sprintf('##### realigning (coregestring)........ and unwarping (and reslicing)........'))
                fprefix = ['u' fprefix];
                spm_jobman('run' , jobs) ;
                clear jobs ;
                
            elseif realign_unwarp_flag==1,
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %    % Realignment NO reslice (Coregister)
                %    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if realign_reslice_flag,
                    %	ESTIMATE & WRITE
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.quality = 0.9;         % default
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.sep = 4;               % default
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.fwhm = 5;              % default
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.rtm = 1;
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.interp = 2;            % default
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.wrap = [0 0 0];        % default
                    jobs{1}.spatial{1}.realign{1}.estwrite.eoptions.weight = {};           % don't weight
                    
                    jobs{1}.spatial{1}.realign{1}.estwrite.roptions.which  = [2 1];        % create mean image only when reslicing
                    jobs{1}.spatial{1}.realign{1}.estwrite.roptions.interp = 4;            % default
                    jobs{1}.spatial{1}.realign{1}.estwrite.roptions.wrap   = [0 0 0];      % no wrap (default)
                    jobs{1}.spatial{1}.realign{1}.estwrite.roptions.mask   = 1;            % enable masking (default)
                else
                    %	ESTIMATE ONLY
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.quality = 0.9;         % default
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.sep = 4;               % default
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.fwhm = 5;              % default
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.rtm = 1;
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.interp = 2;            % default
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.wrap = [0 0 0];        % default
                    jobs{1}.spatial{1}.realign{1}.estimate.eoptions.weight = {};           % don't weight
                end;
                % load images
                for ses=1:nses
                    tmp1 = [] ;
                    tmp2 = [] ;
                    tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                    if isempty(tmp1{1}) % try lower case names for files
                        tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                    end
                    for ttt=1:size(tmp1{1},1)
                        tmp2{1}(ttt,:) = fullfile(folData{ses}, tmp1{1}(ttt,:));
                    end
                    if realign_reslice_flag,
                        jobs{1}.spatial{1}.realign{1}.estwrite.data(ses)  = {tmp2{1}} ;
                    else
                        jobs{1}.spatial{1}.realign{1}.estimate.data(ses)  = {tmp2{1}} ;
                    end
                end
                
                if realign_reslice_flag,%
                    disp(sprintf('##### realigning and reslicing........'))
                    fprefix = ['r' fprefix];
                else
                    disp(sprintf('##### realigning NOT reslicing........'))
                    % fprefix = ['r' fprefix];
                end;
                spm_jobman('run' , jobs) ;
                clear jobs ;
                
            end; % if realign_unwarp_flag == 1
            
            
            if coreg_struct2mean_flag,
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % Coregister structural 2 mean
                % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                jobs{1}.spatial{1}.coreg{1}.estimate.eoptions.cost_fun = defaults.coreg.estimate.cost_fun ;
                jobs{1}.spatial{1}.coreg{1}.estimate.eoptions.sep = defaults.coreg.estimate.sep ;
                jobs{1}.spatial{1}.coreg{1}.estimate.eoptions.tol = defaults.coreg.estimate.tol ;
                jobs{1}.spatial{1}.coreg{1}.estimate.eoptions.fwhm = defaults.coreg.estimate.fwhm ;
                
                % ---------- select mean func image: --------------
                % NOTE: in case you did only realign:estimate, then no mean image was written.
                %		In that case you need to either take first func or make your own mean image.
                %		Problem with making your own mean image: you would have to make it
                %			with non-resliced functional data, ie you end up averaging over subject-movement.
                %		Problem with using first functional image: it's noisier than a mean.
                %		Solution: make mean using the first 10 images: gets rid of noise, and only limited motion.
                
                tmp = spm_select('List',folData{1},['^mean^*' fsuffix]);
                % --- if NO MEAN IMAGE found: ---
                if isempty(tmp),
                    % -- make mean of first ten images:
                    ip  = spm_select('FPList',folData{1},['^' fprefix fsuffix]); % list all images
                    tmp = min( size(ip,1) ,10); % 10 or less imgs
                    ipv = spm_vol(ip( 1:tmp , :)); % map first 10 or less imgs
                    mpv = rmfield(ipv(1),{'pinfo','private'}); % vol info for new mean img
                    mpv.fname = fullfile(folData{1},['mean_first_10_' fprefix(1) '_.' fsuffix(end-2:end)]); % new fname for mean img
                    
                    spm_imcalc(ipv, mpv, 'mean(X)',{1 0 0}); % write new mean img (into first session dir)
                    tmp = spm_select('List',folData{1},['^mean^*' fsuffix]); % name of mean img
                    % tmp = spm_select('List',folData{1},['^' fprefix fsuffix]); % take first functional or make new mean img
                end;
                tmp = tmp(1,:); % take first in case there are more
                
                jobs{1}.spatial{1}.coreg{1}.estimate.ref = {fullfile(folData{1}, tmp )} ;
                struct_image = spm_select('List',folStruct{sub},['^ms^*' fsuffix_struct]) ;
                if isempty(struct_image), struct_image = spm_select('List',folStruct{sub},['^s^*' fsuffix_struct]) ; end
                jobs{1}.spatial{1}.coreg{1}.estimate.source = {fullfile(folStruct{sub}, struct_image)} ;
                jobs{1}.spatial{1}.coreg{1}.estimate.other = {''} ;
                disp(sprintf('##### Coregistering structural volume to mean image.........'))
                spm_jobman('run' , jobs) ;
                clear jobs ;
            end
        end; % if 0 / 1
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Segmentation of the structural
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jobs{1}.spatial{1}.preproc.output.GM = [0 0 0] ;
        jobs{1}.spatial{1}.preproc.output.WM = [0 0 0] ;
        jobs{1}.spatial{1}.preproc.output.CSF = [0 0 0] ;
        jobs{1}.spatial{1}.preproc.output.biascor = 1 ;
        jobs{1}.spatial{1}.preproc.output.cleanup = defaults.segment.write.cleanup ;
        jobs{1}.spatial{1}.preproc.opts.tpm = {defaults.segment.estimate.priors} ;
        jobs{1}.spatial{1}.preproc.opts.ngaus = [3 2 2 5];
        jobs{1}.spatial{1}.preproc.opts.regtype = defaults.segment.estimate.affreg.regtype ;
        jobs{1}.spatial{1}.preproc.opts.warpreg = 1 ;
        jobs{1}.spatial{1}.preproc.opts.warpco = defaults.normalise.estimate.cutoff ;
        jobs{1}.spatial{1}.preproc.opts.biasreg = defaults.segment.estimate.reg ;
        jobs{1}.spatial{1}.preproc.opts.biasfwhm =  75 ;
        jobs{1}.spatial{1}.preproc.opts.samp = defaults.segment.estimate.samp ;
        jobs{1}.spatial{1}.preproc.opts.msk = {''} ;
        jobs{1}.spatial{1}.preproc.data = {fullfile(folStruct{sub}, struct_image)} ;
        disp(sprintf('##### Segmentation of the coregistered structural volume............'))
        spm_jobman('run' , jobs) ;
        clear jobs ;
        
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Normalisation of FUNCTIONAL data
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if normalization == 1,
            jobs{1}.spatial{1}.normalise{1}.write.roptions.preserve = 0 ;
            jobs{1}.spatial{1}.normalise{1}.write.roptions.bb = defaults.normalise.write.bb ;
            jobs{1}.spatial{1}.normalise{1}.write.roptions.vox = defaults.normalise.write.vox ;
            jobs{1}.spatial{1}.normalise{1}.write.roptions.interp = defaults.normalise.write.interp ;
            jobs{1}.spatial{1}.normalise{1}.write.roptions.wrap = defaults.normalise.write.wrap ;
            tmp = struct_image(1); % 1st char of struct_image
            jobs{1}.spatial{1}.normalise{1}.write.subj(1).matname = {fullfile(folStruct{sub}, spm_select('List',folStruct{sub},['^' tmp '.*_seg_sn*.*mat']))} ;
            P = cell(1,nses);
            for ses=1:nses
                tmp1 = [] ;
                tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
                for ttt=1:size(tmp1{1},1)
                    P{ses}(ttt,:) = fullfile(folData{ses}, tmp1{1}(ttt,:));
                end
            end
            % include mean func image:
            tmp = spm_select('List',folData{1},['^mean^*' fsuffix]);
            tmp = tmp(1,:); % only first in case there are more..
            P{nses+1} = fullfile(folData{1}, tmp) ;
            jobs{1}.spatial{1}.normalise{1}.write.subj(1).resample = {P{:}} ;
            disp(sprintf('##### writing normalised functional images................. '))
            fprefix = ['w' fprefix];
            spm_jobman('run' , jobs) ;
            clear jobs ;
        end
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Normalisation of STRUCTURAL data
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jobs{1}.spatial{1}.normalise{1}.write.roptions.preserve = 0 ;
        jobs{1}.spatial{1}.normalise{1}.write.roptions.bb = defaults.normalise.write.bb ;
        jobs{1}.spatial{1}.normalise{1}.write.roptions.vox = [1 1 1] ;
        jobs{1}.spatial{1}.normalise{1}.write.roptions.interp = defaults.normalise.write.interp ;
        jobs{1}.spatial{1}.normalise{1}.write.roptions.wrap = defaults.normalise.write.wrap ;
        tmp = struct_image(1); % 1st char of struct_image
        jobs{1}.spatial{1}.normalise{1}.write.subj(1).matname = {fullfile(folStruct{sub}, spm_select('List',folStruct{sub},['^' tmp '.*_seg_sn*.*mat']))} ;
        P = cell(1,1);
        P{1} = fullfile(folStruct{sub}, ['m' struct_image]) ;   % ms*.img
        jobs{1}.spatial{1}.normalise{1}.write.subj(1).resample = {P{:}} ;
        disp(sprintf('##### writing normalised STRUCT ................. '))
        spm_jobman('run' , jobs) ;
        clear jobs ;
        
        %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Smoothing
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        jobs{1}.spatial{1}.smooth.fwhm = [smooth_fwhm] ;
        jobs{1}.spatial{1}.smooth.dtype = 0 ;
        P = cell(1,nses);
        for ses=1:nses
            tmp1 = [] ;
            tmp1{1} = spm_select('List',folData{ses},['^' fprefix fsuffix]);
            for ttt=1:size(tmp1{1},1)
                P{ses}(ttt,:) = fullfile(folData{ses}, tmp1{1}(ttt,:));
            end
        end
        jobs{1}.spatial{1}.smooth.data = strvcat(P{:}) ;
        disp(sprintf('##### smoothing of the normalised functional images................. '))
        spm_jobman('run' , jobs) ;
        clear jobs ;
    end
    
end
