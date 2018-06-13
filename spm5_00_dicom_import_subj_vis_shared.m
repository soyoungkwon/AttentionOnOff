% batch for spm5
% Batch script for dicom import of single or multiple subjects
% essential batch files: spm_dicom_headers, spm_dicom_convert
% 2010/03/23 Soyoung Kwon
% clear all;
% subfoldList;
% task = 4;
clear all; clc;
spm_read_shared_parameters
% main_dir = '/Users/soyoung/Atten_Corr_BOLD/';
main_dir = strrep(pwd, '00_BatchFile/New_spm', '');

% task_dir = {'/10_Vis_Short/', '/11_Motion_OnOff/', '/12_Vis_Sess/'};
% task_dir = {'/10_Vis_Short/', '/11_Motion_OnOff/', '/Archive/Prev_sm6mm/12_Vis_Sess/'};
task_dir = {'/04_AttenMotionShort/', '/05_MotionOnOff/', '/07_AttenMotionLong/'};
ima_pardir = [main_dir, '03_IMA/'];
dir_final = [main_dir, '00_BatchFile/New_spm/'];
% dir_final = [main_dir, '/00_BatchFiles/SPM/'];
dummy = 0;

for sub = 2:size(subj_dir,2)    % each subject
    for task = 2%1:3%1:size(task_dir,2)%1:3 % each task % edit  
            
        %%-----Each task = which folder? ----%%
        nii_pardir = [main_dir, task_dir{task}];
        mkdir([nii_pardir, subj_dir{sub}], 'NII') % make NIFTI folder
        mkdir([nii_pardir, subj_dir{sub}], 'Stat') % make NIFTI folder
        all_dir{task}{sub} = [subj(sub).func_dir{task}, subj(sub).struct_dir];
        nii_subjdir=[nii_pardir, subj_dir{sub},'/NII/'];
        nfolder = size(all_dir{task}{sub},2);
       
        %%-----Each Task-----%%
        for i=1:nfolder
            % dir: read ima & save nii
            ima_dir = sprintf('%s%s%s%1.0f%s', ima_pardir, subj_dir{sub}, '/',all_dir{task}{sub}(i),'/');
            
            % make nifti folder
            if i~=size(all_dir{task}{sub},2) % functional folder
                mkdir(nii_subjdir, sprintf('%1.0f',all_dir{task}{sub}(i)));
                nii_dir=sprintf('%s%1.0f',nii_subjdir, all_dir{task}{sub}(i));
            elseif i == nfolder  % structural folder
                mkdir(nii_subjdir, 'struct_mean');
                nii_dir=[nii_subjdir, 'struct_mean'];
            end
            
            % save nii
            ima_hdr=spm_select('List',ima_dir,['.*ima']); %header
            full_hdr= [repmat(ima_dir, [size(ima_hdr,1) 1]), ima_hdr]; % header with full directory
            hdr=spm_dicom_headers(full_hdr);
            cd(nii_dir);
            spm_dicom_convert(hdr, 'all','flat','nii');
            
            % delte dummy
            for dum=1:dummy
                delete(sprintf('%s%1.0f%s', '*00',dum,'-01.nii'));
            end
        end
    end
    cd(dir_final);
end