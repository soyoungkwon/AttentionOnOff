% ROI connectivity

% 1. BOLD preprocess (HF cutoff, global mean regress)
% 2. Mean BOLD (sliding, whole)
% 3. BOLD connectivity (sliding, whole)
% time window can vary
% preprocess can vary
% --> save these parameter somewhere

clear all;
% change read_analysis_info_roi
analysis.preproc = 1;
analysis.mergebold = 1;
analysis.boldsubj = 0;
analysis.mean_plot = 0;
analysis.psych = 0;
analysis.regMean = 0;
analysis.psychbold = 0;
analysis.connectivity = 0;

read_group_info_roi
n_subj = size(subj,2);

%=============== Preprocessing ===========%
if analysis.preproc
    for s=3:n_subj%:11
        for sess=1:subj{s}.nsess
            run_BOLDproc(s,sess)
            run_BOLDClean4D(s,sess)
            run_extractBOLD_ROI(s,sess) % here
            run_BOLDTcs_resample(s,sess)
%             run_regressMean_hrf(s,sess) % maybe not used
            run_BOLDTcs_epoch_sort(s,sess)
        end
    end
end
%----- merge bold------%
if analysis.mergebold
    for s=3:n_subj
        run_MergeBOLD(s) % boldresample(sess_length, n_vox) --> boldatten(cond_length, n_vox, n_sess*2)
%         run_regressMean_vox(s) % each voxel regress % probably here!!!
%         run_Regress_ButtonPress_Vox(s) % used for paper!!!!!!
%         run_Regress_StiButton_Vox(s)
%         run_regressMeanTri_vox(s)
    end
end
%----- regress mean (ROI mean base) ------% % Not use now!
if analysis.regMean
    for s=1
        run_Regress_Mean_Subj(s)
    end
%    run_RegressMean_Subj    % boldatten, boldrest(roi,tcs,tri). data structure remains same
%    run_RegressMean_Group   % boldatten, boldrest(roi,tcs,tri)
end
if analysis.boldsubj   % merge ROI--> then each subject has 1 file
%     run_BOLDTcs_Subj       % probably run_BOLDTcs_Subj_Vox replaces it
    run_BOLDTcs_Subj_Vox   % boldatten(tcs,vox,tri) --> boldatten(roi,tcs,tri)
%     run_BOLDTcs_epoch_sort_subj

end
%=============== Main analysis ===========%
%----------- mean tcs --------REST_reg_sub46_franzi_PCC_fwe0001------%
if analysis.mean_plot % load BOLD_sphere_craf_reg1_none_low0_P_N11_TR1....
    run_ROI_tcs_plot
end


if analysis.connectivity
    %-------- Individual Connectivity-----%
    run_connectivity_subj              % corrAtten, corrRest(roi,roi,tri)
    run_connectivity_subj_sliding      % corrAtten, corrRest(roi,roi,tri,n_window)
    run_connectivity_regMean           % corrAtten, corrRest(roi,roi,tri) - regMean!!
    
    %------ corr btw Psych & Connectivity ----%
    run_Corr_Psych_Conn_Subj    % corrPsychConn(roi,roi)
    run_Corr_Psych_Conn_Trial   % corrPsychConn(roi,roi,tri)
    run_Corr_PsychHRF_Conn_Trial % corrPsychConn(roi,roi,tri) ??
end