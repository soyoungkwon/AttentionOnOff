function run_BOLDTcs_resample(s,ses)%,roi)
% resampling ---%
% load each subject roi bold

%-----load group & analysis info-----%
read_group_info_roi
read_analysis_info_roi
read_shared_short

%---- some param
resample = 1;

%----- load bold_roi
for roi=1:n_roi%27:32%n_roi-1:n_roi,%1:n_roi%
    
%     roiBoldTcsContain = sprintf('BOLD_sphere_%ssub%0.2d_sess%d_.*%s.*', commonfix, subj{s}.number,  ses, whichroi{roi});
    roiBoldTcsContain = sprintf('BOLD_%s%0.2d_%ssub%0.2d_sess%d_.*%s.*', roiShape, voxMin, commonfix, subj{s}.number, ses, whichroi{roi});

    roiBoldTcsNames = spm_select('List', Dir.ROIBOLD, roiBoldTcsContain);
    roiBoldTcsName = strrep(roiBoldTcsNames(1,:), ' ', '');
    load([Dir.ROIBOLD, roiBoldTcsName]); % bold2d_roi_subj[s,n_tcs_sess_tr, vox, sess]
    [n_TR, n_vox] = size(bold2d_roi);
    
    %------- resample time info -----%
    TR = subj{s}.TR;
    if resampleOn,
        oriT = 0:TR:TR*n_TR-1;%(sess_length+3);
        newT = 0:1/Fs:sess_length-1;%TR*n_TR-3;%(sess_length);
        postfix = [sprintf('_TR%d', 1/Fs)];
    end

    boldresample = zeros([sess_length, n_vox]);
    for v=1:n_vox,
        boldtcs = bold2d_roi(:,v);
        if isnan(mean(boldtcs)) || mean(boldtcs) == 0 % only resample when there is a value
            boldresample(:,v) = zeros([1, sess_length]);
        else
            boldresample(:,v) = interp1(oriT, boldtcs, newT);
        end
    end

    % save the resampled bold roi
    roiBoldTcsUpsName = strrep(roiBoldTcsName, '.mat', postfix);
    % roiBoldTcsUpsName = sprintf('BOLD_sphere_%ssub%0.2d_TR%d', commonfix, subj{s}.number, 1/Fs);
    save([Dir.ROITcsUps, roiBoldTcsUpsName], 'boldresample');
    
end % roi

end