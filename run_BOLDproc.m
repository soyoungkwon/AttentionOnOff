function run_BOLDproc(s,sess)

read_group_info_roi
read_analysis_info_roi
read_shared_info

h=figure
%---- plot raw bold ------%
bold4d = reshape(bold2d, [n_tsteps, X,Y,Z]);
meanImgs = squeeze(mean(bold4d,1));
% meanImgs = reshape(meanImg2D, [X,Y,Z]);
h1=subplot(6,2,1); imagesc(meanImgs(:,:,30));
h2=subplot(6,2,2); plot(realign_params);

%---- load bold from roi-----%
Ts = 1;
if plot_roi,
    [bold2d_roi] = loadBoldRoi(bold2d, roi_indx, Ts);
    h3=subplot(6,2,3); plot(bold2d_roi); title('one ROI');
    h4=subplot(6,2,4); plot(bold2d_roi(:,sampleV));
    %----- load bold BRIAN -----%
    [bold2d_roi_brain] = loadBoldRoi(bold2d, roi_indx_brain, Ts);
    h5=subplot(6,2,5); plot(bold2d_roi_brain); title('only brain');
    h6=subplot(6,2,6); plot(bold2d_roi_brain(:,sampleV));
end
%------- processing ------%
fprintf('\n');
fprintf('\nProcessing subject %s', sub);
fprintf('\nSess %d', sess);
fprintf('\nYour TR is %3.2f% secs.', TR);
fprintf('\n');

% ===== RAW -->(%) ======= %
if settings.P
    [bold2d] = boldPchange(bold2d, brainIndx);
%     bold2d = bold2d./repmat(mean(bold2d,1), [size(bold2d,1) 1]);
end

% ===== Regress out ====== %
if ~strcmpi(settings.regress, 'None')
    if strcmpi(settings.regress, 'Motion')
        disp('Will regress out movement parameters and keep residuals.');
        [bold2d] = regressTcs(bold2d, realign_params, brainIndx);
    elseif strcmpi(settings.regress, 'MotionGlobal')
        disp('Will regress out movement parameters + global signal and keep residuals');
        brainTcs = bold2d(:, brainIndx);
        globalmean = mean(brainTcs,2);
        globalmean = globalmean - mean(globalmean);
        [bold2d] = regressTcs(bold2d, [realign_params globalmean], brainIndx);
    elseif strcmpi(settings.regress, 'MotionWM')
        disp('Will regress out movement parameters + WM signal and keep residuals');
        roiName = spm_select('List', [Dir.ROIspace, subj{s}.name], ['^rc2.*', '.*.nii']);
        
        %------ load roi index ---------%
        [wm_indx, wm_mask] = loadRoi(roiName, [Dir.ROIspace, subj{s}.name]);
        WMTcs = bold2d(:, wm_indx);
        WMmean = mean(WMTcs,2);
        WMmean = WMmean - mean(WMmean);
        [bold2d] = regressTcs(bold2d, [realign_params WMmean], brainIndx);
    end
end
if plot_roi,
    ts = ([settings.regress]==0);
    [bold2d_roi_brain] = loadBoldRoi(bold2d, roi_indx_brain, ts);
    h7=subplot(6,2,7); plot(bold2d_roi_brain); title('regress');
    h8=subplot(6,2,8); plot(bold2d_roi_brain(:,sampleV));
end
% ====== Filter temporally ======%
if settings.filter
    fprintf('\nWill temporally filter-cutoff: %d secs, high-filter %s \n', HFcutoff, filterOption);
    [bold2d] = filterTcs(bold2d, TR, HFcutoff, filterOption, brainIndx);
end
if plot_roi,
    [bold2d_roi_brain] = loadBoldRoi(bold2d, roi_indx_brain, ts);
    h9=subplot(6,2,9); plot(bold2d_roi_brain); title('filter');
    h10=subplot(6,2,10); plot(bold2d_roi_brain(:,sampleV));
end
% ----- save .nii ----- %%
bold4d = reshape(bold2d', [X,Y,Z,n_tsteps]);
for tt=1:n_tsteps % time
    bold3d = squeeze(bold4d(:,:,:,tt));
    if tt<10, time = sprintf('%d%d%d', 0, 0, tt); elseif tt<100, time = sprintf('%d%d', 0, tt); else time = sprintf('%d', tt); end
%     tmp = sprintf('craf_reg%d_hrf%d_sub%s_%s.nii', settings.regress, settings.filter, sub, time);
%     sprintf('craf_reg%d_hrf%d_sub%s');
    tmp = sprintf('%s_%s.nii', commonfix, time);    
    V.fname = [Dir.Cleaned, tmp];
    V = spm_write_vol(V, bold3d);
end
% ------- Spatial Smoothing
fprintf('\nTime course is getting smoothed with a gaussian kernel of %d mm\n', FWHM);
if FWHM > 0
    spatialSmooth(sub, Dir.Cleaned, FWHM)
    close
end
if save_figure, 
    figName = sprintf('BOLDplot_reg%s_HRF%d_sm%d_sub%s_sess%d', [settings.regress], [settings.filter], FWHM, sub, sess);
    saveas(h, [Dir.Fig,figName], 'fig');
end
close % close figure