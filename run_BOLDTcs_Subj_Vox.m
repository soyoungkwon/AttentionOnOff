% boldatten(tcs, vox, tri) --> boldatten(roi, tcs, tri)
% 19_ROI_Tcs_Subj_Vox_RegMean/ --> 20_ROI_Tcs_Subj_ROI/

% save the BOLD_ROI (each subject base)
% plot the bold tcs (each individual)

% clear all;
%-----load group & analysis info-----%
read_group_info_roi
read_analysis_info_roi
read_shared_short

regMean =0;
if regMean == 1,
    regfix = 'regMeanSubjVox_';%'regButtonMeanSubjVox_';%
    Dir.ROITcsSort = '../19_ROI_Tcs_Subj_Vox_RegMean/';%'../19_ROI_Tcs_Subj_Vox_RegButtonMean/';%'17_ROI_Tcs_Sort_RegMean/';
    endfix = '';
elseif regMean == 0,
    regfix = '';
    Dir.ROITcsSort = '../18_ROI_Tcs_Subj_Vox/';
    endfix = '';
elseif regMean == 2,
    regfix = 'regCorrRMeanSubjVox_';
    Dir.ROITcsSort = '../19_ROI_Tcs_Subj_Vox_RegButtonMean/';%'18_ROI_Tcs_Subj_Vox_RegMean/';
    endfix = '';%'_regMean';
elseif regMean == 4,
    regfix = 'TestregMeanSubjVox_';
    Dir.ROITcsSort = '../19_ROI_Tcs_Subj_Vox_RegMean/';
    endfix = '';
end
saveOn = 1;
plotGroup = 1;
saveGroup = 0;
n_tri = 4; %%%%% change!!!!!!!!
bmax = 0.1;
plot_gROI = 1;
plotsubtract = 0;

% n_tri = 12;

if plotGroup
    boldattenGroup = zeros([n_roi, cond_length, n_tri, n_subj]);
    boldrestGroup = zeros([n_roi, cond_length, n_tri, n_subj]);
end
% save the BOLD_ROI
for s=13:n_subj
    n_sess = subj{s}.nsess;
    if n_sess > 6, n_sess = 6; end % limit the maximum number of session
    n_tri = n_sess*2;
    
    % preassign the boldatten
    boldattensub = zeros([n_roi, cond_length, n_tri]);
    boldrestsub = zeros([n_roi, cond_length, n_tri]);    
    
    % each subject
    for roi=1:n_roi
        roiBoldTcsSortContain = sprintf('BOLD_%s%0.2d_%s%ssub%0.2d_sess.*%s.*TR%d.*%s', roiShape, voxMin, commonfix, regfix, subj{s}.number,whichroi{roi}, 1/Fs, endfix);
        roiBoldTcsSortList = spm_select('List', Dir.ROITcsSort, roiBoldTcsSortContain);
        roiBoldTcsSortName = strrep(roiBoldTcsSortList(1,:), ' ', '');
        load([Dir.ROITcsSort, roiBoldTcsSortName], 'boldatten', 'boldrest'); % boldatten
        
        boldattensub(roi,:,:) = squeeze(nanmean(boldatten(:,:,1:n_tri),2));
        boldrestsub(roi,:,:) = squeeze(nanmean(boldrest(:,:,1:n_tri),2));
    end
    
    % save the tcs (subj)
    boldatten = boldattensub;
    boldrest = boldrestsub;
    if plotGroup
        boldattenGroup(:,:,:,s) = boldattensub;
        boldrestGroup(:,:,:,s) = boldrestsub;
    end
    % plot the tcs
    figure('name', sprintf('sub%0.2d', subj{s}.number));
    for roi=1:n_roi
        boldattentmp = squeeze(boldatten(roi,:,:)); boldresttmp = squeeze(boldrest(roi,:,:));
        if nanmean(boldatten(roi,:),2)>90, boldattentmp = squeeze(boldatten(roi,:,:))-100; boldresttmp = squeeze(boldrest(roi,:,:))-100; end
        subplot(round(n_roi/4),4,roi); plot(nanmean(boldattentmp,2), 'r--');  xlim([0 cond_length]);hold on; ylim([-bmax bmax]);
        subplot(round(n_roi/4),4,roi); plot(nanmean(boldresttmp,2), 'k--');  xlim([0 cond_length]); hold on; ylim([-bmax bmax]);
        if plotsubtract
            subplot(round(n_roi/4),4,roi); plot(nanmean(boldattentmp(roi,:,:),3)-nanmean(boldresttmp(roi,:,:),3), 'b'); xlim([0 cond_length]); %ylim([-bmax bmax]);
        end
    end
    if saveOn
        roiBoldTcsName = sprintf('BOLD_%s%0.2d_%s%ssub%0.2d%s', roiShape, voxMin, commonfix, regfix, subj{s}.number, endfix);
%         roiBoldTcsName = sprintf('BOLD_%s%0.2d_%s%sroi%0.2d_sub%0.2d%s', roiShape, voxMin, commonfix, regfix, n_roi, subj{s}.number, endfix);
        
        save([Dir.ROITcsSubj, roiBoldTcsName], 'boldatten', 'boldrest');
    end
end

if plotGroup
    figure;
    for roi=1:n_roi
%         if nanmean(boldattenGroup(roi,:,:) > 90, boldattenGroup(roi,:,:,:) = boldattenGroup(roi,:,:,:)-100; boldrestGroup(roi,:,:) = boldrestGroup(roi,:,:); end
        if nanmean(boldattenGroup(roi,:),2) > 90, boldattenG = squeeze(boldattenGroup(roi,:,:,:))-100; boldrestG = squeeze(boldrestGroup(roi,:,:,:))-100; 
        else, boldattenG = squeeze(boldattenGroup(roi,:,:,:)); boldrestG = squeeze(boldrestGroup(roi,:,:,:)); end
        subplot(round(n_roi/4),4,roi); plot(nanmean(nanmean(boldattenG(:,:,:),2),3), 'r--'); xlim([0 cond_length]); if ~regMean, ylim([-bmax bmax]); end; hold on;
        subplot(round(n_roi/4),4,roi); plot(nanmean(nanmean(boldrestG(:,:,:),2),3), 'k--'); xlim([0 cond_length]); hold on;
        %    subplot(6,4,roi); plot(nanmean(boldattenGroup(roi,:,:),3)-nanmean(boldrestGroup(roi,:,:),3), 'b');
    end
    if saveGroup
        roiBoldGroupName = sprintf('BOLD_%s%0.2d_%s%sN%d', roiShape, voxMin, commonfix, regfix, n_subj);
        save([Dir.GroupSummary, roiBoldGroupName], 'boldattenGroup', 'boldrestGroup');
    end
    % boldattenGroup(roi, tcs, tri, s)
    if plot_gROI
        figure;
        for g1=1:nG
            subplot(2,round(nG/2),g1); 
            atri = atri(atri<=n_tri);%n_tri); %%%%%%%%%% ======= change ========== %%%%%%
            rtri = rtri(rtri<=n_tri);%n_tri);
            boldattengROI(g1,:) = nanmean(nanmean(nanmean(boldattenGroup(ROIgNrBoth{g1},:,atri,:),3),4),1);
            boldrestgROI(g1,:) = nanmean(nanmean(nanmean(boldrestGroup(ROIgNrBoth{g1},:,rtri,:),3),4),1);
            boldsubtractgROI(g1,:) = boldattengROI(g1,:)-boldrestgROI(g1,:);
            if mean(boldattengROI(g1,:)) > 90, boldattengROI(g1,:) = boldattengROI(g1,:)-100; end
            plot(boldattengROI(g1,:), 'r--');%nanmean(nanmean(nanmean(boldattenGroup(ROIgNr{g1},:,atri,:),1),3),4), 'r--');
            if mean(boldrestgROI(g1,:)) > 90, boldrestgROI(g1,:) = boldrestgROI(g1,:)-100; end
            hold on; subplot(2,round(nG/2),g1); plot(boldrestgROI(g1,:), 'k--');
            title(ROIgName{g1});            xlim([0 120]);          if ~regMean,  ylim([-bmax bmax]); end
        end        
    end
end