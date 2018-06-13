% conn-hierarchy
% after group_stat_LRboth
plot_subj=1;
% n_subj = size(gROIAttenBoth,3);%Both,3);%Both,3);
condname={'vision','atten','diff'};
visA = [1:5]; n_vis=5;
meanReg=1;%mean(visA); % to make V1 as 0!!!! new
for LR=0:1%LR = 1; % 1-left 0-right
    for c=1:3,% 1-Rest, 2-Atten, 3-Atten-Rest
        if LR, subI=[1:20 41:60]; else, subI=[21:60]; end
        if c==1, gROIConn=gROIRestBoth(:,:,subI);
        elseif c==2, gROIConn=gROIAttenBoth(:,:,subI);
        elseif c==3, gROIConn=gROIAttenBoth(:,:,subI)-gROIRestBoth(:,:,subI);
        end
        n_subj = size(gROIConn,3);
        % gROIConn = gROIRestBoth;%-
        
        % ==== DAN - VIS  & DMN - VIS===== %
        i_fef=size(gROIConn,1)-3; 
        if LR,
            FEFVIS = squeeze(gROIConn(visA,i_fef,:)); % left
            IPSVIS = squeeze(gROIConn(visA,i_fef+1,:));
            VISLP = squeeze(gROIConn(visA,i_fef+2,:));
            VISCC = squeeze(gROIConn(visA,i_fef+3,:));
        else,
            FEFVIS = squeeze(gROIConn(i_fef,visA,:)); % right
            IPSVIS = squeeze(gROIConn(i_fef+1,visA,:));
            VISLP = squeeze(gROIConn(i_fef+2,visA,:));
            VISCC = squeeze(gROIConn(i_fef+3,visA,:));
        end
        
        fef{c}=ezfit([1:n_vis]-meanReg, mean(FEFVIS(1:n_vis,:),2), 'a*x+b', [2 0]);
        ips{c}=ezfit([1:n_vis]-meanReg, mean(IPSVIS(1:n_vis,:),2), 'a*x+b', [2 0]);
        lp{c}=ezfit([1:n_vis]-meanReg, mean(VISLP(1:n_vis,:),2), 'a*x+b', [2 0]);
        cc{c}=ezfit([1:n_vis]-meanReg, mean(VISCC(1:n_vis,:),2), 'a*x+b', [2 0]);
        
        % slope
        fefs{c} = nanmean(FEFVIS,2);
        ipss{c} = nanmean(IPSVIS,2);
        lps{c} = nanmean(VISLP,2);
        ccs{c} = nanmean(VISCC,2);
        if LR==1, LRname = 'L'; else, LRname = 'R'; end
        figure('name',[condname{c} LRname]);
        subplot(2,2,1);
        plotCI(FEFVIS); title(sprintf('FEF-VIS y=%0.3fX + %0.3f', fef{c}.m(1), fef{c}.m(2))); ylim([0 0.3]); if c==3, ylim([0 0.3]); end
        subplot(2,2,2);
        % figure('name',condname{c});
        plotCI(IPSVIS); title(sprintf('IPS-VIS y=%0.3fX + %0.3f', ips{c}.m(1), ips{c}.m(2)));  ylim([0 0.3]); if c==3, ylim([0 0.3]); end
        subplot(2,2,3);
        % figure('name',condname{c});
        plotCI(VISLP); title(sprintf('LP-VIS y=%0.3fX + %0.3f', lp{c}.m(1), lp{c}.m(2)));  ylim([-0.5 0.1]);
        subplot(2,2,4);
        % figure('name',condname{c});
        plotCI(VISCC); title(sprintf('CC-VIS y=%0.3fX + %0.3f', cc{c}.m(1), cc{c}.m(2)));  ylim([-0.5 0.1]);
        
        figure('name',condname{c});
        subplot(2,2,1); plot(FEFVIS, '--.'); title('FEF-VIS'); axis([0 6 -0.4 0.8]);
        subplot(2,2,2); plot(IPSVIS, '--.'); title('IPS-VIS'); axis([0 6 -0.4 0.8]);
        subplot(2,2,3); plot(VISLP, '--.'); title('LP-VIS'); axis([0 6 -0.4 0.8]);
        subplot(2,2,4); plot(VISCC, '--.'); title('CC-VIS'); axis([0 6 -0.4 0.8]);
        
        if plot_subj
            %----- subject- stat-----%
            for s=1:n_subj
                F=ezfit([1:n_vis], FEFVIS(:,s), 'a*x+b', [2 0]);
                fef_b1(c,s)=F.m(1);
                fef_b2(c,s)=F.m(2);
                I=ezfit([1:n_vis], IPSVIS(:,s), 'a*x+b', [2 0]);
                ips_b1(c,s)=I.m(1);
                ips_b2(c,s)=I.m(2);
                L=ezfit([1:n_vis], VISLP(:,s), 'a*x+b', [2 0]);
                lp_b1(c,s)=L.m(1);
                lp_b2(c,s)=L.m(2);
                C=ezfit([1:n_vis], VISCC(:,s), 'a*x+b', [2 0]);
                cc_b1(c,s)=C.m(1);
                cc_b2(c,s)=C.m(2);
            end
        end
    end % cond
    
    %------- difference -----------%
    [h_fef_b1,p_fef_b1]=ttest(fef_b1(2,:)-fef_b1(1,:));
    [h_fef_b2,p_fef_b2]=ttest(fef_b2(2,:)-fef_b2(1,:));
    [h_ips_b1,p_ips_b1]=ttest(ips_b1(2,:)-ips_b1(1,:));
    [h_ips_b2,p_ips_b2]=ttest(ips_b2(2,:)-ips_b2(1,:));
    [h_lp_b1,p_lp_b1]=ttest(lp_b1(2,:)-lp_b1(1,:));
    [h_lp_b2,p_lp_b2]=ttest(lp_b2(2,:)-lp_b2(1,:));
    [h_cc_b1,p_cc_b1]=ttest(cc_b1(2,:)-cc_b1(1,:));
    [h_cc_b2,p_cc_b2]=ttest(cc_b2(2,:)-cc_b2(1,:));
    
    if LR==1, % left
        rest_gain_dan(1,:)=fef_b1(1,:); rest_gain_dan(3,:)=ips_b1(1,:);
        atten_gain_dan(1,:)=fef_b1(2,:); atten_gain_dan(3,:)=ips_b1(2,:);
        rest_gain_dmn(1,:)=lp_b1(1,:); rest_gain_dmn(3,:)=cc_b1(1,:);
        atten_gain_dmn(1,:)=lp_b1(2,:); atten_gain_dmn(3,:)=cc_b1(2,:);
        rest_base_dan(1,:)=fef_b2(1,:); rest_base_dan(3,:)=ips_b2(1,:);
        atten_base_dan(1,:)=fef_b2(2,:); atten_base_dan(3,:)=ips_b2(2,:);
        rest_base_dmn(1,:)=lp_b2(1,:); rest_base_dmn(3,:)=cc_b2(1,:);
        atten_base_dmn(1,:)=lp_b2(2,:); atten_base_dmn(3,:)=cc_b2(2,:);
    elseif LR==0, % right
        rest_gain_dan(2,:)=fef_b1(1,:); rest_gain_dan(4,:)=ips_b1(1,:);
        atten_gain_dan(2,:)=fef_b1(2,:); atten_gain_dan(4,:)=ips_b1(2,:);
        rest_gain_dmn(2,:)=lp_b1(1,:); rest_gain_dmn(4,:)=cc_b1(1,:);
        atten_gain_dmn(2,:)=lp_b1(2,:); atten_gain_dmn(4,:)=cc_b1(2,:);
        rest_base_dan(2,:)=fef_b2(1,:); rest_base_dan(4,:)=ips_b2(1,:);
        atten_base_dan(2,:)=fef_b2(2,:); atten_base_dan(4,:)=ips_b2(2,:);
        rest_base_dmn(2,:)=lp_b2(1,:); rest_base_dmn(4,:)=cc_b2(1,:);
        atten_base_dmn(2,:)=lp_b2(2,:); atten_base_dmn(4,:)=cc_b2(2,:);
    end
end
figure('name', 'DAN');
subplot(2,1,1); plotCI(atten_gain_dan-rest_gain_dan)
subplot(2,1,2); plotCI(atten_base_dan-rest_base_dan)

figure('name', 'DMN');
subplot(2,1,1); plotCI(atten_gain_dmn-rest_gain_dmn)
subplot(2,1,2); plotCI(atten_base_dmn-rest_base_dmn)

rest_gain=[fef_b1(1,:); ips_b1(1,:); lp_b1(1,:); cc_b1(1,:)];
atten_gain=[fef_b1(2,:); ips_b1(2,:); lp_b1(2,:); cc_b1(2,:)];

rest_base=[fef_b2(1,:); ips_b2(1,:); lp_b2(1,:); cc_b2(1,:)];
atten_base=[fef_b2(2,:); ips_b2(2,:); lp_b2(2,:); cc_b2(2,:)];

figure; subplot(2,2,1); p_atten_gain=plotCI(atten_gain); title('atten gain');
subplot(2,2,2); p_rest_gain=plotCI(rest_gain); title('rest gain');
subplot(2,2,3); p_diff_gain=plotCI(atten_gain-rest_gain); title('diff gain');
subplot(2,2,4); p_diff_base=plotCI(atten_base-rest_base); title('diff base');

figure('name','model');
subplot(2,2,1); title('FEF');
hold on; plot([1:n_vis+1], fef{1}.m(1)*[1:n_vis+1]+fef{1}.m(2), '-.');
hold on; plot([1:n_vis+1], fef{2}.m(1)*[1:n_vis+1]+fef{2}.m(2), '-.');
subplot(2,2,2); title('IPS');
hold on; plot([1:n_vis+1], ips{1}.m(1)*[1:n_vis+1]+ips{1}.m(2), '-.');
hold on; plot([1:n_vis+1], ips{2}.m(1)*[1:n_vis+1]+ips{2}.m(2), '-.');
subplot(2,2,3); title('LP');
hold on; plot([1:n_vis+1], lp{1}.m(1)*[1:n_vis+1]+lp{1}.m(2), '-.');
hold on; plot([1:n_vis+1], lp{2}.m(1)*[1:n_vis+1]+lp{2}.m(2), '-.');
subplot(2,2,4); title('CC');
hold on; plot([1:n_vis+1], cc{1}.m(1)*[1:n_vis+1]+cc{1}.m(2), '-.');
hold on; plot([1:n_vis+1], cc{2}.m(1)*[1:n_vis+1]+cc{2}.m(2), '-.');
