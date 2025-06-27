%% This code is for presetting parameters in advance

%% Setting paths, parameters & indices
%% HC
clear; clc;
p2p='../All_Preproc_Data/Params_PSTP_HC.mat';    %%%Parameter setting

load(p2p,'nscn','nX','nt','nT','tres','p2B','nP','cth','ncth1','nitr','ITPstp','nsbj','ibY','iG2Y');
IS=1:nsbj;

a ='../Results_individual_tailored/SubjPSTP/';
if ~exist(a,'dir'), mkdir(a); end

p2S=cell(nsbj,nP); for is=1:nsbj, for ip=1:nP % pth2 save SbjPSTPs
        p2S{is,ip}=[a 'S' num2str(is,'%03d') '_' num2str(ip)]; end; end

%%%%Read masks
load('../mask/DMN_TPN_parcels_combine.mat','mask_group_zheng','mask_group_fu');
pcc_corr_indi_zheng=mask_group_zheng;  %%%%DMN
pcc_corr_indi_fu=mask_group_fu;     %%%%TPN
pathpng='../Results_individual_tailored/png_SbjPSTP/'; % dir2 save secondary outputs or temporary variables
if ~exist(pathpng,'dir'), mkdir(pathpng); end
pathpngz='../Results_individual_tailored/png_SbjPSTP_zscore/'; % dir2 save secondary outputs or temporary variables
if ~exist(pathpngz,'dir'), mkdir(pathpngz); end
%%
fprintf('Loading data\n');load(p2B,'B');   %%%Time series for each subject

numcir = 20;
period_np1= cell(nsbj,1);
ratio_all = zeros(nsbj,1);   %%%Save the proportion of each person's window
ratio_r_positve_PSTP_allsub= zeros(nsbj,1);     %%The proportion of PSTP that is positively correlated with the optimal PSTP
ratio_p_sig_PSTP_allsub= zeros(nsbj,1);  %%The proportion of PSTP that is significantly positively correlated with the optimal PSTP
ratio_r_positve_PSTP_all_allsub= zeros(nsbj,1);
ratio_p_sig_PSTP_all_allsub= zeros(nsbj,1);
% nt_b = nt; nT_b = nT;
for is=IS
    
    %     nt = nt_b; nT = nT_b;
    D1=zeros(nX,nT,'single');
    for iscn=1:nscn
        D1(:,(iscn-1)*nt+(1:nt))=B{is,iscn};
    end
    %     nT = 50; nt = 50;
    %     D1 = D1(:,1:nT);
    D=D1;
    TT=cell(nP,1); CT=zeros(nP,nT,'single');
    % For neuromodulation, it is sufficient to set nP = 1, which significantly speeds up the process.
    % Setting nP = 2 is equivalent to regressing out the PSTP and recalculating the PSTP using the residuals.
    % The residuals may disrupt the intrinsic structure of the DMN and TPN, making it difficult to identify a stable PSTP.
    % Therefore, PSTP2 was calculated only once without an iterative process.
    for ip=1:nP
        tic;
        if ip~=1
            if exist([p2S{is,ip-1} '.mat'], 'file')
                DT=load(p2S{is,ip-1});
                if isfield(DT, 'Dr')
                    D = DT.Dr;
                end
            end
        end
        isip=sprintf('Sbj%d-PSTP%d-',is,ip);
        if ip==1; PL0=10; PL=PL0; end
        %         PL0=10; PL=PL0;
        if ip==1
            period_m = zeros(1,numcir+1); period_m(1,1)=PL0;
            for cir=1:numcir
                
                [PLh,PLc,PLe,esg,ITP,ssg,tsh] = parameter_detectingPSTP(PL,nt,nP,nscn);
                
                [PSTP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX,period,xp,yp_DMN,fun_DMN,yp_TPN,threshold_parceil_np1,threshold_parceil_np1_transition,...
                    fun_TPN,c1,c2,yc_DMN,yc_TPN,label_plot,x_L,x_R,x_point,y_point,...
                    ratio,ratio_r_positive_PSTP,ratio_p_sig_PSTP,ratio_r_positive_PSTP_all,ratio_p_sig_PSTP_all]=...
                    Detect_PSTP...
                    (pcc_corr_indi_zheng,pcc_corr_indi_fu,PLc,D,nscn,PL,cth{ip},ncth1(ip),nitr,ssg(ip),ITP{ip},PLh,tres,...
                    [isip 'f1detect'],ITPstp); % >> ITPl{is,ip} for limited segments
                
                if ~isempty(PSTP)
                    if ~isempty(find( period_m == ceil(period)))
                        period_m(1,cir+1)=ceil(period); break;
                    else
                        PL=ceil(period);
                        
                        if PL<6       %%% the period length of the anti-correlation characteristic between DMN and TPN cannot be less than 6 TRs
                            PL=6; %If PSTP does not exist, it indicates that it is too short. For PL+1, recalculate PSTP
                            period_m(1,cir+1)=6+100; %
                        else
                            period_m(1,cir+1)=PL;
                        end
                    end
                else
                    PL=PL+1;
                    period_m(1,cir+1)=PL+100; %100 can be any value that is not a plausible period length, indicating that no optimal period could be identified for this subject.
                end
            end
            if ~isempty(PSTP)
                period_m(period_m==0)=[];
                period_np1{is,1} = period_m;
                ratio_all(is) = ratio;
                ratio_r_positve_PSTP_allsub(is) = ratio_r_positive_PSTP;
                ratio_p_sig_PSTP_allsub(is) = ratio_p_sig_PSTP;
                ratio_r_positve_PSTP_all_allsub(is) = ratio_r_positive_PSTP_all;
                ratio_p_sig_PSTP_all_allsub(is) = ratio_p_sig_PSTP_all;
            end
        else
            if  exist([p2S{is,ip-1} '.mat'], 'file') && isfield(DT, 'Dr')
                [PLh,PLc,PLe,esg,ITP,ssg,tsh] = parameter_detectingPSTP(PL,nt,nP,nscn);
                
                [PSTP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX,period,xp,yp_DMN,fun_DMN,yp_TPN,threshold_parceil_np1,threshold_parceil_np1_transition,...
                    fun_TPN,c1,c2,yc_DMN,yc_TPN,label_plot,x_L,x_R,x_point,y_point,...
                    ratio,ratio_r_positive_PSTP,ratio_p_sig_PSTP,ratio_r_positive_PSTP_all,ratio_p_sig_PSTP_all]=...
                    Detect_PSTP...
                    (pcc_corr_indi_zheng,pcc_corr_indi_fu,PLc,D,nscn,PL,cth{ip},ncth1(ip),nitr,ssg(ip),ITP{ip},PLh,tres,...
                    [isip 'f1detect'],ITPstp); % >> ITPl{is,ip} for limited segments
            end
        end
        
        if ~isempty(PSTP) && cir < numcir
            plot_figure(xp,yp_DMN,yp_TPN,fun_DMN,fun_TPN,PL,...
                c1,c2, period,pathpng,pathpngz,ip,is,yc_DMN,yc_TPN,label_plot,x_point,x_L,x_R,y_point)
            
            fprintf([isip 'f2regscn\n']);
            TT{ip}=PSTP; CT(ip,:)=C;
            [Dr,Cr,FCr]=PSTPf2regscn...
                (D1,TT(1:ip),CT(1:ip,:),nscn,PL,PLc,ibY,iG2Y,0);
            
            save(p2S{is,ip},'PSTP','TMX','C','MET','ITER','TMPL','TMXTMPL',...
                'CTMPL','SCMX','Dr','Cr','FCr','threshold_parceil_np1','threshold_parceil_np1_transition','ratio','ratio_r_positive_PSTP',...
                'ratio_p_sig_PSTP','ratio_r_positive_PSTP_all','ratio_p_sig_PSTP_all');
            if ip==1, save(p2S{is,ip},'period_m','-append'); end
            fprintf('Sbj%d-PSTP%d %dsec\n\n',is,ip,round(toc));
        end
    end
end; clear D1 D

save('../Results_individual_tailored/period_np1.mat','period_np1');
save('../Results_individual_tailored/ratio.mat','ratio_all','ratio_r_positve_PSTP_allsub','ratio_p_sig_PSTP_allsub',...
    'ratio_r_positve_PSTP_all_allsub','ratio_p_sig_PSTP_all_allsub')