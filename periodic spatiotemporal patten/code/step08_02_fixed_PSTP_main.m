
%% 
clear; clc;
PL=10; %%%% the window width of the Fixed PSTP

p2p='../All_Preproc_Data/Params_PSTP_HC.mat';    %%%Parameter setting
load(p2p,'nscn','nX','nt','nT','tres','p2B','nP','cth','ncth1','nitr','ITPstp','nsbj','ibY','iG2Y');   
IS=1:nsbj; 

a ='../Results_fixed/SubjPSTP/'; 
if ~exist(a,'dir'), mkdir(a); end

p2S=cell(nsbj,nP); for is=1:nsbj, for ip=1:nP % pth2 save SbjPSTPs
        p2S{is,ip}=[a 'S' num2str(is,'%03d') '_' num2str(ip)]; end; end

%%%%Read masks
load('../mask/DMN_TPN_parcels_combine.mat','mask_group_zheng','mask_group_fu');
pcc_corr_indi_zheng=mask_group_zheng;  %%%%DMN
pcc_corr_indi_fu=mask_group_fu;     %%%%TPN
pathpng='../Results_fixed/png_SbjPSTP/'; %%save the fit curves of the DMN and TPN 
if ~exist(pathpng,'dir'), mkdir(pathpng); end
pathpngz='../Results_fixed/png_SbjPSTP_zscore/';
if ~exist(pathpngz,'dir'), mkdir(pathpngz); end
%%
fprintf('Loading data\n'); load(p2B,'B');   %%%Time series for each subject
ratio_all = zeros(nsbj,1);   %%%Save the proportion of each person's window
ratio_r_positve_PSTP_allsub= zeros(nsbj,1);     %%The proportion of PSTP that is positively correlated with the optimal PSTP
ratio_p_sig_PSTP_allsub= zeros(nsbj,1);  %%The proportion of PSTP that is significantly positively correlated with the optimal PSTP
ratio_r_positve_PSTP_all_allsub= zeros(nsbj,1);   
ratio_p_sig_PSTP_all_allsub= zeros(nsbj,1);   

for is=IS

    D1=zeros(nX,nT,'single');
    for iscn=1:nscn
        D1(:,(iscn-1)*nt+(1:nt))=B{is,iscn};
    end
    D=D1;
    TT=cell(1,1); CT=zeros(1,nT,'single');
  
    for ip=1:nP    %%
    
        tic;
        
        if ip~=1; D=load(p2S{is,ip-1},'Dr'); D=D.Dr;end
        
        isip=sprintf('Sbj%d-PSTP%d-',is,ip);
        
        [PLh,PLc,PLe,esg,ITP,ssg,tsh] = parameter_detectingPSTP(PL,nt,nP,nscn);

        [PSTP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX,period,xp,yp_DMN,fun_DMN,yp_TPN,threshold_parceil_np1,threshold_parceil_np1_transition,...
                    fun_TPN,c1,c2,yc_DMN,yc_TPN,label_plot,x_L,x_R,x_point,y_point,...
                    ratio,ratio_r_positive_PSTP,ratio_p_sig_PSTP,ratio_r_positive_PSTP_all,ratio_p_sig_PSTP_all]=...
                    Detect_PSTP...
                    (pcc_corr_indi_zheng,pcc_corr_indi_fu,PLc,D,nscn,PL,cth{ip},ncth1(ip),nitr,ssg(ip),ITP{ip},PLh,tres,...
                    [isip 'f1detect'],ITPstp); % >> ITPl{is,ip} for limited segments
                
          if ip==1
            ratio_all(is) = ratio;
            ratio_r_positve_PSTP_allsub(is) = ratio_r_positive_PSTP;   
            ratio_p_sig_PSTP_allsub(is) = ratio_p_sig_PSTP;   
            ratio_r_positve_PSTP_all_allsub(is) = ratio_r_positive_PSTP_all;
            ratio_p_sig_PSTP_all_allsub(is) = ratio_p_sig_PSTP_all;  
          end
        plot_figure(xp,yp_DMN,yp_TPN,fun_DMN,fun_TPN,PL,...
            c1,c2, period,pathpng,pathpngz,ip,is,yc_DMN,yc_TPN,label_plot,x_point,x_L,x_R,y_point)   
       
        fprintf([isip 'f2regscn\n']);
        TT{ip}=PSTP; CT(ip,:)=C;
        [Dr,Cr,FCr]=PSTPf2regscn...
            (D1,TT(1:ip),CT(1:ip,:),nscn,PL,PLc,ibY,iG2Y,0);        
   
        
       save(p2S{is,ip},'PSTP','TMX','C','MET','ITER','TMPL','TMXTMPL',...
            'CTMPL','SCMX','Dr','Cr','FCr','threshold_parceil_np1','threshold_parceil_np1_transition','ratio','ratio_r_positive_PSTP',...
            'ratio_p_sig_PSTP','ratio_r_positive_PSTP_all','ratio_p_sig_PSTP_all');
       
       fprintf('Sbj%d-PSTP%d %dsec\n\n',is,ip,round(toc));
    end
end; clear D1 D
save('../Results_fixed/ratio.mat','ratio_all','ratio_r_positve_PSTP_allsub','ratio_p_sig_PSTP_allsub',...
    'ratio_r_positve_PSTP_all_allsub','ratio_p_sig_PSTP_all_allsub');
