function [PSTP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX,period,xp,yp_DMN,fun_DMN,yp_TPN,threshold_parceil_np1,threshold_parceil_np1_transition,...
    fun_TPN,c1,c2,yc_DMN,yc_TPN,label_plot,x_L,x_R,x_point,y_point,...
    ratio,ratio_r_positive_PSTP,ratio_p_sig_PSTP,ratio_r_positive_PSTP_all,ratio_p_sig_PSTP_all]=...
    Detect_PSTP...
    (pcc_corr_indi_zheng,pcc_corr_indi_fu,PLc,D,nscn,PL,cth,ncth1,nitr,ssg,ITP,PLh,tres,s,stp)
%%%
[nX,nT]=size(D);
nt=nT/nscn;  
esg=nt-PL+1;  
nXL=nX*PL; %flattening
PLe=PL+sum(PLh); 

%% Flattening all segments & zscoring
SGf=zeros(nXL,nT,'single'); SGfn=zeros(nXL,nT,'single');
SGf_ex=zeros(nXL,nT,'single'); SGfn_ex=zeros(nXL,nT,'single');
for iscn=1:nscn
    for isg=ssg:esg   
        t=(isg:isg+PL-1)+(iscn-1)*nt;
        S=D(:,t);
        
        %%To load the DMN mask, see if the DMN is in the first moment
        %put the largest moment point of the DMN into the first moment point, and the rest of the delay
        S_DMN=(pcc_corr_indi_zheng*S)/length(find(pcc_corr_indi_zheng==1));
        [~,max_idx] = max(S_DMN);
        S_ex = [S(:,max_idx:end), S(:,1:max_idx-1)];
        S_ex = S_ex(:);
        SGf_ex(:,(iscn-1)*nt+isg)=S_ex; 
        S_ex=S_ex-sum(S_ex)/nXL;   
        S_ex=S_ex/sqrt(S_ex'*S_ex);  
        SGfn_ex(:,(iscn-1)*nt+isg)=S_ex;
        %%%%%SGfn, as a sliding window, is related to template computation
        S=S(:);
        SGf(:,(iscn-1)*nt+isg)=S; 
        S=S-sum(S)/nXL;  
        S=S/sqrt(S'*S);  
        SGfn(:,(iscn-1)*nt+isg)=S;
    end
end
clear S

%% Running the PSTP algorithm for initial segments of ITP
nITP=length(ITP); TMPL=cell(nITP,1); TMXTMPL=cell(nITP,1);
CTMPL=zeros(nITP,nT,'single'); SCMX=zeros(nITP,1,'single');
ITR=zeros(nITP,1,'single'); e=0.9999;
threshold_parceil_np1_transition =zeros(length(ITP),nitr+1);  %Save the threshold for each iteration
%%%%%%%%%%%%%% tmxsub_num=zeros(nITP+1,nitr);
%%%%%%%%%%%%%% itr_all=zeros(1,nITP);
for itp=1:nITP
    c=SGfn_ex(:,ITP(itp))'*SGfn;  
  
    cc=c;
    tmx=[];
    for numsub=1:nscn
        subcc=cc((numsub-1)*nt+1:numsub*nt);
        location=nt-PL+2:nt;
        subcc(location)=[];
        bb=prctile(subcc,cth(1));
        threshold_parceil_np1_transition(itp,1) = bb;
        [~,tmxsub]=findpeaks(subcc,'MinPeakHeight',bb,'MinPeakDistance',PL);  
        tmxsub=(numsub-1)*nt+tmxsub;
        tmx=[tmx tmxsub];
    end
    for iscn=1:nscn
        tmx( tmx==ssg+(iscn-1)*nt | tmx==esg+(iscn-1)*nt)=[];
    end
    
    nmx=length(tmx);
    cn=c-sum(c)/nT; cn=cn'/sqrt(cn*cn');  
    
    cn1=cn; cn2=cn; cn3=cn; itr=1;
    while itr<=nitr   
        if nmx<=1, break; end
        
        T=sum(SGf(:,tmx),2)/nmx;  
        T=T-sum(T)/nXL;
        T=T/sqrt(T'*T);  %Re-evaluate zscore, where T is equivalent to the updated template
        T = reshape(T, [nX, PL]);
        
        %%%%Update the template to ensure that the DMN-TPN maximum is at the first point in time
        T_DMN=(pcc_corr_indi_zheng*T)/length(find(pcc_corr_indi_zheng==1));
        [~,max_idx] = max(T_DMN);
        T_ex = [T(:,max_idx:end), T(:,1:max_idx-1)];
        T = T_ex;
        T = T(:);
        c=T'*SGfn;

        if itr<=ncth1-1, ith=1; else, ith=2; end
        cc=c;
        tmx=[];
        for numsub=1:nscn
            subcc=cc((numsub-1)*nt+1:numsub*nt);
            location=nt-PL+2:nt;
            subcc(location)=[];
            bb=prctile(subcc,cth(ith));
            threshold_parceil_np1_transition(itp,itr+1)=bb;
            [~,tmxsub2]=findpeaks(subcc,'MinPeakHeight',bb,'MinPeakDistance',PL);  
            tmxsub2=(numsub-1)*250+tmxsub2;
            
            tmx=[tmx tmxsub2];
            
        end
        
        for iscn=1:nscn
            tmx( tmx==ssg+(iscn-1)*nt | tmx==esg+(iscn-1)*nt)=[]; end
        nmx=length(tmx);
        
        cn=c-sum(c)/nT; cn=cn/sqrt(cn*cn');  
        if cn*cn1>e || cn*cn2>e || cn*cn3>e, break; end  
        cn3=cn2; cn2=cn1; cn1=cn'; itr=itr+1;
    end; clear cn cn1 cn2 cn3
    
    if nmx>1
        T=zeros(nX,PLe,'single');
        tS=tmx-PLh(1); 
        tE=tmx+PL-1+PLh(2); 
        for i=1:nmx, ts=tS(i); te=tE(i);
            zs=[];
            if ts<=0
                zs=zeros(nX,abs(ts)+1,'single');
                ts=1;
            end
            ze=[];
            if te>nT
                ze=zeros(nX,te-nT,'single');
                te=nT;
            end
            T=T+[zs D(:,ts:te) ze];  
        end
        TMPL{itp}=T/nmx; clear T   %Any starting point, the final template after each iteration.
        TMXTMPL{itp}=single(tmx);  %The time point in which the local maximum value is stored
        CTMPL(itp,:)=single(c);    %c here is the correlation timecourse of the final template
        SCMX(itp)=single(sum(c(tmx)));  %The sum of the local maximum correlation values of the final template's correlation timecourse
        ITR(itp)=itr;  %Number of iterations
        
    end
    if ~rem(itp,stp), fprintf([s '-ITP%d\n'],itp); end
end
if rem(itp,stp), fprintf([s '-ITP%d\n'],itp); end
clear SGf SGfn c
%% Delete the Windows greater than nitr+1
% calculate the proportion of Windows that are not exceeded, and then perform the following calculation
lo_sup_nitr = find(ITR>nitr);  
%% 
%%%%All the RSTPP satisfying the conditions are calculated by flatten correlation
% and then the maximum value of the sum of the rows is found
%%%%%The RSTP corresponding to the maximum value is used as the optimal RRP

stacod = zeros(1,nITP);  %%Used to store RSTPs that meet the conditions

DMN_f_all = cell(nITP,1);
TPN_f_all = cell(nITP,1);
pc_DMN_all = cell(nITP,1);
pc_TPN_all = cell(nITP,1);
yp_DMN_all = zeros(nITP,PL);
yp_TPN_all = zeros(nITP,PL);
yc_DMN_all = zeros(nITP,PL);
yc_TPN_all = zeros(nITP,PL);
c1_all = zeros(nITP,1);    %How many times the DMN in each RSTP is fitted
c2_all = zeros(nITP,1);    %How many times the TPN in each RSTP is fitted

for kk=1:nITP    
    PSTP_detect=TMPL{kk}; 
    PSTP_detect_PL=PSTP_detect(:,PLc);
    PSTP_DMN_series=(pcc_corr_indi_zheng*PSTP_detect_PL)/length(find(pcc_corr_indi_zheng==1));
    PSTP_TPN_series=pcc_corr_indi_fu*PSTP_detect_PL/length(find(pcc_corr_indi_fu==1));
    %% conditions
    %%%%%%%%%%Determine whether RSTP meets the first condition%%%%%%%%%%%%%%
    xp=1:PL;
    yp_DMN=PSTP_DMN_series;
    yp_TPN=PSTP_TPN_series;
    %%%%%%%%%%%Fitting function
    warning('off')
    [xp,yc_DMN,pc_DMN,c1] = plot_polyfit(xp,yp_DMN,2,PL-1);  
    [xp,yc_TPN,pc_TPN,c2] = plot_polyfit(xp,yp_TPN,2,PL-1);  
    warning('on')
    [max_val,max_idx] = max(PSTP_DMN_series);
    
    if max_idx==1 && max_val>0  && (yp_DMN(1)>yp_TPN(1))
       
        yp_DMN_all(kk,:) = yp_DMN;
        yp_TPN_all(kk,:) = yp_TPN;
        yc_DMN_all(kk,:) = yc_DMN;
        yc_TPN_all(kk,:) = yc_TPN;
        pc_DMN_all{kk} = pc_DMN;
        pc_TPN_all{kk} = pc_TPN;
        
        c1_all(kk,1) = c1;
        c2_all(kk,1) = c2;
    
        %%
        %%%Find the extreme point of the function (two networks can only have one extreme point on the left and one extreme point on the right) 
        %%%and the maximum and minimum values of the function
        
        %%%%%%%%%%%Create function%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        syms DMN_f(x)
        DMN_f(x)=pc_DMN;
        syms TPN_f(x)
        TPN_f(x)=pc_TPN;
        
        %%%%赋值
        DMN_f_all{kk} = DMN_f(x);
        TPN_f_all{kk} = TPN_f(x);
        stacod(1,kk)=kk;    %%%Save the RSTPs that meets the conditions
    end
end
clear kk
%% %% Correlation between all windows
parwintim_allwin = zeros(nXL,nITP);
for jj=1:nITP   
    parwin=zeros(nX,PL);
    for ii=1:PL
        parwin(:,ii) = TMPL{jj}(:,PLc(ii));    %360*PL
    end
    parwin = parwin(:);  %%
    parwin = parwin-sum(parwin)/nXL;   
    parwin = parwin/sqrt(parwin'*parwin);   
    parwintim_allwin(:,jj) = parwin;   
end
clear jj ii
[r_all,p_all]=corr(parwintim_allwin);  

%% %%%%%%%%%%%%%%%%%%%%the correlation between RSTPs that meet the conditions%%%%%%%%%%%%%%%%%%%%
stacod(lo_sup_nitr)=0; ratio = (nITP-length(find(stacod==0)))/nITP;
stacod(stacod==0)=[];
num_stacod=length(stacod);parwintim = zeros(nXL,num_stacod);

for jj=1:num_stacod    
    parwin=zeros(nX,PL);
    for ii=1:PL
        parwin(:,ii) = TMPL{stacod(jj)}(:,PLc(ii));   
        
    end
    
    parwin = parwin(:);  
    parwin = parwin-sum(parwin)/nXL;   
    parwin = parwin/sqrt(parwin'*parwin);   
    parwintim(:,jj) = parwin;    
    
end
clear jj ii
%%
%%%%%%%%%%（3600*num_stacod)%%%%%%%%%%%%
[r,p]=corr(parwintim); 
rr_mean=mean(r,2);    %%%%%%Sum the correlation values by row
[~,isscmx]=sort(rr_mean,'descend');
Label = 0;    %%%If the number of intersection points is less than 1, choose the second best RSTP, and so on
%% Preset parameters
PSTP=[];period=[];label_plot=[];x_L=[];
TMX=[];C=[];ITER=[];MET=[];x_R=[];threshold_parceil_np1=[];
for ss=1:num_stacod
    close all;
    ITP1 = stacod(isscmx(ss));
    
    [c1,c2,pc_DMN,pc_TPN,DMN_f,yp_DMN,TPN_f,yp_TPN,yc_DMN,yc_TPN] = ...
        iter_optimal_PSTP(ITP1,DMN_f_all,...
        yp_DMN_all,pc_DMN_all,pc_TPN_all,TPN_f_all,yp_TPN_all,...
        yc_DMN_all,yc_TPN_all,c1_all,c2_all);

    fun_DMN=matlabFunction(DMN_f);
    fun_TPN=matlabFunction(TPN_f);
    
    %%%%%%%%%%%%%%%%Finding Intersection Points%%%%%%%%%%%%%%%%%%%%%%%%%
    [x_point,y_point] = fsolve_xy(pc_DMN,pc_TPN,PL);
    
    %%%%%%%%%%%%%%%Find the extremum point%%%%%%%%%%%%
    [x_DMN,y_DMN] = inflection_point(DMN_f,PL,x);   
    [x_TPN,y_TPN] = inflection_point(TPN_f,PL,x);
%%%%Distinguish the situation where there is no intersection point between 1~PL, where the period is the distance between the two (optimal) peaks
    num_med = length(find(x_point>1 & x_point<PL));
    %%%%Calculate the mean and standard deviation of the signal
    A(1:PL)=yp_TPN;
    A(PL+1:PL+PL)=yp_DMN;
    B=[mean(A) std(A)];
    D=B(1)+3*B(2);
    E=B(1)-3*B(2);
    %%
    if y_point<D & y_point>E
        if length(y_DMN)>1 && length(y_TPN)>1
            if num_med == 0   
                if length(y_DMN)>2 && length(y_TPN)>2  
                    Label=1; label_plot = 1;
                    x1 = x_DMN(1); x2= x_DMN(3);
                    x3 = x_TPN(1); x4= x_TPN(3);
                    
                    y1 = fun_DMN(x1)-fun_TPN(x1);
                    y2 = fun_DMN(x2)-fun_TPN(x2);
                    y3 = fun_DMN(x3)-fun_TPN(x3);
                    y4 = fun_DMN(x4)-fun_TPN(x4);
                    %%%%%Find the point with the farthest distance between the two extreme points of DMN and TPN
                    if y1>y3;x_L = x1;else;x_L=x3; end
                    if y2>y4;x_R = x2;else;x_R=x4;end                    
                else
                    x5 = x_DMN(2); x6= x_TPN(2);
                    y5 = fun_DMN(x5)-fun_TPN(x5);
                    y6 = fun_DMN(x6)-fun_TPN(x6);
                    %%%%%Find the point with the closest distance between the two extreme points of DMN and TPN
                    if y5<y6;x_R = x5;else;x_R=x6;end
              %%The middle intersection point is empty and there are no more than two extreme points, but there are intersection points outside of 1~PL
                    if ~isempty(x_point)  
              %%If there is no intersection point in the middle, it is impossible to have three additional intersection points outside of 1~PL, at most two
                        if length(x_point)==1    
                            x_L = x_point;
                        else
                            x_1=abs(x_point(1)-1);
                            x_2=abs(x_point(2)-PL);
                            if x_1<x_2
                                x_L = x_point(1);
                            else
                                x_L = x_point(2);
                            end
                        end
            %%The middle intersection point is empty and there are no more than two extreme points, and there are no intersection points outside of 1~PL             
                    else 
                      
                        x_1 = fit_extend_R(PL,fun_DMN,fun_TPN,D);
                        x_2 = fit_extend_L(fun_DMN,fun_TPN,D);
                        if  (x_1-PL)>(x_2-1)
                            x_L=x_2;
                        else
                            x_L=x_1;
                        end
                    end
                    
                end
                period = abs(double(x_R - x_L));  
      
            elseif num_med>3   
                Label=1; label_plot = 2;
                x_point(x_point<1)=[];x_point(x_point>PL)=[];
                period = double(x_point(4)-x_point(2));
                x_L=x_point(2);x_R=x_point(4);     
            elseif num_med==3
                Label=1; label_plot=3;
                x_point(x_point<1)=[];x_point(x_point>PL)=[];
                period = double(x_point(3)-x_point(1));
         %%%%There are two intersection points between two fitting curves, but the fitting curve has two or three intersection points
            elseif num_med==2  
        %%If there are two intersection points in the middle, determine the distance from the first intersection point to 1 
        %%and the distance from the second intersection point to the last time point
        %%%and retain the larger one
               
                Label=1;  label_plot=4;
                x_pointf = x_point;
                x_pointf(x_pointf<1)=[];x_pointf(x_pointf>PL)=[];
                x_diff(1,1) = abs( x_pointf(1)-1);     
                x_diff(1,2) =  abs(x_pointf(end)-PL);
                
                [~,local_x] = max(x_diff);
                if local_x(1)==1 
                    x_L = fit_extend_L(fun_DMN,fun_TPN,D);
                    x_R = x_pointf(2);
                else          
                    x_R = fit_extend_R(PL,fun_DMN,fun_TPN,D);
                    x_L = x_pointf(1);
                end
                period = double(x_R-x_L);   
            elseif num_med==1                
                Label=1;
                if length(x_point)==3
                    label_plot=5;  
                    period = double(x_point(end)-x_point(1));
                    x_L=x_point(1);x_R=x_point(end); 
                elseif length(x_point)==1 || length(x_point)==2
                    label_plot=6;
                    x_pointf=x_point;
                    loca_med_point = find((x_pointf>1)&(x_pointf<PL)==1); 
                    x_pointf(loca_med_point)=[];
                    if ~isempty(x_pointf)
                        
                        if x_pointf<1   
                            x_L = x_pointf(1);
                            x_R = fit_extend_R(PL,fun_DMN,fun_TPN,D);
                        else
                            x_R = x_pointf(1);
                            x_L = fit_extend_L(fun_DMN,fun_TPN,D);
                        end
                    else 
                        
                        x_R = fit_extend_R(PL,fun_DMN,fun_TPN,D);
                        x_L = fit_extend_L(fun_DMN,fun_TPN,D);
                    end
                    if ~isempty(x_R) && ~isempty(x_L)
                        period = double(x_R-x_L);
                    end
                    
                end
            end
        end
    end
    if isempty(x_R) || isempty(x_L)
        Label=0;
    end
    if Label==1
        %%%%%%%%%%%%%%%Selecting PSTP%%%%%%%%%%%%%%%%%%%
        [PSTP,TMX,C,ITER,MET,~,threshold_parceil_np1] = optimal_detected_PSTP(TMPL,ITP1,TMXTMPL,CTMPL,ITR,tres,threshold_parceil_np1_transition);
        %%%
        lo_stacod_ITP1 = find(stacod==ITP1);
        r_positive_PSTP = r(lo_stacod_ITP1,:);
        ratio_r_positive_PSTP = length(find(r_positive_PSTP>0))/num_stacod;
        p_sig_PSTP = p(lo_stacod_ITP1,:);
        p_sig_PSTP(find(r_positive_PSTP<0)) = 1;
        ratio_p_sig_PSTP = length(find(p_sig_PSTP<0.05))/num_stacod;
        %%%% 
        r_positive_PSTP_all = r_all(ITP1,:);
        ratio_r_positive_PSTP_all = length(find(r_positive_PSTP_all>0))/nITP;
        p_sig_PSTP_all = p_all(ITP1,:);
        p_sig_PSTP_all(find(r_positive_PSTP_all<0)) = 1;
        ratio_p_sig_PSTP_all = length(find(p_sig_PSTP_all<0.05))/nITP;
        break;
    else
        ratio_r_positive_PSTP=[];
        ratio_p_sig_PSTP=[];
        ratio_r_positive_PSTP_all=[];
        ratio_p_sig_PSTP_all=[];
    end

end


