
%At the group level
%the top 10% of positive connections positively correlated with the PCC are defined as the DMN
%the negative connections with the largest absolute values (most strongly negatively correlated with the PCC) are defined as the TPN.
%% 
nmD='HC'; % dataset name

% addpath './cifti-matlab-master/'
p2u='./Utils/'; addpath(genpath(p2u)); % pth2 cifti utilities
p2prep='../All_Preproc_Data/';  addpath(genpath(p2prep));
load([p2prep 'Prep_Params_' nmD '.mat'],'iXG','nt','nsbj','nscn','ivx','ixG','p2bgwcr','nTg','nVX','nvx');   
ilpcc5=[210 213 215 341]; %LPCC
I=[];
for i=1:length(ilpcc5)
    I=cat(1,I,ixG{ilpcc5(i)});  
end
iVxlpcc5=sort(I);  

%% For each subject, compute the correlation between the PCC and all brain parcels 
%%(triangles), select the top 10% of the strongest positive and negative connections respectively, 
%%and then take the union across subjects.
nG = length(ixG);
meanTime=zeros(nG,nt);
meanR = zeros(1,nG);
for is=1:nsbj 
    for i=1:nscn  
        load(p2bgwcr{is,i},'bpr');
        
        TimePCC = bpr(iVxlpcc5,:);
        meanTimePCC=mean(TimePCC);
        %%%%%%%%%%%%%%%%%%%%%The mean time series of each subject's parcels.%%%%%%%%%%%%%%%%%%%
        for j=1:nG
            meanTime(j,:)=mean(bpr(ixG{j},:));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%Compute the correlation between the PCC and the other parcels.%%%%%%%%%%%%%%%%%%%
        [r,~] = corr(meanTimePCC',meanTime'); 
        clear bpr
        meanR = meanR+r;
    end  
end
meanR = meanR/(nsbj*nscn);
thresh_upper = prctile(meanR,90);
thresh_lower = prctile(meanR,10);
r_upper = meanR .* (meanR > thresh_upper);
r_lower = meanR .* (meanR < thresh_lower);
rr = r_upper + r_lower;
rr(rr>0)=1;
rr(rr<0)=-1;

%%
O=nan(nVX,1); 
for nn=1:nG
 O(iXG{nn},1)=rr(1,nn);
end
%%%%%%%%%%%%%%%%%save
e=ft_read_cifti('empty.dtseries.nii'); 
e.dtseries=O;
ft_write_cifti('../mask/DMN_TPN_parcels_L',e,'parameter','dtseries'); 

%% save mat
mask_group=rr;
mask_group_zheng=rr;
mask_group_fu=rr;

mask_group_zheng(mask_group_zheng~=1)=0;   %%DMN

mask_group_fu(mask_group_fu~=-1)=0;    %%TPN
mask_group_fu(mask_group_fu~=0)=1;     %%TPN

save('../mask/DMN_TPN_parcels_L.mat','mask_group','mask_group_zheng','mask_group_fu');




