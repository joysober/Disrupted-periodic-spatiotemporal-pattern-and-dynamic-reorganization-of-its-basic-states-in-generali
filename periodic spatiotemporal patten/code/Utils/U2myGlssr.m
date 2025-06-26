clear; clc; addpath './cifti-matlab-master/'

%% Indices of Glasser's parcels
G=ft_read_cifti('Glasser.dlabel.nii');
G1=G.indexmax; g1=G1(~isnan(G1));

nG=max(G1); iXG=cell(nG,1); ixG=cell(nG,1);
for i=1:nG, iXG{i}=single(find(G1==i)); ixG{i}=single(find(g1==i)); end
save('myGlssr.mat','nG','iXG','ixG');

ilpcc5=[210 213 214 215 341]; %LPCC种子点，属于不同的parcels？？？这个可能需要看原文，Glasser是怎么划分的。这五个数是代表的五个parcels
ilv2=184; %可能是左侧早期视觉
ilm1=188; %可能是左侧M1（运动皮层）
ilpf=328; %可能是LPFC（左侧前额叶皮层）
I=[]; 
for i=1:5
    I=cat(1,I,ixG{ilpcc5(i)});  %将ilpcc5对应的各三角片的位置提取出来（去除了NaN后的位置）
end
ivxlpcc5=sort(I);  %对位置进行排序
ivxlv2=ixG{ilv2}; ivxlm1=ixG{ilm1}; ivxlpf=ixG{ilpf};
save('myGlssr.mat','ilpcc5','ilv2','ilm1','ilpf',...,
    'ivxlpcc5','ivxlv2','ivxlm1','ivxlpf','-append');
%
%% Yeo's 7 RSNs
%%% >> needs workbench installed in C:/workbench
% addpath(genpath('./cifti-gifti/'));
% Y=ciftiopen('RSN.dlabel.nii','C:/workbench/bin_windows64/wb_command');
% save('RSN.mat','Y');
load('RSN.mat','Y');  %四列值各代表什么？？？第一列最小37，最大44，后面会把37-44划分为7个网络
Y=Y.cdata;

load('myHCPcft.mat','ivx','nvxc','nVX'); ivxc=ivx(1:nvxc);  %原始图像中左右侧大脑皮层的三角片的位置（没有去除NaN）
A=nan(nVX,1); A(ivxc)=Y(ivxc,1);

Y=nan(nVX,1);
Y(A==41)=1;  %VN
Y(A==43)=2;  %SMN
Y(A==38)=3;  %DAN
Y(A==44)=4;  %VAN
Y(A==42)=5;  %LN
Y(A==39)=6;  %FPN
Y(A==40|A==37)=7;  %DMN
YLB={'V','SM','DA','VA','L','FP','DM'}; nY=7; 
save('myGlssr.mat','Y','YLB','nY','-append');

YM=Y; % for checking
YM(A==38)=4; % for compatible colors with Yeo et al 2011
YM(A==44)=3;
YM(A==39)=-7;  %FPN,为了区别DMN
e=ft_read_cifti('empty.dtseries.nii'); e.dtseries=YM;
ft_write_cifti('./yeo_FPN_DMN_wrong',e,'parameter','dtseries'); 
clear A

%% Each Glasser's parcel belongs to which Yeo's RSNs?
O=cell(nG,1); G2Y=nan(nG,1,'single'); h=nan(nG,nY); %每个parcel属于哪个网络，每个网络包含的parcel的个数
for i=1:nG
    n=length(iXG{i});
    O{i}=zeros(n,1);
    for k=1:n
        O{i}(k)=Y(iXG{i}(k));
    end
    h(i,:)=hist(O{i},1:nY); % number of Yeo's RSN per Glasser's parcel
    [~,G2Y(i)]=max(h(i,:));  %某一个parcel有10个三角片属于网络1，但是有20个三角片属于网络2，所以找哪一个网络包含的这个parcel里的三角片个数最多，这个parcel就属于哪一个网络
end
save('myGlssr.mat','G2Y','-append');

%% How to reorder Glasser's parcels based on Yeo's RSNs (i.e., YLB order)?
iG2Y=cell(nY,1); % indices for reordering (see PLTT funtion for an example)
ibY=zeros(nY+1,1); % indices of parcels at the borders of RSNs
for i=1:nY
    iG2Y{i}=single(find(G2Y==i)); 
    ibY(i+1)=length(iG2Y{i});
end

ibY=cumsum(ibY);  %求和[1 2 3] cunmsum-<[1 3 6]

ihY=zeros(nY,1); % indices in between ibY for plot labels
for i=1:nY
    ihY(i)=ibY(i)+(ibY(i+1)-ibY(i))/2; 
end

ihY=round(ihY);
save('myGlssr.mat','iG2Y','ibY','ihY','-append');
