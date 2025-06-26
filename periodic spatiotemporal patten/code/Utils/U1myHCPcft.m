clear; clc
%划分皮层和皮层下区域
%%
addpath('./cifti-matlab-master');
e=ft_read_cifti('empty.dtseries.nii');
es=e.brainstructure; e1=e.dtseries; en=e.brainstructurelabel';

nVX=size(es,1); 
ivx=single(find(~isnan(e1))); 
nvx=length(ivx);
nVXC=find(es==3);  %左侧伏隔核，此时的位置是在有nan的情况下
nVXC=nVXC(1)-1;
es1=es(ivx);  %不为nan的位置对应的值,相当于从原始图像的数值中去掉NAN
nvxc=find(es1==3);   %%左侧伏隔核，此时的位置是在没有nan的情况下
nvxc=nvxc(1)-1;
ivxc=ivx(1:nvxc);   %只有值为1和2 的位置，这两个值分别对应左侧大脑皮层和右侧大脑皮层
save('myHCPcft.mat','nVX','ivx','nvx','nVXC','ivxc','nvxc');

%%
sc{1}=[10 11];   nmsc{1}='Cerebellum';     nmsc1{1}='C'; 
sc{2}=[20 21];   nmsc{2}='Thalamus';       nmsc1{2}='T';
sc{3}=[14 15];   nmsc{3}='Hippocampus';    nmsc1{3}='H';
sc{4}=[5 6];     nmsc{4}='Amygdala';       nmsc1{4}='A';
sc{5}=[7 12 13 16 17]; nmsc{5}='Brainstem & Deep brain';  nmsc1{5}='BD';  
sc{6}=[8 9 18 19 3 4]; nmsc{6}='Striatum'; nmsc1{6}='S';
sc=sc'; nmsc=nmsc'; nmsc1=nmsc1';

ep1=e.pos(ivx,:);   %pos里左右大脑皮层的位置为NaN,皮层下的（小脑、丘脑、海马旁回、杏仁核、脑干和深部小脑）的位置有值
                    %ep1其实是原始图像（empty.dtseries.nii）中值不为NAN的地方，有点像坐标，但不知道是矩阵坐标还是皮层坐标？？？
                    %感觉像是矩阵坐标，因为是整数，surface坐标有小数
N=length(sc); ivxsc=cell(N,1); pvxsc=cell(N,1); 
for i=1:N
    ind=[]; 
    for j=1:length(sc{i})
        ind=[ind; single(find(es1==sc{i}(j)))];
    end
    ivxsc{i}=ind; %找到每个区域的位置（比如小脑、丘脑等）
    pvxsc{i}=single(ep1(ind,:));  %对应的矩阵坐标
end
save('myHCPcft.mat','sc','nmsc','nmsc1','ivxsc','pvxsc','-append')

% for i=1:N % for checking
%     m=nan(nvx,1); m(ivxsc{i})=1; M=nan(nVX,1); M(ivx)=m; e.dtseries=M;
%     ft_write_cifti(['./U3myPrcls/HCPMSK_' num2str(i) nmsc1{i}],e,...
%         'parameter','dtseries');
% end
%% 加入皮层
nmrgn=[{'Cortex'}; nmsc]; nmrgn1=[{'Cx'}; nmsc1]; 
irgn=[single(1:nvxc)'; ivxsc]; nrgn=N+1; 
save('myHCPcft.mat','nmrgn','nmrgn1','irgn','nrgn','-append');

%% 将各区域的左右分开
sclr{1}=10; nmsclr1{1}='C-L'; 
sclr{2}=11; nmsclr1{2}='C-R'; 
sclr{3}=20; nmsclr1{3}='T-L'; 
sclr{4}=21; nmsclr1{4}='T-R';
sclr{5}=14; nmsclr1{5}='H-L'; 
sclr{6}=15; nmsclr1{6}='H-R';
sclr{7}=5; nmsclr1{7}='A-L'; 
sclr{8}=6; nmsclr1{8}='A-R';
sclr{9}=[7 12 16]; nmsclr1{9}='BD-L';  %脑干不分左右，但是还是对脑干进行了分组（以矩阵坐标的第一个值是否为正）
                                       %如果为正，则把这个脑干归为BD-R，如果为负，则归为BD-L。其实也许在矩阵坐标里面，正为右侧，负为左侧
sclr{10}=[7 13 17]; nmsclr1{10}='BD-R'; 
sclr{11}=[8 18 3]; nmsclr1{11}='S-L';
sclr{12}=[9 19 4]; nmsclr1{12}='S-R';
sclr=sclr'; nmsclr1=nmsclr1';

ep1=e.pos(ivx,:); 
N=length(sclr); ivxsclr=cell(N,1); pvxsclr=cell(N,1); 
for i=1:N
    ind=[]; 
    for j=1:length(sclr{i})
        s=sclr{i}(j); I=find(es1==s);
        if s==7, x=ep1(:,1); 
            if i==9, I(x(I)>=0)=[]; else, I(x(I)<0)=[]; end
        end
        ind=[ind; single(I)]; 
    end
    ivxsclr{i}=ind; pvxsclr{i}=single(ep1(ind,:));
end
save('myHCPcft.mat','sclr','nmsclr1','-append')

ivxclr={single(find(es1==1)); single(find(es1==2))};
irgnlr=[ivxclr; ivxsclr]; nmrgnlr1=[{'Cx-L'}; {'Cx-R'}; nmsclr1]; 
save('myHCPcft.mat','irgnlr','nmrgnlr1','-append');
