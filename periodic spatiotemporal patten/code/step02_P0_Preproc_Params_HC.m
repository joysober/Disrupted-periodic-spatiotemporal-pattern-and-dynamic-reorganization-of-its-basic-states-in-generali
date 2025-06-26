%% This code is for presetting parameters in advance
%% The original code comes from Behnaz Yousefi,NeuroImage, 2016
%% Modified by Junxia Chen, Email: joysober@163.com.

%% Setting preprocessing parameters, paths, etc, compatible with QPPv0420
%% Dataset name, sbj IDs, etc
clear; clc; nmD='HC'; % dataset name
p2prep=... % path-name to (pth2/p2) SAVE all params, paths, etc, set here
    ['../All_Preproc_Data/Prep_Params_' nmD '.mat']; %% >> KEEP THIS FORMAT OR CHANGE P1_P~.m

p2data='../data/';%note that this is require to modify   
addpath(genpath(p2data)); % pth2 cifti-/nii-related utilities
ID=importdata([p2data 'ID_' nmD '.txt']); 

nsbj=size(ID,1); % IDs & consequent #sbj
nD=1; nPD=2; % #days, #scans/day (>>HCP-BASED)   nD: number of scanning days; nPD: number of sessions per day.
iD=1:nD;   %
iPD=1:nPD; % indices of day&scan to process (all scans iPD=1:nPD;)

nscn=length(iD)*length(iPD); % nscn = number of scanning days per subject * number of sessions per day.
nscng=nsbj*nscn; % conseq. #scans/grp, subjects*scans=3*2

p2u='./Utils/';%note that this is require to modify  
addpath(genpath(p2u)); % pth2 cifti-/nii-related utilities

%% Pth2 read scans & WM&CSF signals (nuisance regressors) & get params/scan
p1={'day01/'};
p2={'session01/','session02/'};
nmb={'_task-restclose_space-fsLR_den-91k_bold.dtseries.nii','_task-restclose_space-fsRL_den-91k_bold.dtseries.nii'}; 

p2confound='../confounds_WM_CSF_Global/';%note that this is require to modify 
addpath(genpath(p2confound)); % pth2 cifti-/nii-related utilities
nmrg={'WM','CSF'};
nmc='_task-restclose_desc-confounds_timeseries_';

nt=250; tres=2; % another way to set params/scan    %时间点和TR
load('myHCPcft.mat','nVX','ivx','nvx'); % see [p2u 'U1myHCPcft.m']

nT=nt*nscn;  
nTg=nt*nscng; % conseq. # all-timepoints/sbj 

%% Directory to save preprocessed scans in grayordinates
d2prep={['../All_Preproc_Data/WCR_' nmD '/'],... % WM&CSF Regressed scans
    ['../All_Preproc_Data/GWCR_' nmD '/']};  % Gray-Matter(GM)&WM&CSF Regressed scans
for i=1:2
    if ~exist(d2prep{i},'dir')
        mkdir(d2prep{i});
    end
end

%% Putting path-names to read&save grayordinate scans in cell arrays
c=cell(nsbj,nscn);
p2b=c; p2rg=c; p2bwcr=c; p2bgwcr=c; IDF=c; clear c
for is=1:nsbj
    iscn=1;
    %     IDS=num2str(ID(is));
    IDS=ID{is};

    for id=iD   
        for ipd=iPD  
            pb=[p1{iD} p2{ipd}];
            nmbb=[IDS nmb{ipd}];  %%
            p2b{is,iscn}=[p2data pb nmbb]; % pth2 read scans
            prg=[p2confound p1{iD} p2{ipd}];
            for i=1:length(nmrg)
                p2rg{is,iscn}{i}=[prg nmrg{i} '/' IDS nmc nmrg{i} '.txt'];
            end
            % pth2 read regressors
            a=[IDS '_day_' num2str(id,'%02d') '_session_' num2str(ipd,'%02d')]; 
            IDF{is,iscn}=a;
            a=[a '.mat'];
            p2bwcr{is,iscn}=[d2prep{1} a];
            p2bgwcr{is,iscn}=[d2prep{2} a]; % p2 save
            iscn=iscn+1;
        end
    end
end
clear pb nmbb a ID IDS d2prep i is  ...  
    id ipd iscn nD nPD
INFO.IDF=IDF; INFO.iD=iD; INFO.iPD=iPD; clear IDF iD iPD

%% Bandpass filter & zeropad size
ALowPass_HighCutoff=0.1;    %the cut-off frequence of lowpass
AHighPass_LowCutoff=0.01;   %the cut-off frequence of highpass
AAddMeanBack='no';    %no need mask when filter
AMaskFilename=[];     
%% Indices of Glasser's parcellation
load('myGlssr.mat','ixG','iXG'); nX=length(ixG); % see [p2u 'U2myGssr.m']
%% Pth2 save parcellated scans
p2BWCR=['B_WCR_' nmD '.mat']; p2BGWCR=['B_GWCR_' nmD '.mat'];
%% Clearing extras & SAVING all
clear nmD  
save(p2prep);

