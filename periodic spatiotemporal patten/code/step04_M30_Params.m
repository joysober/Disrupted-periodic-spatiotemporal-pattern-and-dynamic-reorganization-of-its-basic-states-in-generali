
%% Setting paths, parameters & indices
%% 
clear; clc;
nmD='HC'; % dataset name 
p2prep=... % path to (pth2/p2) preproc params, see ../../All_Pre~/P0_Pre~.m 
    ['../All_Preproc_Data/Prep_Params_' nmD '.mat'];

%% For saving & reading path-names
p2p=['../All_Preproc_Data/Params_PSTP_' nmD]; % p2 SAVE SETTINGS HERE >> SHOULD BEGIN WITH Params_

load(p2prep,'p2BGWCR','p2bgwcr','nsbj','nscn'); a=strfind(p2prep,'/'); 
p2B=[p2prep(1:a(end)) p2BGWCR]; % pth2 parcellated scans
p2bpr=cell(nsbj,nscn); % pth2 preproc-ed grayordinate scans

for is=1:nsbj, for i=1:nscn
        p2bpr{is,i}=[p2prep(1:2) p2bgwcr{is,i}(3:end)];  end; end

nP=2; %The first one is the PSTP, and the second one is PSTP2, which is recalculated after regressing out the original PSTP.
cth=cell(nP,1); % correlation threshold in the main algorithm
nitr=15; % max #iterations in the main algorithm
ncth1=3*ones(nP,1); if nP>=4, ncth1(4:end)=3; end % #iters with lower cth
%%%%%%%%%%%%%%%%%%¸Ä%%%%%%%%%%%%%%%%%5
% cth(1:min(nP,3))={[0.1 0.2]}; cth(4:2end)={[0.1 0.2]};
cth(1:min(nP,3))={[75 85]}; cth(4:end)={[75 85]};
load(p2prep,'nt','nscng');
ITPstp=50; % step to show progress of algorithm when running ITP times
%% For functional connectivity (FC)
load('myGlssr.mat',... % indices to reorder Glasser's parcels based on ...
   'nY','G2Y','ibY','iG2Y','YLB'); %  Yeo's 7 RSNs, see [p2u 'U2myGlssr.m']
fcbn=-1:0.01:1; % histogram bins for FC
fcth=0.1; % FC-values with magnitude > fcth are found after regressing QPPs
%% Clearing extras, loading needed from prep.mat & saving all set here
clear a i ip is nmD  p2BGWCR p2bgwcr
load(p2prep,'nT','nTg','ivx','nvx','nVX','nX','tres','INFO'); % nsbj,nscn,nt
save(p2p);
