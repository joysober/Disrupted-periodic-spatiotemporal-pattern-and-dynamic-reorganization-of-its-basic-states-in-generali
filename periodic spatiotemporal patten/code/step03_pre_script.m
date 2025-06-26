%% processing
clc;clear;
nmD='HC';
p2prep=['../All_Preproc_Data/Prep_Params_' nmD '.mat'];
addpath(genpath(p2prep));
load(p2prep,'nsbj');
prfxin='';
ISin=nsbj;
P1_Preproc_Run_HC(nmD,prfxin,ISin);
 
