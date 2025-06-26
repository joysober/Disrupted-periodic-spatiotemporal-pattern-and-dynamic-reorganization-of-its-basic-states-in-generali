
function P1_Preproc_Run_HC(nmD,prfxin,ISin)

%% Preprocessing
% Per scan the following steps are performed: demeaning & filtering,
% regressing WM&CSF & parcellating , regressing GM&WM&CSF & parcellating
% > after regressing WM&CSF or GM&WM&CSF, each scan is saved seperately
% > after parcellating, all scans are saved in a mat file, seperately for
% WCR & GWCR
% > Preprocessing of subjects can be divided across cpus/servers and the
% cell arrays containing parcellated scans later combined by P2_Bcmbn
% IS: indices of sbjs to be preprocessed per cpu/server, [] for all sbjs
% prfx: temporary prefix to save the mat file of IS subjects  
% nmD: dataset name
% > Example: P1_Preproc_Run_HCP('HC','',1:4)

%%
p2prep='../All_Preproc_Data/';
addpath(genpath(p2prep));
load([p2prep 'Prep_Params_' nmD '.mat'],'nsbj','nscn','ivx','nt','nX','ixG',...
    'p2b','p2rg','p2bwcr','p2bgwcr','p2BWCR','p2BGWCR','p2u','INFO',...
    'tres','ALowPass_HighCutoff','AHighPass_LowCutoff','AAddMeanBack','AMaskFilename');  

addpath(genpath(p2u));
prfx=''; IS=1:nsbj;
if nargin>1
    prfx=prfxin;
    IS=ISin;
end

BWCR=cell(nsbj,nscn); BGWCR=cell(nsbj,nscn);
for is=1:IS
    for iscn=1:nscn
        fprintf('Sbj%dScn%d ',is,iscn)
        tic;
        
        b=ft_read_cifti(p2b{is,iscn});
        Data=b.dtseries(ivx,:);
        % Perform mean removal. This step is commented out here because mean removal is already performed during filtering.
%       Data=b-repmat(mean(b,2),1,nt); 
        y = y_bandpass_Surf(Data, ...
               tres, ALowPass_HighCutoff, AHighPass_LowCutoff, ...
               AAddMeanBack, ...
               AMaskFilename);
        fprintf('\nIdeal rectangular filter:\t"%s"', p2b{is,iscn});
        bft=y-repmat(mean(y,2),1,nt);    
        bft=bft';
        rg0=[dlmread(p2rg{is,iscn}{1}) dlmread(p2rg{is,iscn}{2})];   
        rg=zeros(nt,2,'single');
        for i=1:2   %Perform filtering on the nuisance (confounding) variables.
            y0= y_bandpass_Surf(rg0(:,i)', ...
                 tres, ALowPass_HighCutoff, AHighPass_LowCutoff, ...
                 AAddMeanBack, ...
                 AMaskFilename,1);
            rg(:,i)=y0;  
        end
        rg=zscore(rg);   
        beta=(rg'*rg)\rg'*bft;   
        bpr=(bft-rg*beta)';   %Residuals after regressing out white matter and cerebrospinal fluid signals
        save(p2bwcr{is,iscn},'bpr','-v7.3');
        BWCR{is,iscn}=zeros(nX,nt,'single');
        for i=1:nX
            BWCR{is,iscn}(i,:)=zscore(mean(bpr(ixG{i},:)));
        end
       
        rg=[zscore(mean(bft,2)) rg];  
        beta=(rg'*rg)\rg'*bft; bpr=(bft-rg*beta)';  
        save(p2bgwcr{is,iscn},'bpr','-v7.3');
        BGWCR{is,iscn}=zeros(nX,nt,'single');
        for i=1:nX, BGWCR{is,iscn}(i,:)=zscore(mean(bpr(ixG{i},:))); end    
        
        clear bpr bft; t=toc; t=round(t); fprintf('%d sec\n',t)
    end
end
B=BWCR; save([p2prep prfx p2BWCR],'B','INFO','-v7.3');
B=BGWCR; save([p2prep prfx p2BGWCR],'B','INFO','-v7.3');
end
