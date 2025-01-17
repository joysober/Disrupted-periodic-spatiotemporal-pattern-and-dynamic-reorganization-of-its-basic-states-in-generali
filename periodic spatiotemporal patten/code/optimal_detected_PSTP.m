function [PSTP,TMX,C,ITER,MET,ring,threshold_parceil_np1] = optimal_detected_PSTP(TMPL,ITP1,TMXTMPL,CTMPL,ITR,tres,threshold_parceil_np1_transition)

%% Selecting PSTP
PSTP=TMPL{ITP1}; %Select the most suitable starting point and the corresponding final template
TMX=TMXTMPL{ITP1}; %Select the most suitable starting point and the time point corresponding to the local maximum value
C=CTMPL(ITP1,:); %correlation timecourse
ITER=ITR(ITP1);  %the number of iterations
MET=[median(C(TMX)) median(diff(TMX))*tres length(TMX)];
ring=ITP1;
threshold_parceil_np1 = threshold_parceil_np1_transition(ITP1,:);
end

