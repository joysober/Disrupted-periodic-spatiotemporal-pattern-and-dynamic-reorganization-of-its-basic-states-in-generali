function [PLh,PLc,PLe,esg,ITP,ssg,tsh] = parameter_detectingPSTP(PL,nt,nP,nscn)

PLh=round(PL/2)+[0 -rem(PL,2)]; % pad length for temporally extending a PSTP
PLc=round(PL/2)+(1:PL); % range of interest in an extended PSTP, matches PLh
PLe=PL+sum(PLh); % length of an extended PSTP (derivable but saved/read)
esg=nt-PL+1; % 241 end segment/scan when correlating a PSTP template

ssg=ones(nP,1); ssg(2:end)=PL; % starting segment/scan ~
ITP=cell(nP,1); % all possible initial-segments/sbj for robust pstpdetection
for ip=1:nP
    for i=1:nscn
       a=(i-1)*nt+(ssg(ip):esg); 
	   ITP{ip}=[ITP{ip}; single(a')]; 
    end
end
%%%The above shows all participants traversing from all time points as starting points
tsh=floor(PL/4); % max timeshift when comparing PSTPs or any 2 templates
ssg=ones(nP,1); ssg(2:end)=PL; % starting segment/scan ~






