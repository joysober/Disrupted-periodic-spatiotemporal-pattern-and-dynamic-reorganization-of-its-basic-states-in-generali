function AllVolumeBrain = y_bandpass_Surf(Data, ...
    ASamplePeriod, ALowPass_HighCutoff, AHighPass_LowCutoff, ...
    AAddMeanBack, ...
    AMaskFilename,...
    CUTNUMBER)
%Ideal Band pass filter
%Algorithm re-written to call y_IdealFilter by YAN Chao-Gan
% FORMAT [AllVolumeBrain, GHeader] = y_bandpass_Surf(InFile, ...
%     ASamplePeriod, ALowPass_HighCutoff, AHighPass_LowCutoff, ...
%     AAddMeanBack, ...
%     AMaskFilename,...
%     AResultFilename,...
%     CUTNUMBER)
% Use Ideal rectangular filter to filter a time series
% Input:
% 	InFile              The input surface time series file
% 	ASamplePeriod		TR, or like the variable name
% 	ALowPass_HighCutoff			low pass, high cutoff of the band, eg. 0.08
% 	AHighPass_LowCutoff			high pass,  low cutoff of the band, eg. 0.01
%	AAddMeanBack			'Yes' or 'No'. 	if yes, then add the mean back after filtering
% 	AMaskFilename		the mask file name, only compute the point within the mask
%	AResultFilename		the output filename 
%   CUTNUMBER           cut the data into pieces if small RAM memory e.g. 4GB is available on PC. It can be set to 1 on server with big memory (e.g., 50GB).
%                       default: 10
% Output:
%	AllVolumeBrain      The filted data
%   GHeader             The GIfTI Header
%	 Filtered data saved as AResultFilename
%-----------------------------------------------------------
% Inherited from y_bandpass.m
% Revised by YAN Chao-Gan 181117.
% Key Laboratory of Behavioral Science and Magnetic Resonance Imaging Research Center, Institute of Psychology, Chinese Academy of Sciences, Beijing, China
% ycg.yan@gmail.com


if ~exist('CUTNUMBER','var')
    CUTNUMBER = 10;
end

tic;
% fprintf('\nIdeal rectangular filter:\t"%s"', InFile);

% GHeader=ft_read_cifti(InFile);
% AllVolume=GHeader.dtseries;
% location=find(~isnan(AllVolume));
% location_nan=find(isnan(AllVolume(:,1)));
% AllVolume(location_nan,:)=[];
% RAWdata=AllVolume;

[nDimVertex nDimTimePoints]=size(Data);

fprintf('\nLoad mask "%s".\n', AMaskFilename);
if ~isempty(AMaskFilename)
    MaskData=gifti(AMaskFilename);
    MaskData=MaskData.cdata;
    if size(MaskData,1)~=nDimVertex
        error('The size of Mask (%d) doesn''t match the required size (%d).\n',size(MaskData,1), nDimVertex);
    end
    MaskData = double(logical(MaskData));
else
    MaskData=ones(nDimVertex,1);
end
MaskDataOneDim=reshape(MaskData,1,[]);


% First dimension is time
Data=Data';
Data=Data(:,find(MaskDataOneDim));

fprintf('\n\t Band Pass Filter working.\tWait...');

%Remove the mean. YAN Chao-Gan 100401.
theMean=mean(Data);
Data=Data-repmat(theMean,[nDimTimePoints,1]);

SegmentLength = ceil(size(Data,2) / CUTNUMBER);
for iCut=1:CUTNUMBER
    if iCut~=CUTNUMBER
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(Data,2);
    end
    
    Data(:,Segment) = y_IdealFilter(Data(:,Segment), ASamplePeriod, [AHighPass_LowCutoff, ALowPass_HighCutoff]);
    
    fprintf('.');
end

% Add the mean back after filter.
if strcmpi(AAddMeanBack, 'Yes')
    Data=Data+repmat(theMean,[nDimTimePoints,1]);
end

% Get the brain back
AllVolumeBrain = zeros(nDimTimePoints,length(MaskDataOneDim));
AllVolumeBrain(:,find(MaskDataOneDim)) = Data;
AllVolumeBrain = AllVolumeBrain';

%Save all images to disk
fprintf('\n\t Band pass filter finished.\n\t');
toc;
