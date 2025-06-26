
%Extract cerebrospinal fluid and white matter signals for each subject.
clc;clear;

inputDir='../confounds_timeseries';
outputDir='../confounds_WM_CSF_Global';
listing=dir(inputDir);
listing(1:2)=[];
num_day=length(listing);
filename={'WM','CSF','Global'};
num_signal=length(filename);

for d=1:num_day

    session=[inputDir filesep listing(d).name];
    list_se=dir(session);
    list_se(1:2)=[];
    num_se=length(list_se);

    for s=1:num_se
        for j=1:num_signal
            myRoot1 = [outputDir filesep listing(d).name filesep list_se(s).name filesep filename{j}];
            if ~isfolder(myRoot1) 
                mkdir(myRoot1);
            end
        end
    end
    clear myRoot1 s j

    for s=1:num_se
        G=[inputDir filesep listing(d).name filesep list_se(s).name filesep '*.tsv'];
        list=dir(G);
        num=length(list);

        for i=1:num

            data=importdata([inputDir filesep listing(d).name filesep list_se(s).name filesep list(i).name]);
            CSF=data.textdata(2:end,5);
            WM=data.textdata(2:end,9);
            GL=data.textdata(2:end,1);

            save_CSF=[outputDir filesep listing(d).name filesep list_se(s).name filesep filename{2} filesep list(i).name(1:end-4) '_CSF.txt'];
            save_WM=[outputDir filesep listing(d).name filesep list_se(s).name filesep filename{1} filesep list(i).name(1:end-4) '_WM.txt'];
            save_GL=[outputDir filesep listing(d).name filesep list_se(s).name filesep filename{3} filesep list(i).name(1:end-4) '_Global.txt'];

            writecell(CSF,save_CSF);
            writecell(WM,save_WM);
            writecell(GL,save_GL);
        end

    end
end





