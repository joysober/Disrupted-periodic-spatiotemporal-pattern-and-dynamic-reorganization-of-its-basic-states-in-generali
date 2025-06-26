
%提取各被试的脑脊液和白质信号
%这个代码是基于已经将预处理后的所有被试的“_desc-confounds_timeseries.tsv”单独全部放进一个文件夹
%这里只有一个被试只有一个扫描，如果一个被试有多天，每天有多个session，则需要多几层文件夹
%多个session 的情况

clc;clear;

inputDir='F:\work2_相位依赖\pipeline\work_PSTP\periodic spatiotemporal patten_V2\confounds_timeseries';
outputDir='F:\work2_相位依赖\pipeline\work_PSTP\periodic spatiotemporal patten_V2\confounds_WM_CSF_Global';
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
            if ~isfolder(myRoot1) %判断路径是否存在
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





