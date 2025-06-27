function x_R = fit_extend_R(PL,fun_DMN,fun_TPN,D)
%%%%D is the mean plus three times the standard deviation of the mean
xppo = PL:0.0001:PL+5;
y1 = fun_DMN(xppo);
y2 = fun_TPN(xppo);
y_diff = y2-y1;  %%%TPN-DMN
[min_va,lo_min]=min(abs(y_diff));

if min_va~=abs(y_diff(end)) && (abs(y1(lo_min))<D)
    x_R = PL+0.0001*(lo_min-1);   
else
    xppo = PL:0.0001:PL+10;
    y1 = fun_DMN(xppo);
    y2=fun_TPN(xppo);
    y_diff = y2-y1;  %%%TPN-DMN
    [~,lo_min]=min(abs(y_diff));
    if abs(y1(lo_min))<D
        x_R = PL+0.0001*(lo_min-1);   %
    else
        x_R=[];
    end
end
end



