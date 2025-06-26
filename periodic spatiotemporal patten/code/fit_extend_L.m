function x_L = fit_extend_L(fun_DMN,fun_TPN,D)
%Expand the left side
xpne=-2:0.0001:1;
y3=fun_DMN(xpne);
y4=fun_TPN(xpne);
y_diff = y3-y4;   %%%DMN-TPN
[min_va,lo_min] = min(abs(y_diff));
if min_va~= abs(y_diff(1)) && (abs(y3(lo_min)) < D)
    x_L = -2+0.0001*(lo_min-1);     
else
    xpne=-10:0.0001:1;
    y3=fun_DMN(xpne);
    y4=fun_TPN(xpne);
    y_diff = y3-y4;   %%%DMN-TPN
    [~,lo_min] = min(abs(y_diff));
    if abs(y3(lo_min))< D
        x_L = -10+0.0001*(lo_min-1);     
    else
        x_L=[];
    end
end
end

