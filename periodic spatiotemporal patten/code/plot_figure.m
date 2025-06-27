function plot_figure(xp,yp_DMN,yp_TPN,fun_DMN,fun_TPN,PL,...
    c1,c2, period,pathpng,pathpngz,ip,is,yc_DMN,yc_TPN,label_plot,x_point,x_L,x_R,y_point)

set(0,'DefaultFigureVisible', 'off');
%  
if label_plot==1   %The situation where there is no intersection point between 1 and PL

    y_DMN_L = fun_DMN(x_L);y_DMN_R = fun_DMN(x_R);
    y_TPN_L = fun_TPN(x_L);y_TPN_R = fun_TPN(x_R);
    plot(x_L,y_DMN_L,'ko');hold on;plot(x_R,y_DMN_R,'ko');hold on;
    plot(x_L,y_TPN_L,'ko');hold on;plot(x_R,y_TPN_R,'ko');hold on;
    
    
elseif label_plot==2    %%%The intersection point in the middle is greater than 3
    
    plot(x_point(4),fun_DMN(x_point(4)),'ko');hold on;plot(x_point(2),fun_DMN(x_point(2)),'ko');hold on;
    
elseif label_plot==3  %%%There are two intersection points in the middle, but the fitted curve has two or three intersection points
    
    plot(x_L,fun_DMN(x_L),'ko'); hold on;  %%%
    plot(x_R,fun_DMN(x_R),'ko'); hold on;
    %%%%%%%%%%%%%%%Draw the extended area (intersection on both sides)%%%%%%%%%%%%%%%%%%%%
    if  x_L<1
        xpne=x_L-1:0.01:1;  plot(xpne,fun_DMN(xpne),'b:');hold on; %
        plot(xpne,fun_TPN(xpne),'r:');hold on;
    end
    
    if x_R>PL
        xppo=PL:0.01:x_R+1;       
        plot(xppo,fun_DMN(xppo),'b:');hold on;
        plot(xppo,fun_TPN(xppo),'r:');hold on;
    end
    %
elseif label_plot==4     %%There is only one intersection point in the middle, but the fitted curve has three intersection points (one in the middle and the other two outside of 1~PL)
    
%     plot(x_point,y_point,'ko'); hold on;
    
    plot(x_L,fun_DMN(x_L),'ko'); hold on; 
    plot(x_R,fun_DMN(x_R),'ko'); hold on;
    %%%%%%%%%%%%%%%Draw the extended area (intersection on both sides)%%%%%%%%%%%%%%%%%%%%
    
    xpne=x_L:0.01:1;  plot(xpne,fun_DMN(xpne),'b:');hold on
    plot(xpne,fun_TPN(xpne),'r:');hold on
    
    xppo=PL:0.01:x_R(end);
    plot(xppo,fun_DMN(xppo),'b:');
    hold on
    plot(xppo,fun_TPN(xppo),'r:');
    
elseif label_plot==5   
 %There is only one intersection point in the middle, but the curve may have two intersection points (except for the middle, there is one outside) or one (which is the intersection point itself)
    plot(x_L,fun_DMN(x_L),'ko'); hold on;  
    plot(x_R,fun_DMN(x_R),'ko'); hold on;
    %%%%%%%%%%%%%%%Draw the extended area (intersection on both sides)%%%%%%%%%%%%%%%%%%%%
    xpne=x_L-0.5:0.01:1;  plot(xpne,fun_DMN(xpne),'b:');hold on; 
    plot(xpne,fun_TPN(xpne),'r:');hold on;
    
    xppo=PL:0.01:x_R+0.5;      
    plot(xppo,fun_DMN(xppo),'b:');hold on;
    plot(xppo,fun_TPN(xppo),'r:');hold on;
end


%%%%%%%%%%%%PLOT%%%%%%%%%%%%
plot(xp,yp_DMN,'b*');hold on;
fplot(fun_DMN,[1,PL],'b');hold on
plot(xp,yp_TPN,'r*');hold on
fplot(fun_TPN,[1,PL],'r');hold on
grid on

xlabel('time');ylabel('signal value');
title(['c1=',num2str(c1),'; c2=',num2str(c2),'; period=',num2str(ceil(period))]);
%%%
saveas(gcf,[pathpng 'PSTP' num2str(ip) '_sub' num2str(is,'%03d') '.png']);
close;
figure;yc_DMN=zscore(yc_DMN);yc_TPN=zscore(yc_TPN);plot(xp,yc_DMN,'b'); hold on
plot(xp,yc_TPN,'r');grid on
xlabel('time');ylabel('zscore signal value');legend('DMN','TPN')
saveas(gcf,[pathpngz 'PSTP' num2str(ip) '_sub' num2str(is,'%03d') '_zscore.png']);
close;
end











