function [xp,yc,pc,c] = plot_polyfit(xp,yp,lowlevel,highlevel)

c=[];
r2_all = zeros(2,highlevel-1);  %The order starts from 2, so here we subtract 1; The second line stores labels (whether overfitting occurs)
for i=lowlevel:highlevel        %Find the most suitable fit within the range of 4-12, and use the sum of squared errors to determine if it is appropriate
    %     y2=polyfit(xp,yp,i);
    [y2,lable]=polyfit_modify_cjx(xp,yp,i);
    Y=polyval(y2,xp);%Calculate the value of the fitting function at x.
    [r2,~]=rsquare(yp,Y);
    r2_all(1,i-1)=r2;      
    r2_all(2,i-1)=lable;
    if r2>0.99  %Sum of squared errors accuracy
        c=i;
        break;
    end
end
clear i
if isempty(c)
    [~,loca]=max(r2_all(1,:));  
    for j=loca:-1:1
        if r2_all(2,j)~=1 %If the label corresponding to the maximum value of r2 is 1, it indicates a warning. Then use the previous label of loca, and so on until there is no warning
            c=j+1;    %j is the position, and adding 1 is because the first position is fitted from the 2nd order onwards
            break;
        end
    end
    
end
pc=polyfit(xp,yp,c);
yc=polyval(pc,xp);
pc=vpa(poly2sym(pc),c);

end

