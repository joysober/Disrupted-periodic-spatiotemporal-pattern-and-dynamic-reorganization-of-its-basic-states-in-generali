function [x0,y1] = inflection_point(fun,PL,x)
    d2 = diff(fun,x,1); %Find the first derivative of function 
    x0 = solve(d2,0); %Find the point where the first derivative is 0 (i.e. the extremum point)
    x0(x0<1)=0;
    x0(x0>PL)=0;
    for rr=1:length(x0)
       complex=isreal(x0(rr));
            if complex==0
                x0(rr)=0;
            end
    end
    x0(x0==0)=[];  
    y0 = subs(fun,x,x0); %The y-value at the extremum point
    para=matlabFunction(y0);
    y1=vpa(para,3);
end

