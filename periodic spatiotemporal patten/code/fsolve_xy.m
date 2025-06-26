function [x_pointf,y_point] = fsolve_xy(x1,x2,PL)

%%%%At least two or three intersection points will be output
%%%%The two intersection points are when there is no intersection point between the middle 1 and PL; 
%%Three intersection points refer to having one, two, or three intersection points in the middle
syms func1(x)
func1(x)=x1;
syms func2(x)
func2(x)=x2;
p=sym2poly(func1-func2);
x_pointf=roots(p);

      for realx=1:length(x_pointf)
          returnValue=isreal(x_pointf(realx));
          if returnValue==0
              x_pointf(realx)=-100;    %Set a complex number to a number that cannot be an intersection point
          end
      end
      x_pointf(x_pointf==-100)=[];
%%%%%%%%%Only keep the two intersection points in the middle and next to it
x_pointf=single(x_pointf);
x_pointf=floor(x_pointf*10000)/10000;
x_pointf=unique(x_pointf);

x_pointf=sort(x_pointf);

med1=find(x_pointf<1 );
med2=find(x_pointf>PL);
if ~isempty(med1) && length(med1)~=1
    x_pointf(1:med1(end)-1)= -100;
end
if ~isempty(med2) && length(med2)~=1
    x_pointf(med2(1)+1:end)= -100;
end
x_pointf(x_pointf==-100)=[];   %%%%Ensure that only one intersection point on each side of 1 and PL is preserved

num_med = length(find((x_pointf>1 & x_pointf<PL)==1));

x_point_med = x_pointf;
x_point_med(x_point_med<1)=-100;
x_point_med(x_point_med>PL)=-100;
x_point_med(x_point_med==-100)=[];
if num_med==3
    x_pointf = x_point_med;
elseif num_med==2
    if ~isempty(find(x_pointf<1)) && ~isempty(find(x_pointf>PL))
        %%There must be no intersection on the left hand side, because the first condition (DMN(1) is greater than TPN(1) has already been decided)
        x_pointf(end)=[];  
    end
    
end
y_point=func1(x_pointf);
end

