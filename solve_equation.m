function [ x ] = solve_equation( a,b,c )
%SOLVE_EQUATION 此处显示有关此函数的摘要
%   此处显示详细说明
delt=b^2-4*a*c;
if(delt<0)
    x=-100;
    disp('二次等式无解');
elseif (delt==0)
    disp('有唯一解');
    x=-b/(2*a);
else
    x_1=(-b+sqrt(delt))/(2*a);
    x_2=(-b-sqrt(delt))/(2*a);
    if(x_1>=0&&x_1<=1)
        x=x_1;
    
    elseif(x_2>=0&&x_2<=1)
        x=x_2;
    else
        x=-100;
    end
end

end

