function [ x ] = solve_equation( a,b,c )
%SOLVE_EQUATION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
delt=b^2-4*a*c;
if(delt<0)
    x=-100;
    disp('���ε�ʽ�޽�');
elseif (delt==0)
    disp('��Ψһ��');
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

