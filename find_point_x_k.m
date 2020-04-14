function [ x_k,error ] = find_point_x_k(solution_NLP,solution_LP,Q,r,l)
%FIND_POINT_X_K 此处显示有关此函数的摘要
%   此处显示详细说明


deta_x=solution_NLP-solution_LP;
a=deta_x'*Q*deta_x;
b=solution_LP'*Q*deta_x+deta_x'*Q*solution_LP+l'*deta_x;
c=solution_LP'*Q*solution_LP+l'*solution_LP-r;
[lambda]=solve_equation(a,b,c);

if(lambda~=-100)
    
    x_k=lambda*deta_x+solution_LP;
    error=0;
    max_g=x_k'*Q*x_k+l'*x_k-r;
    
 
else
    error=1;
    disp('λ错误');
end
end

