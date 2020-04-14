function [ new_model ] = obj_linear( model,L,ThPimin,ThPimax,N,T )
%OBJ_LINEAR 此处显示有关此函数的摘要
%   此处显示详细说明
variation_num=size(model.ctype,2);
secondary_vari_location=find(diag(model.H)~=0);
secondary_vari_num=size(secondary_vari_location,1);
new_model_variation_num=variation_num+secondary_vari_num;
Aineq=[];
bineq=[];
si1=1:T:N*T;
for i=1:secondary_vari_num
    k=ceil(i/T);
    for j=0:L
        pil=ThPimin(k)+(j/L)*(ThPimax(k)-ThPimin(k));
        gradiant=sparse(1,new_model_variation_num);
        gradiant(1,i)=-0.5*model.H(secondary_vari_location(i),secondary_vari_location(i))*(pil)^2;
        gradiant(1,secondary_vari_location(i))=model.H(secondary_vari_location(i),secondary_vari_location(i))*(pil);
        gradiant(1,variation_num+i)=-1;
%          if(isempty(find(i==si1)))
            gradiant(1,2*N*T+i)=-0.5*model.H(secondary_vari_location(i),secondary_vari_location(i))*(pil)^2;
%          end
        Aineq=[Aineq;gradiant];
        bineq=[bineq;0];
    end
end
ui1=1:T:N*T;
si1=si1+2*N*T;
new_model.H=[];
% alpha=model.f(1:secondary_vari_num,1);
% model.f(1:secondary_vari_num,1)=0;
% model.f(secondary_vari_location,1)=0;
% coefficient=model.f(si1);
% model.f(2*N*T+1:3*N*T)=model.f(2*N*T+1:3*N*T)-alpha;
% model.f(si1)=coefficient;
new_model.f=[model.f;ones(secondary_vari_num,1)];
new_model.Aineq=[model.Aineq,sparse(size(model.Aineq,1),secondary_vari_num)];
new_model.bineq=model.bineq;
new_model.Aineq=[new_model.Aineq;Aineq];
new_model.bineq=[new_model.bineq;bineq];
new_model.Aeq=[model.Aeq,sparse(size(model.Aeq,1),secondary_vari_num)];
new_model.beq=model.beq;
new_model.lb=[model.lb;-inf*ones(secondary_vari_num,1)];
new_model.ub=[model.ub;inf*ones(secondary_vari_num,1)];
add_type(1,1:secondary_vari_num)='C';
new_model.ctype=strcat(model.ctype,add_type);
end

