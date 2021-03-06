function [new_model]=obj_linear_general_tangent(model,L,ThPimin,ThPimax,N,T,CET)
variation_num=size(model.ctype,2);
secondary_vari_location=find(diag(model.H)~=0);
secondary_vari_num=size(secondary_vari_location,1);
new_model_variation_num=variation_num+secondary_vari_num;
Aineq=[];
bineq=[];

J=0:L;
for j=1:L+1
    pil=reshape(repmat(ThPimin+(J(j)/L)*(ThPimax-ThPimin),1,T)',secondary_vari_num,1);
    gradiant=sparse(secondary_vari_num,new_model_variation_num);
    gradiant(1:secondary_vari_num,secondary_vari_location)=sparse(diag(model.H(secondary_vari_location,secondary_vari_location)*(pil)));
    gradiant(1:secondary_vari_num,variation_num+1:new_model_variation_num)=sparse(diag(-1*ones(1,secondary_vari_num)));
    Aineq=[Aineq;gradiant];
    bineq=[bineq;0.5*model.H(secondary_vari_location,secondary_vari_location)*power(pil,2)];
end


new_model.H=[];
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
if(CET==1)
    new_model.Q=[model.Q,sparse(size(model.Q,1),secondary_vari_num);sparse(secondary_vari_num,size([model.Q,sparse(size(model.Q,1),secondary_vari_num)],2))];
    new_model.l=[model.l;sparse(secondary_vari_num,1)];
    new_model.r=model.r;
end
end