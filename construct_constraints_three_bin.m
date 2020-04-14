function [model]=construct_constraints_three_bin(dataUC)
dataUC=readdataUC('UC_AF/10_std.mod');    
Alpha = dataUC.alpha;                           %�����鷢�纯��ϵ��Alpha--N*1����
Beta = dataUC.beta;                             %�����鷢�纯��ϵ��Beta--N*1����
Gama = dataUC.gamma;                            %�����鷢�纯��ϵ��Gama--N*1����
ThPimin = dataUC.p_low;                         %�����鷢�繦���½�--N*1����
ThPimax = dataUC.p_up;                          %�����鷢�繦���Ͻ�--N*1����
Piup = dataUC.p_rampup;                         %���������¹����Ͻ�--N*1����
Pidown = dataUC.p_rampdown;                     %���������¹����Ͻ�--N*1����
Dt = dataUC.PD;                                 %��������--T*1����
N = dataUC.N;                                   %����������--1*1����
T = dataUC.T;                                   %ʱ�����--1*1����
Spin = dataUC.spin;                             %��ת�ȱ���--T*1����
ThTime_on_off_init = dataUC.time_on_off_ini;    %�������ڳ�ʼ״̬ǰ�Ѿ�����/ͣ����ʱ��--N*1����
Thon_off_init = dataUC.p_initial;               %����������ʼ����--N*1����
ThTime_on_min = dataUC.time_min_on;             %��������С����ʱ��--N*1����
ThTime_off_min = dataUC.time_min_off;           %��������Сͣ��ʱ��--N*1����
ThCold_cost_start = dataUC.Cold_cost;           %����������������--N*1����
ThHot_cost_start = dataUC.Hot_cost;             %����������������--N*1����
ThCold_time_start = dataUC.Cold_hour;           %������������ʱ��--N*1����
Pistartup = dataUC.p_startup;                   %�����鿪������--N*1����
Pishutdown = dataUC.p_shutdown;                 %������ػ�����--N*1����



%��������
uit=sparse(N*T,1);
sit=sparse(N*T,1);
dit=sparse(N*T,1);
pit=sparse(N*T,1);
Sit=sparse(N*T,1);

%�ϲ�������x
x=[uit;sit;dit;pit;Sit];
varaible_num=length(x);%��ȡ�ܵı�������

                                        %����Լ��
%������֪��
Ui0 = full(spones(Thon_off_init));    
p0=Thon_off_init;
Ui = max(0,min(ones(N,1) * T,Ui0 .* (ThTime_on_min - ThTime_on_off_init)));                    %--N*1����
Li = max(0,min(ones(N,1) * T,(ones(N,1) - Ui0) .*  (ThTime_off_min + ThTime_on_off_init)));     %--N*1����

%% initial status
Aeq_constraints_initial_statues=[];
beq_constraints_initial_statues=[];

for i=1:N
    if(Ui(i) + Li(i) >= 1)
         for t=1:Ui(i)+Li(i)
             
             uit=sparse(1,N*T);
             
             uit(1,(i-1)*T+t)=1;
             
             constraints_initial_statues=sparse(1,varaible_num);
             constraints_initial_statues(1,1:N*T)=uit;
             
             Aeq_constraints_initial_statues=[Aeq_constraints_initial_statues;constraints_initial_statues];
             beq_constraints_initial_statues=[beq_constraints_initial_statues;Ui0(i)];
         end
    end
end
        
%% state constraint
%sit-dit=uit- uit-1
Aeq_constraint_state=[];
beq_constraint_state=-1*Ui0;

%t=1ʱ
for i=1:N
    constraint_state=sparse(1,varaible_num);

    uit=sparse(1,N*T);
    sit=sparse(1,N*T);
    dit=sparse(1,N*T);
    
    uit(1,(i-1)*T+1)=-1;
    sit(1,(i-1)*T+1)=1;
    dit(1,(i-1)*T+1)=-1;
    
    constraint_state(1,1:N*T)=uit;
    constraint_state(1,N*T+1:2*N*T)=sit;
    constraint_state(1,2*N*T+1:3*N*T)=dit;
    
    Aeq_constraint_state=[Aeq_constraint_state;constraint_state];
end
%t>1
for i=1:N
    for t=2:T
         constraint_state=sparse(1,varaible_num);
         
         uit=sparse(1,N*T);
         sit=sparse(1,N*T);
         dit=sparse(1,N*T);
         
         uit(1,(i-1)*T+t)=-1;
          uit(1,(i-1)*T+t-1)=1;
         sit(1,(i-1)*T+t)=1;
         dit(1,(i-1)*T+t)=-1;
         
         constraint_state(1,1:N*T)=uit;
         constraint_state(1,N*T+1:2*N*T)=sit;
         constraint_state(1,2*N*T+1:3*N*T)=dit;
         
         Aeq_constraint_state=[Aeq_constraint_state;constraint_state];
    end
end
beq_constraint_state=[beq_constraint_state;sparse(N*(T-1),1)];
%%  power balance constraint
Aeq_constraint_power_balance=[];
beq_constraint_power_balance=Dt(1:T-1);
pit=[];
    constraint_power_balance=sparse(T,varaible_num);
 for t=1:N     
    pit=[pit,eye(T)];
 end   
    constraint_power_balance(1:T,3*N*T+1:4*N*T)=pit;
    
    Aeq_constraint_power_balance=[Aeq_constraint_power_balance;constraint_power_balance(1:T-1,:)];

%% start_cost_constraint
%Sit>=Chot,i,t*sit
Aineq_constraint_start_cost=[];

for i=1:N
    for t=1:T
        constraint_start_cost=sparse(1,varaible_num);
        
        sit=sparse(1,N*T);
        Sit=sparse(1,N*T);
        
        sit(1,(i-1)*T+t)=ThHot_cost_start(i);
        Sit(1,(i-1)*T+t)=-1;
        
        constraint_start_cost(1,N*T+1:2*N*T)=sit;
        constraint_start_cost(1,4*N*T+1:5*N*T)=Sit;
        
        Aineq_constraint_start_cost=[Aineq_constraint_start_cost;constraint_start_cost];
    end
end
bineq_constraint_start_cost=sparse(N*T,1);

%Sit>=Ccold,i,t*[sit-sum(dit)-finint]
for i=1:N
    for t=1:T
        constraint_start_cost=sparse(1,varaible_num);
        
        sit=sparse(1,N*T);
        dit=sparse(1,N*T);
        Sit=sparse(1,N*T);
        
        sit(1,(i-1)*T+t)=ThCold_cost_start(i);
        dit(1,(i-1)*T+max(t-ThTime_off_min(i)-ThCold_time_start(i),1):(i-1)*T+t-1)=-1*ThCold_cost_start(i);
        Sit(1,(i-1)*T+t)=-1;
        
        constraint_start_cost(1,N*T+1:2*N*T)=sit;
         constraint_start_cost(1,2*N*T+1:3*N*T)=dit;
          constraint_start_cost(1,4*N*T+1:5*N*T)=Sit;
          
          Aineq_constraint_start_cost=[Aineq_constraint_start_cost;constraint_start_cost];
          
        if(t-ThTime_off_min(i)-ThCold_time_start(i))<=0&& max(0,-ThTime_on_off_init(i))<abs(t-ThTime_off_min(i)-ThCold_time_start(i)-1)+1
            bineq_constraint_start_cost=[bineq_constraint_start_cost;ThCold_cost_start(i)];
        else
            bineq_constraint_start_cost=[bineq_constraint_start_cost;0];
        end
    end
end
%% unit generation limit
Aineq_constraint_generation_limit=[];
bineq_constraint_generation_limit=sparse(2*N*T,1);

for i=1:N
    for t=1:T
        constraint_generation_limit=sparse(1,varaible_num);
        
        uit=sparse(1,N*T);
        pit=sparse(1,N*T);
        
        uit(1,(i-1)*T+t)=ThPimin(i);
        pit(1,(i-1)*T+t)=-1;
        
        constraint_generation_limit(1,1:N*T)=uit;
        constraint_generation_limit(1,3*N*T+1:4*N*T)=pit;
        
        Aineq_constraint_generation_limit=[Aineq_constraint_generation_limit;constraint_generation_limit];
    end
end


for i=1:N
    for t=1:T
        constraint_generation_limit=sparse(1,varaible_num);
        
        uit=sparse(1,N*T);
        pit=sparse(1,N*T);
        
        uit(1,(i-1)*T+t)=-ThPimax(i);
        pit(1,(i-1)*T+t)=1;
        
        constraint_generation_limit(1,1:N*T)=uit;
        constraint_generation_limit(1,3*N*T+1:4*N*T)=pit;
        
        Aineq_constraint_generation_limit=[Aineq_constraint_generation_limit;constraint_generation_limit];
    end
end
%% spinning reserve requirment
Aineq_constraints_spinning_reserve_requirment=[];
bineq_constraints_spinning_reserve_requirment=-Dt-Spin;
uit=[];
 constraints_spinning_reserve_requirment=sparse(1,varaible_num);
for i=1:N
    uit=[uit,-1*eye(T)*ThPimax(i)];
end
Aineq_constraints_spinning_reserve_requirment=[uit,sparse(T,varaible_num-N*T)];
%% ramp up constraints
Aineq_constraints_ramp_up=[];
bineq_constraints_ramp_up=[];
%t=1
for i=1:N
        constraints_ramp_up=sparse(1,varaible_num);
        
        uit=sparse(1,N*T);
        sit=sparse(1,N*T);
        pit=sparse(1,N*T);
        
        uit(1,(i-1)*T+1)=-1*(Piup(i)+ThPimin(i));
        sit(1,(i-1)*T+1)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));
         pit(1,(i-1)*T+1)=1;
        
        
         constraints_ramp_up(1,1:N*T)=uit;
         constraints_ramp_up(1,N*T+1:2*N*T)=sit;
         constraints_ramp_up(1,3*N*T+1:4*N*T)=pit;
         
         Aineq_constraints_ramp_up=[Aineq_constraints_ramp_up;constraints_ramp_up];
         bineq_constraints_ramp_up=[bineq_constraints_ramp_up;Thon_off_init(i)-Ui0(i)*ThPimin(i)];
         
end

%t>1
for i=1:N
    for t=2:T
        constraints_ramp_up=sparse(1,varaible_num);
        
        uit=sparse(1,N*T);
        sit=sparse(1,N*T);
        pit=sparse(1,N*T);
        
        uit(1,(i-1)*T+t)=-1*(Piup(i)+ThPimin(i));
        uit(1,(i-1)*T+t-1)=ThPimin(i);
        sit(1,(i-1)*T+t)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));
         pit(1,(i-1)*T+t)=1;
          pit(1,(i-1)*T+t-1)=-1;
        
        
         constraints_ramp_up(1,1:N*T)=uit;
         constraints_ramp_up(1,N*T+1:2*N*T)=sit;
         constraints_ramp_up(1,3*N*T+1:4*N*T)=pit;
         
         Aineq_constraints_ramp_up=[Aineq_constraints_ramp_up;constraints_ramp_up];
    end
end
bineq_constraints_ramp_up=[bineq_constraints_ramp_up;sparse(N*(T-1),1)];

%% ramp down constraints
Aineq_constraints_ramp_down=[];
bineq_constraints_ramp_down=[];
%t=1
for i=1:N
     constraints_ramp_down=sparse(1,varaible_num);
     
     uit=sparse(1,N*T);
     dit=sparse(1,N*T);
     pit=sparse(1,N*T);
     
     uit(1,(i-1)*T+1)=ThPimin(i);
     dit(1,(i-1)*T+1)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));
     pit(1,(i-1)*T+1)=-1;
     
     constraints_ramp_down(1,1:N*T)=uit;
     constraints_ramp_down(1,2*N*T+1:3*N*T)=dit;
     constraints_ramp_down(1,3*N*T+1:4*N*T)=pit;
     
     
     Aineq_constraints_ramp_down=[Aineq_constraints_ramp_down;constraints_ramp_down];
end
bineq_constraints_ramp_down=Ui0.*(Pidown+ThPimin)-p0;
%t>2
for i=1:N
    for t=2:T
        constraints_ramp_down=sparse(1,varaible_num);
        
        uit=sparse(1,N*T);
        dit=sparse(1,N*T);
        pit=sparse(1,N*T);
        
        uit(1,(i-1)*T+t)=ThPimin(i);
        uit(1,(i-1)*T+t-1)=-1*(Pidown(i)+ThPimin(i));
        dit(1,(i-1)*T+t)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));
        pit(1,(i-1)*T+t)=-1;
        pit(1,(i-1)*T+t-1)=1;
        
            constraints_ramp_down(1,1:N*T)=uit;
     constraints_ramp_down(1,2*N*T+1:3*N*T)=dit;
     constraints_ramp_down(1,3*N*T+1:4*N*T)=pit;
        
        
        Aineq_constraints_ramp_down=[Aineq_constraints_ramp_down;constraints_ramp_down];
    end
end
bineq_constraints_ramp_down=[bineq_constraints_ramp_down;sparse(N*(T-1),1)];
%% minim up/down time constraint
Aineq_constraint_min_upordown_time=[];
bineq_constraint_min_upordown_time=[];
for i=1:N
    for t=Ui(i) + 1:T
       
            constraint_min_upordown_time=sparse(1,varaible_num);
            
            uit=sparse(1,N*T);
            sit=sparse(1,N*T);
            
            uit(1,(i-1)*T+t)=-1;
            sit(1,(i-1)*T+max(0,t-ThTime_on_min(i))+1:(i-1)*T+t)=1;
            
            constraint_min_upordown_time(1,1:N*T)=uit;
            constraint_min_upordown_time(1,N*T+1:2*N*T)=sit;
            
            Aineq_constraint_min_upordown_time=[Aineq_constraint_min_upordown_time;constraint_min_upordown_time];
            bineq_constraint_min_upordown_time=[bineq_constraint_min_upordown_time;0];
      
    end
end

for i=1:N
    for t=Li(i)+1:T
       
            constraint_min_upordown_time=sparse(1,varaible_num);
            
            uit=sparse(1,N*T);
            dit=sparse(1,N*T);
            
            uit(1,(i-1)*T+t)=1;
            dit(1,(i-1)*T+max(0,t-ThTime_off_min(i))+1:(i-1)*T+t)=1;
            
            constraint_min_upordown_time(1,1:N*T)=uit;
            constraint_min_upordown_time(1,2*N*T+1:3*N*T)=dit;
            
            Aineq_constraint_min_upordown_time=[Aineq_constraint_min_upordown_time;constraint_min_upordown_time];
            bineq_constraint_min_upordown_time=[bineq_constraint_min_upordown_time;1];
    end
end
%% carbon emission constraints


%% objective function
%������
Pwan_coefficient=sparse(N*T,varaible_num);
Pwan_coefficient(1:N*T,3*N*T+1:4*N*T)=2*diag(reshape(repmat(Gama,1,T)',N*T,1));
H=sparse(varaible_num,varaible_num);
H(3*N*T+1:4*N*T,:)=Pwan_coefficient;

%һ����
f=sparse(varaible_num,1);
f(1:N*T,1)=reshape(repmat(Alpha,1,T)',N*T,1);
f(3*N*T+1:4*N*T,1)=reshape(repmat(Beta,1,T)',N*T,1);
f(4*N*T+1:5*N*T,1)=ones(N*T,1);

%����ʽԼ��
Aineq=[Aineq_constraints_ramp_up;Aineq_constraints_ramp_down;Aineq_constraint_generation_limit;Aineq_constraint_min_upordown_time;Aineq_constraint_start_cost;Aineq_constraints_spinning_reserve_requirment];%];%];%];
bineq=[bineq_constraints_ramp_up;bineq_constraints_ramp_down;bineq_constraint_generation_limit;bineq_constraint_min_upordown_time;bineq_constraint_start_cost;bineq_constraints_spinning_reserve_requirment];%];%];%];
%��ʽԼ��
Aeq=[Aeq_constraint_power_balance;Aeq_constraint_state;Aeq_constraints_initial_statues];%];
beq=[beq_constraint_power_balance;beq_constraint_state;beq_constraints_initial_statues];%];
%�������½�
lb=sparse(varaible_num,1);
ub=inf*ones(varaible_num,1);
ub(1:N*T,1)=ones(N*T,1);
ub(N*T+1:2*N*T,1)=ones(N*T,1);
ub(2*N*T+1:3*N*T,1)=ones(N*T,1);
 ub(3*N*T+1:4*N*T,1)=reshape(repmat(ThPimax,1,T)',N*T,1);

%��������
ctype='';
ctype(1,1:varaible_num)='C';
ctype(1,1:N*T)='B';
ctype(1,N*T+1:2*N*T)='B';
ctype(1,2*N*T+1:3*N*T)='B';

model.H=H;
model.f=f;
model.Aineq=Aineq;
model.bineq=bineq;
model.Aeq=Aeq;
model.beq=beq;
model.lb=lb;
model.ub=ub;
model.ctype=ctype;

%% ����cplex���

%���ýӿ�cplexmiqp
options = cplexoptimset;
options.Display='on';
% options.mip.tolerances.mipgap=5e-3;
[x,fval,exitflag,output]=cplexmiqp(H,f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
fprintf('\n solution status=%s\n',output.cplexstatusstring);
fprintf('Solution value = %f \n', fval);

%������
% cplex_milp = Cplex('Milp for HTC');
% cplex_milp.Model.sense = 'minimize';
% cplex_milp.Model.obj = f;
% cplex_milp.Model.lb = lb;
% cplex_milp.Model.ub = ub;
% cplex_milp.Model.A = [Aineq;Aeq];
% cplex_milp.Model.lhs = [-Inf.*ones(size(bineq,1),1);beq];
% cplex_milp.Model.rhs = [bineq;beq];
% cplex_milp.Model.ctype = ctype;
% cplex_milp.Model.Q = H;
% cplex_milp.Param.mip.tolerances.mipgap.Cur = 5e-3;
% % cplex_milp.Param.mip.pool.relgap.Cur = 0.9;
% % cplex_milp.Param.timelimit.Cur=2;
% % cplex_milp.Param.mip.limits.solutions.Cur=1;
%  cplex_milp.Param.mip.tolerances.uppercutoff.Cur=1.584e7;
% cplex_milp.solve();
% cplex_milp.writeSolution('solution.sol');
% numsoln = size (cplex_milp.Solution.pool.solution, 1);
% fprintf ('The solution pool contains %d solutions\n', numsoln);
% F = cplex_milp.Solution.objval;
% disp('F=');
% disp(F);
end