function model=model_produce(dataUC,CET)
% %% init
%     dataUC=readdataUC('UC_AF/5_std.mod');
Alpha = dataUC.alpha;                           %火电机组发电函数系数Alpha--N*1矩阵
Beta = dataUC.beta;                             %火电机组发电函数系数Beta--N*1矩阵
Gama = dataUC.gamma;                            %火电机组发电函数系数Gama--N*1矩阵
ThPimin = dataUC.p_low;                         %火电机组发电功率下界--N*1矩阵
ThPimax = dataUC.p_up;                          %火电机组发电功率上界--N*1矩阵
Piup = dataUC.p_rampup;                         %火电机组上坡功率上界--N*1矩阵
Pidown = dataUC.p_rampdown;                     %火电机组下坡功率上界--N*1矩阵
Dt = dataUC.PD;                                 %负载需求--T*1矩阵
N = dataUC.N;                                   %火电机组数量--1*1矩阵
T = dataUC.T;                                   %时间段数--1*1矩阵
Spin = dataUC.spin;                             %旋转热备用--T*1矩阵
ThTime_on_off_init = dataUC.time_on_off_ini;    %火电机组在初始状态前已经开机/停机的时间--N*1矩阵
Thon_off_init = dataUC.p_initial;               %火电机组机组初始功率--N*1矩阵
ThTime_on_min = dataUC.time_min_on;             %火电机组最小开机时间--N*1矩阵
ThTime_off_min = dataUC.time_min_off;           %火电机组最小停机时间--N*1矩阵
ThCold_cost_start = dataUC.Cold_cost;           %火电机组冷启动费用--N*1矩阵
ThHot_cost_start = dataUC.Hot_cost;             %火电机组热启动费用--N*1矩阵
ThCold_time_start = dataUC.Cold_hour;           %火电机组冷启动时间--N*1矩阵
Pistartup = dataUC.p_startup;                   %火电机组开机功率--N*1矩阵
Pishutdown = dataUC.p_shutdown;                 %火电机组关机功率--N*1矩阵
Ui0 = full(spones(Thon_off_init));                                   %火电机组机组初始开/停机状态
ThPit_wani0 = (Thon_off_init - Ui0 .* ThPimin) ./ (ThPimax - ThPimin);      %火电Pit弯弯初始值
Ui = max(0,min(ones(N,1) * T,Ui0 .* (ThTime_on_min - ThTime_on_off_init)));                    %--N*1矩阵
Li = max(0,min(ones(N,1) * T,(ones(N,1) - Ui0) .*  (ThTime_off_min + ThTime_on_off_init)));     %--N*1矩阵

% ThTime_on_min=2*ones(length(ThTime_on_min),1);
% ThHot_cost_start=8*ones(length(ThCold_cost_start),1);
% Pishutdown=1.01*Piup+ThPimin;
% Pistartup=1.01*Pidown+ThPimin;
% location_Pishutdown_besmall=find((Pishutdown-Pistartup)<0);
% location_Pishutdown_begreater=find((Pishutdown-Pistartup)>0);
% if(isempty(location_Pishutdown_besmall))
%     disp('Pishutdown>Pistartup');
% end
% if(isempty(location_Pishutdown_begreater))
%      disp('Pishutdown<Pistartup');
% end

%project
Alpha_wan = Alpha + Beta .* ThPimin + Gama .* ThPimin .* ThPimin;           %Alpha弯弯--N*1矩阵
Beta_wan = (ThPimax - ThPimin) .* (Beta + 2 * Gama .* ThPimin);             %Beta弯弯--N*1矩阵
Gama_wan = Gama .* (ThPimax - ThPimin) .* (ThPimax - ThPimin);              %Gama弯弯--N*1矩阵

Piup_wan = Piup ./ (ThPimax - ThPimin);                                     %Piup弯弯--N*1矩阵
Pidown_wan = Pidown ./ (ThPimax - ThPimin);                                 %Pidown弯弯--N*1矩阵
Pistartup_wan = (Pistartup - ThPimin) ./ (ThPimax - ThPimin);               %Pistartup弯弯--N*1矩阵
Pishutdown_wan = (Pishutdown - ThPimin) ./ (ThPimax - ThPimin);             %Pishutdown弯弯--N*1矩阵

% ThTime_on_min=1*ones(size(ThTime_on_min,1),1);


%variables
variable={'uit/oit','B',N*T;'Pit/PitWan','C',N*T;'sit','B',N*T;'Sit/SitWan','C',N*T;'dit','B',N*T;'T3','B',N*T};


%% startup cost constraint
Aineq_constraint_start_cost=[];
bineq_constraint_start_cost=[];

Aineq_constraint_start_cost_project=[];
bineq_constraint_start_cost_project=[];
%Sit>=C_hot,i*sit
Aineq_constraint_start_cost=[sparse(diag(reshape(repmat(ThHot_cost_start,1,T)',N*T,1))),sparse(-1*eye(N*T)),sparse(N*T,N*T)];
bineq_constraint_start_cost=sparse(N*T,1);
%Sit>=Ccoldi*[sit-sum（dit）-f_inint]
dit_orgin=[];
dit_project=[];
for i=1:N
    for t=1:T
        dit_part_orgin=sparse(1,N*T);
        dit_part_project=sparse(1,N*T);
        dit_part_orgin(1,(i-1)*T+max(t-ThTime_off_min(i)-ThCold_time_start(i),1):(i-1)*T+t-1)=-1*ThCold_cost_start(i);
        dit_part_project(1,(i-1)*T+max(t-ThTime_off_min(i)-ThCold_time_start(i),1):(i-1)*T+t-1)=-1*(ThCold_cost_start(i)-ThHot_cost_start(i));
        dit_orgin=[dit_orgin;dit_part_orgin];
        dit_project=[dit_project;dit_part_project];
        if(t-ThTime_off_min(i)-ThCold_time_start(i))<=0&& max(0,-ThTime_on_off_init(i))<abs(t-ThTime_off_min(i)-ThCold_time_start(i)-1)+1
            bineq_constraint_start_cost=[bineq_constraint_start_cost;ThCold_cost_start(i)];
            bineq_constraint_start_cost_project=[bineq_constraint_start_cost_project;ThCold_cost_start(i)-ThHot_cost_start(i)];
        else
            bineq_constraint_start_cost=[bineq_constraint_start_cost;0];
            bineq_constraint_start_cost_project=[bineq_constraint_start_cost_project;0];
        end
    end
end

Aineq_constraint_start_cost=[Aineq_constraint_start_cost;sparse(diag(reshape(repmat(ThCold_cost_start,1,T)',N*T,1))),sparse(-1*eye(N*T)),dit_orgin];%sit,Sit,dit
Aineq_constraint_start_cost_project=[sparse(diag(reshape(repmat(ThCold_cost_start-ThHot_cost_start,1,T)',N*T,1))),sparse(-1*eye(N*T)),dit_project];

Aineq_constraint_start_cost_2P_Co=[sparse(size(Aineq_constraint_start_cost,1),2*N*T),Aineq_constraint_start_cost];
bineq_constraint_start_cost_2P_Co=bineq_constraint_start_cost;

Aineq_constraint_start_cost_2P_Ti=[sparse(size(Aineq_constraint_start_cost,1),2*N*T),Aineq_constraint_start_cost];
bineq_constraint_start_cost_2P_Ti=bineq_constraint_start_cost;

Aineq_constraint_start_cost_3P_Ti=[sparse(size(Aineq_constraint_start_cost,1),2*N*T),Aineq_constraint_start_cost];
bineq_constraint_start_cost_3P_Ti=bineq_constraint_start_cost;

Aineq_constraint_start_cost_3P_Ti_ST=[sparse(size(Aineq_constraint_start_cost_project,1),2*N*T),Aineq_constraint_start_cost_project];
bineq_constraint_start_cost_3P_Ti_ST=bineq_constraint_start_cost_project;

Aineq_constraint_start_cost_3P_HD=[sparse(size(Aineq_constraint_start_cost_project,1),2*N*T),Aineq_constraint_start_cost_project,sparse(size(Aineq_constraint_start_cost_project,1),N*T)];
bineq_constraint_start_cost_3P_HD=bineq_constraint_start_cost_project;

Aineq_constraint_start_cost_3P_HD_Pr=[sparse(size(Aineq_constraint_start_cost_project,1),2*N*T),Aineq_constraint_start_cost_project];
bineq_constraint_start_cost_3P_HD_Pr=bineq_constraint_start_cost_project;
%% state variables and logical constraints
Aeq_constraint_state=[];
beq_constraint_state=[];

Aineq_constraint_state_3P_HD=[];
bineq_constraint_state_3P_HD=[];
%t=1时
for i=1:N
    
    uit=sparse(1,N*T);
    sit=sparse(1,N*T);
    dit=sparse(1,N*T);
    
    uit(1,(i-1)*T+1)=1;
    sit(1,(i-1)*T+1)=-1;
    dit(1,(i-1)*T+1)=1;
    
    Aeq_constraint_state=[Aeq_constraint_state;uit,sparse(1,N*T),sit,sparse(1,N*T),dit];
end
beq_constraint_state=Ui0;

Aeq_constraint_state_ST=Aeq_constraint_state;
beq_constraint_state_ST=beq_constraint_state;



%t>1

for i=1:N
    uit=[sparse(diag(-1*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(ones(1,T-1)))];
    sit=[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];
    dit=[sparse(T-1,1),sparse(diag(ones(1,T-1)))];
    
    constraint_state=sparse(T-1,5*N*T);
    constraint_state(1:T-1,(i-1)*T+1:i*T)=uit;
    constraint_state(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=sit;
    constraint_state(1:T-1,(i-1)*T+4*N*T+1:(i-1)*T+4*N*T+T)=dit;
    
    Aeq_constraint_state=[Aeq_constraint_state;constraint_state];
    %3P_Ti_ST
    %ui1-oi2-di2=0
    %oit-1 - oit + sit-1 -dit=0;
    constraint_state_ST=sparse(T-1,5*N*T);
    constraint_state_ST(1:T-1,(i-1)*T+1:(i-1)*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%oit
    constraint_state_ST(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=[sparse(1,T);sparse(T-2,1),sparse(diag(ones(1,T-2))),sparse(T-2,1)];%because sit=0 when t=2
    constraint_state_ST(1:T-1,(i-1)*T+4*N*T+1:(i-1)*T+4*N*T+T)=[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%dit
    
    Aeq_constraint_state_ST=[Aeq_constraint_state_ST;constraint_state_ST];
    beq_constraint_state_ST=[beq_constraint_state_ST;sparse(T-1,1)];
    
    %3P_HD
    %T3>=sit+ dit+1 -uit
    constraint_state_3P_HD=sparse(T-1,6*N*T);
    constraint_state_3P_HD(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit
    constraint_state_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(ones(1,T-1)));%sit
    constraint_state_3P_HD(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(ones(1,T-1)));%dit
    constraint_state_3P_HD(1:T-1,(i-1)*T+5*N*T+1:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*ones(1,T-1)));%T3
    Aineq_constraint_state_3P_HD=[Aineq_constraint_state_3P_HD;constraint_state_3P_HD];
    bineq_constraint_state_3P_HD=[bineq_constraint_state_3P_HD;sparse(T-1,1)];
    
    %T3<=sit
    constraint_state_3P_HD=sparse(T,6*N*T);
    constraint_state_3P_HD(1:T,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=sparse(diag(-1*ones(1,T)));%sit
    constraint_state_3P_HD(1:T,(i-1)*T+5*N*T+1:(i-1)*T+5*N*T+T)=sparse(diag(ones(1,T)));%T3
    Aineq_constraint_state_3P_HD=[Aineq_constraint_state_3P_HD;constraint_state_3P_HD];
    bineq_constraint_state_3P_HD=[bineq_constraint_state_3P_HD;sparse(T,1)];
    
    %T3<=dit+1
    constraint_state_3P_HD=sparse(T-1,6*N*T);
    constraint_state_3P_HD(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*ones(1,T-1)));%dit
    constraint_state_3P_HD(1:T-1,(i-1)*T+5*N*T+1:(i-1)*T+5*N*T+T-1)=sparse(diag(ones(1,T-1)));%T3
    Aineq_constraint_state_3P_HD=[Aineq_constraint_state_3P_HD;constraint_state_3P_HD];
    bineq_constraint_state_3P_HD=[bineq_constraint_state_3P_HD;sparse(T-1,1)];
end

beq_constraint_state=[beq_constraint_state;sparse(N*(T-1),1)];

Aeq_constraint_state_2P_Co=Aeq_constraint_state;
beq_constraint_state_2P_Co=beq_constraint_state;

Aeq_constraint_state_2P_Ti=Aeq_constraint_state;
beq_constraint_state_2P_Ti=beq_constraint_state;

Aeq_constraint_state_3P_Ti=Aeq_constraint_state;
beq_constraint_state_3P_Ti=beq_constraint_state;

Aeq_constraint_state_3P_Ti_ST=Aeq_constraint_state_ST;
beq_constraint_state_3P_Ti_ST=beq_constraint_state_ST;

Aeq_constraint_state_3P_HD=[Aeq_constraint_state,sparse(size(Aeq_constraint_state,1),N*T)];
beq_constraint_state_3P_HD=beq_constraint_state;

Aeq_constraint_state_3P_HD_Pr=Aeq_constraint_state;
beq_constraint_state_3P_HD_Pr=beq_constraint_state;

%% initial status of units
Aeq_constraints_initial_statues=[];
beq_constraints_initial_statues=[];
Aeq_constraints_initial_statues_3P_HD=[];
beq_constraints_initial_statues_3P_HD=[];
Aeq_constraints_initial_statues_3P_Ti_ST=[];
beq_constraints_initial_statues_3P_Ti_ST=[];
Aeq_constraints_initial_statues_3P_HD_Pr=[];
beq_constraints_initial_statues_3P_HD_Pr=[];

for i=1:N
    if(Ui(i) + Li(i) >= 1)
        for t=1:Ui(i)+Li(i)
            
            if(t==1)
                uit=sparse(1,N*T);
                uit(1,(i-1)*T+1)=1;
                
                constraints_initial_statues=sparse(1,5*N*T);
                constraints_initial_statues(1,1:N*T)=uit;
                
                Aeq_constraints_initial_statues_3P_Ti_ST=[Aeq_constraints_initial_statues_3P_Ti_ST;constraints_initial_statues];
                beq_constraints_initial_statues_3P_Ti_ST=[beq_constraints_initial_statues_3P_Ti_ST;Ui0(i)];
            else
                oit=sparse(1,N*T);
                sit=sparse(1,N*T);
                oit(1,(i-1)*T+t)=1;
                sit(1,(i-1)*T+t)=1;
                
                constraints_initial_statues=sparse(1,5*N*T);
                constraints_initial_statues(1,1:N*T)=oit; %oit
                constraints_initial_statues(1,2*N*T+1:3*N*T)=sit; %sit
                
                Aeq_constraints_initial_statues_3P_Ti_ST=[Aeq_constraints_initial_statues_3P_Ti_ST;constraints_initial_statues];
                beq_constraints_initial_statues_3P_Ti_ST=[beq_constraints_initial_statues_3P_Ti_ST;Ui0(i)];
            end
            
            uit=sparse(1,N*T);
            uit(1,(i-1)*T+t)=1;
            
            constraints_initial_statues=sparse(1,5*N*T);
            constraints_initial_statues(1,1:N*T)=uit;
            
            Aeq_constraints_initial_statues=[Aeq_constraints_initial_statues;constraints_initial_statues];
            beq_constraints_initial_statues=[beq_constraints_initial_statues;Ui0(i)];
            
            
        end
    end
end
Aeq_constraints_initial_statues_2P_Co=Aeq_constraints_initial_statues;
beq_constraints_initial_statues_2P_Co=beq_constraints_initial_statues;

Aeq_constraints_initial_statues_2P_Ti=Aeq_constraints_initial_statues;
beq_constraints_initial_statues_2P_Ti=beq_constraints_initial_statues;

Aeq_constraints_initial_statues_3P_Ti=Aeq_constraints_initial_statues;
beq_constraints_initial_statues_3P_Ti=beq_constraints_initial_statues;

if(size(Aeq_constraints_initial_statues,1)~=0)
    Aeq_constraints_initial_statues_3P_HD=[Aeq_constraints_initial_statues,sparse(size(Aeq_constraints_initial_statues,1),N*T)];
    beq_constraints_initial_statues_3P_HD=beq_constraints_initial_statues;
    
    Aeq_constraints_initial_statues_3P_HD_Pr=Aeq_constraints_initial_statues;
    beq_constraints_initial_statues_3P_HD_Pr=beq_constraints_initial_statues;
end


%% minimum up/dowm time constraints
Aineq_constraint_min_upordown_time=[];
bineq_constraint_min_upordown_time=[];
Aineq_constraint_min_upordown_time_substitude_o=[];
bineq_constraint_min_upordown_time_substitude_o=[];


for i=1:N
    uit=sparse(1,N*T);
    uit(1,(i-1)*T+Ui(i) + 1:(i-1)*T+T)=-1;
    uit=sparse(diag(uit));
    [nonzerosRow,nonzeroscol]=find(uit~=0);
    uit=uit(nonzerosRow,:);
    cons_row_num=size(uit,1);
    sit=[];
    
    for t=Ui(i) + 1:T
        sit_row=sparse(1,N*T);
        sit_row(1,(i-1)*T+max(0,t-ThTime_on_min(i))+1:(i-1)*T+t)=1;
        sit=[sit;sit_row];
    end
    Aineq_constraint_min_upordown_time=[Aineq_constraint_min_upordown_time;uit,sparse(cons_row_num,N*T),sit,sparse(cons_row_num,N*T),sparse(cons_row_num,N*T)];
    bineq_constraint_min_upordown_time=[bineq_constraint_min_upordown_time;sparse(cons_row_num,1)];
    
    uit=sparse(1,N*T);
    uit(1,(i-1)*T+Li(i)+1:(i-1)*T+T)=1;
    uit=sparse(diag(uit));
    [nonzerosRow,nonzeroscol]=find(uit~=0);
    uit=uit(nonzerosRow,:);
    cons_row_num=size(uit,1);
    dit=[];
    for t=Li(i)+1:T
        dit_row=sparse(1,N*T);
        dit_row(1,(i-1)*T+max(0,t-ThTime_off_min(i))+1:(i-1)*T+t)=1;
        dit=[dit;dit_row];
    end
    
    Aineq_constraint_min_upordown_time=[Aineq_constraint_min_upordown_time;uit,sparse(cons_row_num,N*T),sparse(cons_row_num,N*T),sparse(cons_row_num,N*T),dit];
    bineq_constraint_min_upordown_time=[bineq_constraint_min_upordown_time;sparse(ones(cons_row_num,1))];
    
    
    %substitude uit with oit
    sit=[];
    if(Ui(i) ==0)
        
        if(max(0,1-ThTime_on_min(i))<=0)  %t=1
            uit=sparse(1,N*T);
            sit=sparse(1,N*T);
            
            uit(1,(i-1)*T+1)=-1;
            sit(1,(i-1)*T+1)=1;
            Aineq_constraint_min_upordown_time_substitude_o=[Aineq_constraint_min_upordown_time_substitude_o;uit,sparse(1,N*T),sit,sparse(1,N*T),sparse(1,N*T)];
            bineq_constraint_min_upordown_time_substitude_o=[bineq_constraint_min_upordown_time_substitude_o;0];
        end
        
        %t>1
        oit=sparse(1,N*T);
        oit(1,(i-1)*T+2:(i-1)*T+T)=-1;
        oit=sparse(diag(oit));
        [nonzerosRow,nonzeroscol]=find(oit~=0);
        oit=oit(nonzerosRow,:);
        cons_row_num=size(oit,1);
        sit=[];
        for t= 2:T
            sit_row=sparse(1,N*T);
            sit_row(1,(i-1)*T+max(0,t-ThTime_on_min(i))+1:(i-1)*T+t-1)=1;
            sit=[sit;sit_row];
        end
        Aineq_constraint_min_upordown_time_substitude_o=[Aineq_constraint_min_upordown_time_substitude_o;oit,sparse(cons_row_num,N*T),sit,sparse(cons_row_num,N*T),sparse(cons_row_num,N*T)];
        bineq_constraint_min_upordown_time_substitude_o=[bineq_constraint_min_upordown_time_substitude_o;sparse(cons_row_num,1)];
    else
        oit=sparse(1,N*T);
        oit(1,(i-1)*T+Ui(i) + 1:(i-1)*T+T)=-1;
        oit=sparse(diag(oit));
        [nonzerosRow,nonzeroscol]=find(oit~=0);
        oit=oit(nonzerosRow,:);
        cons_row_num=size(oit,1);
        sit=[];
        for t=Ui(i) + 1:T
            sit_row=sparse(1,N*T);
            sit_row(1,(i-1)*T+max(0,t-ThTime_on_min(i))+1:(i-1)*T+t-1)=1;
            sit=[sit;sit_row];
        end
        Aineq_constraint_min_upordown_time_substitude_o=[Aineq_constraint_min_upordown_time_substitude_o;oit,sparse(cons_row_num,N*T),sit,sparse(cons_row_num,N*T),sparse(cons_row_num,N*T)];
        bineq_constraint_min_upordown_time_substitude_o=[bineq_constraint_min_upordown_time_substitude_o;sparse(cons_row_num,1)];
    end
    
    
    
    if(Li(i)==0)
        if(max(0,1-ThTime_off_min(i))<=0)  %t=1
            uit=sparse(1,N*T);
            dit=sparse(1,N*T);
            
            uit(1,(i-1)*T+1)=1;
            dit(1,(i-1)*T+1)=1;
            Aineq_constraint_min_upordown_time_substitude_o=[Aineq_constraint_min_upordown_time_substitude_o;uit,sparse(1,N*T),sparse(1,N*T),sparse(1,N*T),dit];
            bineq_constraint_min_upordown_time_substitude_o=[bineq_constraint_min_upordown_time_substitude_o;1];
        end
        
        %t>1
        oit=sparse(1,N*T);
        oit(1,(i-1)*T+2:(i-1)*T+T)=1;
        oit=sparse(diag(oit));
        [nonzerosRow,nonzeroscol]=find(oit~=0);
        oit=oit(nonzerosRow,:);
        sit=oit;
        cons_row_num=size(oit,1);
        dit=[];
        for t=2:T
            dit_row=sparse(1,N*T);
            dit_row(1,(i-1)*T+max(0,t-ThTime_off_min(i))+1:(i-1)*T+t)=1;
            dit=[dit;dit_row];
        end
        Aineq_constraint_min_upordown_time_substitude_o=[Aineq_constraint_min_upordown_time_substitude_o;oit,sparse(cons_row_num,N*T),sit,sparse(cons_row_num,N*T),dit];
        bineq_constraint_min_upordown_time_substitude_o=[bineq_constraint_min_upordown_time_substitude_o;sparse(ones(cons_row_num,1))];
    else
        oit=sparse(1,N*T);
        oit(1,(i-1)*T+Li(i)+1:(i-1)*T+T)=1;
        oit=sparse(diag(oit));
        [nonzerosRow,nonzeroscol]=find(oit~=0);
        oit=oit(nonzerosRow,:);
        sit=oit;
        cons_row_num=size(oit,1);
        dit=[];
        for t=Li(i)+1:T
            dit_row=sparse(1,N*T);
            dit_row(1,(i-1)*T+max(0,t-ThTime_off_min(i))+1:(i-1)*T+t)=1;
            dit=[dit;dit_row];
        end
        Aineq_constraint_min_upordown_time_substitude_o=[Aineq_constraint_min_upordown_time_substitude_o;oit,sparse(cons_row_num,N*T),sit,sparse(cons_row_num,N*T),dit];
        bineq_constraint_min_upordown_time_substitude_o=[bineq_constraint_min_upordown_time_substitude_o;sparse(ones(cons_row_num,1))];
    end
    
    
end

Aineq_constraint_min_upordown_time_2P_Co=Aineq_constraint_min_upordown_time;
bineq_constraint_min_upordown_time_2P_Co=bineq_constraint_min_upordown_time;

Aineq_constraint_min_upordown_time_2P_Ti=Aineq_constraint_min_upordown_time;
bineq_constraint_min_upordown_time_2P_Ti=bineq_constraint_min_upordown_time;

Aineq_constraint_min_upordown_time_3P_Ti=Aineq_constraint_min_upordown_time;
bineq_constraint_min_upordown_time_3P_Ti=bineq_constraint_min_upordown_time;

Aineq_constraint_min_upordown_time_3P_Ti_ST=Aineq_constraint_min_upordown_time_substitude_o;
bineq_constraint_min_upordown_time_3P_Ti_ST=bineq_constraint_min_upordown_time_substitude_o;

Aineq_constraint_min_upordown_time_3P_HD=[Aineq_constraint_min_upordown_time,sparse(size(Aineq_constraint_min_upordown_time,1),N*T)];
bineq_constraint_min_upordown_time_3P_HD=bineq_constraint_min_upordown_time;

Aineq_constraint_min_upordown_time_3P_HD_Pr=Aineq_constraint_min_upordown_time;
bineq_constraint_min_upordown_time_3P_HD_Pr=bineq_constraint_min_upordown_time;
%% unit generation limits and upper bound
%init
Aineq_constraint_generation_limit_3P_Ti_ST=[];
bineq_constraint_generation_limit_3P_Ti_ST=[];
Aineq_constraint_generation_limit_2P_Ti=[];
bineq_constraint_generation_limit_2P_Ti=[];
Aineq_constraint_generation_limit_3P_Ti=[];
bineq_constraint_generation_limit_3P_Ti=[];
Aineq_constraint_generation_limit_3P_HD=[];
bineq_constraint_generation_limit_3P_HD=[];
Aineq_constraint_generation_limit_3P_HD_Pr=[];
bineq_constraint_generation_limit_3P_HD_Pr=[];

%uit*Pmin<=pit<=uit*Pmax
Aineq_constraint_generation_limit=[sparse(diag(reshape(repmat(ThPimin,1,T)',N*T,1))),sparse(diag(-1*ones(1,N*T)));...%uit*Pmin<=pit
    sparse(diag(reshape(repmat(-1*ThPimax,1,T)',N*T,1))),sparse(diag(ones(1,N*T)))];%pit<=uit*Pmax
bineq_constraint_generation_limit=sparse(2*N*T,1);

%uit*Pmin<=pit
Aineq_constraint_generation_limit_3P_Ti=[sparse(diag(reshape(repmat(ThPimin,1,T)',N*T,1))),sparse(diag(-1*ones(1,N*T)))];
bineq_constraint_generation_limit_3P_Ti=sparse(N*T,1);

Aineq_constraint_generation_limit_2P_Ti=[sparse(diag(reshape(repmat(ThPimin,1,T)',N*T,1))),sparse(diag(-1*ones(1,N*T)))];
bineq_constraint_generation_limit_2P_Ti=sparse(N*T,1);

%uit*Pmin<=pit t=1
location=1:T:N*T;
Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;Aineq_constraint_generation_limit(location,:),sparse(size(Aineq_constraint_generation_limit(location,:),1),3*N*T)];
bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;sparse(N,1)];

%pit<=uit*Pmax t=1
location=N*T+1:T:2*N*T;
Aineq_constraint_generation_limit_3P_Ti=[Aineq_constraint_generation_limit_3P_Ti;Aineq_constraint_generation_limit(location,:)];
bineq_constraint_generation_limit_3P_Ti=[bineq_constraint_generation_limit_3P_Ti;sparse(N,1)];

Aineq_constraint_generation_limit_2P_Ti=[Aineq_constraint_generation_limit_2P_Ti;Aineq_constraint_generation_limit(location,:)];
bineq_constraint_generation_limit_2P_Ti=[bineq_constraint_generation_limit_2P_Ti;sparse(N,1)];

Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;Aineq_constraint_generation_limit(location,:),sparse(size(Aineq_constraint_generation_limit(location,:),1),3*N*T)];
bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;sparse(N,1)];

%pit<=uit*Pmax t=T
location=N*T+T:T:2*N*T;
Aineq_constraint_generation_limit_2P_Ti=[Aineq_constraint_generation_limit_2P_Ti,sparse(size(Aineq_constraint_generation_limit_2P_Ti,1),3*N*T);Aineq_constraint_generation_limit(location,:),sparse(size(Aineq_constraint_generation_limit(location,:),1),3*N*T)];
bineq_constraint_generation_limit_2P_Ti=[bineq_constraint_generation_limit_2P_Ti;sparse(N,1)];

Aineq_constraint_generation_limit_3P_Ti=[Aineq_constraint_generation_limit_3P_Ti,sparse(size(Aineq_constraint_generation_limit_3P_Ti,1),3*N*T);Aineq_constraint_generation_limit(location,:),sparse(size(Aineq_constraint_generation_limit(location,:),1),3*N*T)];
bineq_constraint_generation_limit_3P_Ti=[bineq_constraint_generation_limit_3P_Ti;sparse(N,1)];

%3P_Ti_ST   ui1*Pmin<=pi1 t=1   oit*Pmin+sit*Pmin<=pit t>1   pit<=oit*Pmax+sit*Pmax  t=T
for i=1:N
    
    %oit*Pmin+sit*Pmin<=pit t>1
    constraint_generation_limit_3P_Ti_ST=sparse(1,5*N*T);
    constraint_generation_limit_3P_Ti_ST(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(ThPimin(i)*ones(1,T-1)));%oit
    constraint_generation_limit_3P_Ti_ST(1:T-1,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=sparse(diag(-1*ones(1,T-1)));%pit
    constraint_generation_limit_3P_Ti_ST(1:T-1,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=sparse(diag(ThPimin(i)*ones(1,T-1)));%sit
    Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
    bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;sparse(T-1,1)];
    
    %pit<=oit*Pmax+sit*Pmax  t=T
    constraint_generation_limit_3P_Ti_ST=sparse(1,5*N*T);
    constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+T)=-1*ThPimax(i);%oit
    constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+N*T+T)=1;%pit
    constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+2*N*T+T)=-1*ThPimax(i);%sit
    Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
    bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;0];
end


for i=1:N
    if(ThTime_on_min(i)>1)%ThTime_on_min>=2
        %3P_Ti model    pit<=uit*Pmax-sit(Pmax-Pstart)- dit+1 *(Pmax-Pshut)
        constraint_generation_limit_3P_Ti=sparse(T-1,5*N*T);
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-1)));%uit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag((ThPimax(i)-Pistartup(i))*ones(1,T-1)));%sit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((ThPimax(i)-Pishutdown(i))*ones(1,T-1)));%dit+1
        
        Aineq_constraint_generation_limit_3P_Ti=[Aineq_constraint_generation_limit_3P_Ti;constraint_generation_limit_3P_Ti];
        bineq_constraint_generation_limit_3P_Ti=[bineq_constraint_generation_limit_3P_Ti;sparse(T-1,1)];
        
        %2P_Ti    pit<=uit*Pmax
        constraint_generation_limit_2P_Ti=sparse(T,5*N*T);
        constraint_generation_limit_2P_Ti(1:T,(i-1)*T+1:(i-1)*T+T)=sparse(diag(-1*ThPimax(i)*ones(1,T)));%uit
        constraint_generation_limit_2P_Ti(1:T,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=sparse(diag(ones(1,T)));%pit
        
        Aineq_constraint_generation_limit_2P_Ti=[Aineq_constraint_generation_limit_2P_Ti;constraint_generation_limit_2P_Ti];
        bineq_constraint_generation_limit_2P_Ti=[bineq_constraint_generation_limit_2P_Ti;sparse(T,1)];
        
        %3P_Ti_ST
        %pit<=uit*Pmax-sit*(Pmax-Pstart)- dit+1 *(Pmax-Pshut)   t=1
        constraint_generation_limit_3P_Ti_ST=sparse(1,5*N*T);
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+1)=-1*ThPimax(i);%uit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+N*T+1)=1;%pit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+2*N*T+1)=ThPimax(i)-Pistartup(i);%sit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+4*N*T+2)=ThPimax(i)-Pishutdown(i);%dit+1
        Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
        bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;0];
        %pit<=oit*Pmax+sit*Pstar- dit+1 *(Pmax-Pshut)        t>=2
        constraint_generation_limit_3P_Ti_ST=sparse(T-2,5*N*T);
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-2)));%oit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-2)));%pit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*Pistartup(i)*ones(1,T-2)));%sit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((ThPimax(i)-Pishutdown(i))*ones(1,T-2)));%dit+1
        Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
        bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;sparse(T-2,1)];
        
        %3P_HD
        %pitwan<=uit -sit*(1-Pwan,start)-dit+1 *(1-Pwan,shut)
        %t>=1
        constraint_generation_limit_3P_HD=sparse(T-1,6*N*T);
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pitwan
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)));%sit
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((1-Pishutdown_wan(i))*ones(1,T-1)));%dit+1
        Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
        bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;sparse(T-1,1)];
        
         %pit-1wan<= uit-1 +dit*(Pwan,shut-1)+dit+1*(Pwan,down+Pwan,shut-1)
         %t=1
         constraint_generation_limit_3P_HD=sparse(1,6*N*T);
         constraint_generation_limit_3P_HD(1,(i-1)*T+4*N*T+1)=1-Pishutdown_wan(i);%di1
         constraint_generation_limit_3P_HD(1,(i-1)*T+4*N*T+2)=1-Pishutdown_wan(i)-Pidown_wan(i);%di2
         Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
         bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;Ui0(i)-ThPit_wani0(i)];
         %t>1
         constraint_generation_limit_3P_HD=sparse(T-2,6*N*T);
         constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-1*ones(1,T-2)));%uit
         constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-2)=sparse(diag(ones(1,T-2)));%pitwan
         constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag((1-Pishutdown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((1-Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
         Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
         bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;sparse(T-2,1)];
         
         
         %pit+1wan<= uit+1 +sit*(Pwan,start + Pwan,up - 1)+ sit+1 *(Pwan,start-1)
         %t>=1
         constraint_generation_limit_3P_HD=sparse(T-1,6*N*T);
         constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(-1*ones(1,T-1)));%uit
         constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=sparse(diag(ones(1,T-1)));%pitwan
         constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)+Piup_wan(i)-1)*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)))];%sit
         Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
         bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;sparse(T-1,1)];
        
        %3P_HD_Pr
        %pitwan<=uit -sit*(1-Pwan,start)-dit+1 *(1-Pwan,shut)
        constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)));%sit
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((1-Pishutdown_wan(i))*ones(1,T-1)));%dit+1
        Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
        bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];
        
        %pit-1wan <= uit-1 + dit * (Pwan,shut - 1) +dit+1 * (Pwan,down + Pwan,shut -1)
        %t=1
        constraint_generation_limit_3P_HD_Pr=sparse(1,5*N*T);
        constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=(1-Pishutdown_wan(i));%dit
        constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=(1-Pishutdown_wan(i)-Pidown_wan(i));%dit+1
        Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
        bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;Ui0(i)-ThPit_wani0(i)];        
        
        %t>1
        constraint_generation_limit_3P_HD_Pr=sparse(T-2,5*N*T);
        constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-1*ones(1,T-2)));%uit
        constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-2)=sparse(diag(ones(1,T-2)));%pit
        constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag((1-Pishutdown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((1-Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
        Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
        bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-2,1)];        
        
        %pit+1wan <= uit+1 + sit *(Pwan,start +Pwan,up -1) +sit+1 *(Pwan,start - 1)
        %t>=1
        constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(-1*ones(1,T-1)));%uit
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=sparse(diag(ones(1,T-1)));%pit
        constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=[sparse(diag((1-Pistartup_wan(i)-Piup_wan(i))*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)))];%sit
        Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
        bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];          
        
        
    else%ThTime_on_min=1
        %3P_Ti
        %pit<=uit*Pmax-dit+1 *(Pmax-Pshut)-sit*(max(Pshut-Pstart,0))
        constraint_generation_limit_3P_Ti=sparse(T-1,5*N*T);
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-1)));%uit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max((Pishutdown(i)-Pistartup(i)),0)*ones(1,T-1)));%sit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((ThPimax(i)-Pishutdown(i))*ones(1,T-1)));%dit+1
        Aineq_constraint_generation_limit_3P_Ti=[Aineq_constraint_generation_limit_3P_Ti;constraint_generation_limit_3P_Ti];
        bineq_constraint_generation_limit_3P_Ti=[bineq_constraint_generation_limit_3P_Ti;sparse(T-1,1)];
        
        %pit<=uit*Pmax-sit*(Pmax-Pstart)-dit+1*(max(Pstart-Pshut,0))
        constraint_generation_limit_3P_Ti=sparse(T-1,5*N*T);
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-1)));%uit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag((ThPimax(i)-Pistartup(i))*ones(1,T-1)));%sit
        constraint_generation_limit_3P_Ti(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(max((Pistartup(i)-Pishutdown(i)),0)*ones(1,T-1)));%dit+1
        Aineq_constraint_generation_limit_3P_Ti=[Aineq_constraint_generation_limit_3P_Ti;constraint_generation_limit_3P_Ti];
        bineq_constraint_generation_limit_3P_Ti=[bineq_constraint_generation_limit_3P_Ti;sparse(T-1,1)];
        
        %2P_Ti
        %pit<=uit*Pmax-sit*(Pmax-Pstart)
        constraint_generation_limit_2P_Ti=sparse(T,5*N*T);
        constraint_generation_limit_2P_Ti(1:T,(i-1)*T+1:(i-1)*T+T)=sparse(diag(-1*ThPimax(i)*ones(1,T)));%uit
        constraint_generation_limit_2P_Ti(1:T,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=sparse(diag(ones(1,T)));%pit
        constraint_generation_limit_2P_Ti(1:T,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=sparse(diag((ThPimax(i)-Pistartup(i))*ones(1,T)));%sit
        
        Aineq_constraint_generation_limit_2P_Ti=[Aineq_constraint_generation_limit_2P_Ti;constraint_generation_limit_2P_Ti];
        bineq_constraint_generation_limit_2P_Ti=[bineq_constraint_generation_limit_2P_Ti;sparse(T,1)];
        
        %pit<=uit*Pmax-dit+1 (Pmax-Pshut)
        constraint_generation_limit_2P_Ti=sparse(T-1,5*N*T);
        constraint_generation_limit_2P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-1)));%uit
        constraint_generation_limit_2P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
        constraint_generation_limit_2P_Ti(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((ThPimax(i)-Pishutdown(i))*ones(1,T-1)));%dit+1
        
        Aineq_constraint_generation_limit_2P_Ti=[Aineq_constraint_generation_limit_2P_Ti;constraint_generation_limit_2P_Ti];
        bineq_constraint_generation_limit_2P_Ti=[bineq_constraint_generation_limit_2P_Ti;sparse(T-1,1)];
        
        %3P_Ti_ST
        %pit<=uit*Pmax- dit+1 *(Pmax-Pshut)-sit*(max(Pshut-Pstart,0))  t=1
        constraint_generation_limit_3P_Ti_ST=sparse(1,5*N*T);
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+1)=-1*ThPimax(i);%uit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+N*T+1)=1;%pit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+2*N*T+1)=max(ThPimax(i)-Pistartup(i),0);%sit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+4*N*T+2)=ThPimax(i)-Pishutdown(i);%dit+1
        Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
        bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;0];
        
        %pit<=uit*Pmax-sit*(Pmax-Pstart)- dit+1 *(max(Pstart-Pshut,0))  t=1
        constraint_generation_limit_3P_Ti_ST=sparse(1,5*N*T);
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+1)=-1*ThPimax(i);%uit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+N*T+1)=1;%pit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+2*N*T+1)=ThPimax(i)-Pistartup(i);%sit
        constraint_generation_limit_3P_Ti_ST(1,(i-1)*T+4*N*T+2)=max(Pistartup(i)-Pishutdown(i),0);%dit+1
        Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
        bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;0];
        
        %pit<=oit*Pmax+sit*(Pmax-max(Pshut-Pstar,0))-dit+1*(Pmax-Pshut)  t>=2
        constraint_generation_limit_3P_Ti_ST=sparse(T-2,5*N*T);
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-2)));%oit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-2)));%pit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(ThPimax(i)-max(Pishutdown(i)-Pistartup(i),0))*ones(1,T-2)));%sit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((ThPimax(i)-Pishutdown(i))*ones(1,T-2)));%dit+1
        Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
        bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;sparse(T-2,1)];
        
        %pit<=oit*Pmax+sit*Pstart-dit+1*max(Pstart-Pshut,0)  t>=2
        constraint_generation_limit_3P_Ti_ST=sparse(T-2,5*N*T);
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*ThPimax(i)*ones(1,T-2)));%oit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-2)));%pit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*Pistartup(i)*ones(1,T-2)));%sit
        constraint_generation_limit_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag(max(Pistartup(i)-Pishutdown(i),0)*ones(1,T-2)));%dit+1
        Aineq_constraint_generation_limit_3P_Ti_ST=[Aineq_constraint_generation_limit_3P_Ti_ST;constraint_generation_limit_3P_Ti_ST];
        bineq_constraint_generation_limit_3P_Ti_ST=[bineq_constraint_generation_limit_3P_Ti_ST;sparse(T-2,1)];
        
        %3P_HD
        %pitwan<=uit -sit*(1-Pwan,start)-dit+1 *(1-Pwan,shut)+T3{1-max(Pwan,start,Pwan,shut)}
        %t>=1
        constraint_generation_limit_3P_HD=sparse(T-1,6*N*T);
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pitwan
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)));%sit
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((1-Pishutdown_wan(i))*ones(1,T-1)));%dit+1
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+5*N*T+1:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(1-max(Pistartup_wan(i),Pishutdown_wan(i)))*ones(1,T-1)));%T3
        Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
        bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;sparse(T-1,1)];
        
        %pit-1wan<=T3*(1-Pwan,down-Pwan,shut)+ uit-1 +dit*(Pwan,shut-1)+dit+1*(Pwan,down+Pwan,shut-1)
        %t=1
        constraint_generation_limit_3P_HD=sparse(1,6*N*T);
        constraint_generation_limit_3P_HD(1,(i-1)*T+4*N*T+1)=1-Pishutdown_wan(i);%di1
        constraint_generation_limit_3P_HD(1,(i-1)*T+4*N*T+2)=1-Pishutdown_wan(i)-Pidown_wan(i);%di2
        constraint_generation_limit_3P_HD(1,(i-1)*T+5*N*T+1)=-1*(1-Pidown_wan(i)-Pishutdown_wan(i));%T3
        Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
        bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;Ui0(i)-ThPit_wani0(i)];
        %t>1
        constraint_generation_limit_3P_HD=sparse(T-2,6*N*T);
        constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-1*ones(1,T-2)));%uit
        constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-2)=sparse(diag(ones(1,T-2)));%pitwan
        constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag((1-Pishutdown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((1-Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
        constraint_generation_limit_3P_HD(1:T-2,(i-1)*T+5*N*T+2:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(1-Pidown_wan(i)-Pishutdown_wan(i))*ones(1,T-2)));%T3
        Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
        bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;sparse(T-2,1)];
%         
%         
        %pit+1wan<=T3*(1-Pwan,start - Pwan,up)+ uit+1 +sit*(Pwan,start + Pwan,up - 1)+ sit+1 *(Pwan,start-1)
        %t>=1
        constraint_generation_limit_3P_HD=sparse(T-1,6*N*T);
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(-1*ones(1,T-1)));%uit
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=sparse(diag(ones(1,T-1)));%pitwan
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)+Piup_wan(i)-1)*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)))];%sit
        constraint_generation_limit_3P_HD(1:T-1,(i-1)*T+5*N*T+1:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(1-Pistartup_wan(i)-Piup_wan(i))*ones(1,T-1)));%T3
        Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
        bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;sparse(T-1,1)];
 
        
            %3P_HD_Pr
            %pitwan<=uit -sit*(1-Pwan,start)-dit+1 *(max(PstartWan,PshutWan)-PshutWan)
            constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)));%sit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((max(Pistartup_wan(i),Pishutdown_wan(i))-Pishutdown_wan(i))*ones(1,T-1)));%dit+1
            Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
            bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];
            
            %pitwan<=uit +sit*(PstartWan-max(PshutWan,PstartWan))-dit+1 *(1-Pwan,shut)
            constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup_wan(i)-max(Pishutdown_wan(i),Pistartup_wan(i)))*ones(1,T-1)));%sit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((1-Pishutdown_wan(i))*ones(1,T-1)));%dit+1
            Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
            bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];
            
%             %pit-1wan <= uit-1 + dit*(Pwan,shut -1)
%             %t=1
%             constraint_generation_limit_3P_HD_Pr=sparse(1,5*N*T);
%             constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=1-Pishutdown_wan(i);%di1
%             Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
%             bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;Ui0(i)-ThPit_wani0(i)];
%             
%             %t>1
%             constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
%             constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*ones(1,T-1)));%uit-1
%             constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=sparse(diag(ones(1,T-1)));%pit-1
%             constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((1-Pishutdown_wan(i))*ones(1,T-1)));%dit
%             Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
%             bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];
            
            
            %pit-1wan <= uit-1 +sit*(1 - Pwan,down - Pwan,shut)+dit*(Pwan,shut -1)+ dit+1*(Pwan,down+Pwan,shut -1)
            %t=1
            constraint_generation_limit_3P_HD_Pr=sparse(1,5*N*T);
            constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(1-Pidown_wan(i)-Pishutdown_wan(i));%sit
            constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=(1-Pishutdown_wan(i));%dit
            constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=(1-Pishutdown_wan(i)-Pidown_wan(i));%dit+1
            Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
            bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;Ui0(i)-ThPit_wani0(i)];
            
            %t>1
            constraint_generation_limit_3P_HD_Pr=sparse(T-2,5*N*T);
            constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-1*ones(1,T-2)));%uit
            constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-2)=sparse(diag(ones(1,T-2)));%pit
            constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(1-Pidown_wan(i)-Pishutdown_wan(i))*ones(1,T-2)));%sit
            constraint_generation_limit_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag((1-Pishutdown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((1-Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
            Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
            bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-2,1)];            
            
%             %pit+1wan <= uit+1 + sit+1*(Pwan,start -1)
%             %t>=1
%             constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
%             constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(-1*ones(1,T-1)));%uit
%             constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=sparse(diag(ones(1,T-1)));%pit
%             constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)));%sit
%             Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
%             bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];            

            
            
            %pit+1wan <= uit+1 + dit+1*(1-Pwan,start -Pwan,up) + sit*(Pwan,start +Pwan,up -1) + sit+1*(Pwan,start -1)
            %t>=1
            constraint_generation_limit_3P_HD_Pr=sparse(T-1,5*N*T);
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(-1*ones(1,T-1)));%uit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=sparse(diag(ones(1,T-1)));%pit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T)=[sparse(diag((1-Pistartup_wan(i)-Piup_wan(i))*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag((1-Pistartup_wan(i))*ones(1,T-1)))];%sit
            constraint_generation_limit_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag((Pistartup_wan(i)+Piup_wan(i)-1)*ones(1,T-1)));%dit
            Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
            bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;sparse(T-1,1)];            
           
    end
   
    
        %pit-1wan <= uit-1 - dit*(1-Pwan,shut)  t=2
%     constraint_generation_limit_3P_HD=sparse(1,6*N*T);
%     constraint_generation_limit_3P_HD(1,(i-1)*T+1)=-1;%uit
%     constraint_generation_limit_3P_HD(1,(i-1)*T+N*T+1)=1;%pitwan
%     constraint_generation_limit_3P_HD(1,(i-1)*T+4*N*T+2)=1-Pishutdown_wan(i);%dit
%     Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
%     bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;0];
    
    %pit+1wan <= uit+1 - sit+1 *(1-Pwan,start)  t=T-1
%     constraint_generation_limit_3P_HD=sparse(1,6*N*T);
%     constraint_generation_limit_3P_HD(1,(i-1)*T+T)=-1;%uit
%     constraint_generation_limit_3P_HD(1,(i-1)*T+N*T+T)=1;%pitwan
%     constraint_generation_limit_3P_HD(1,(i-1)*T+2*N*T+T)=1-Pistartup_wan(i);%sit
%     Aineq_constraint_generation_limit_3P_HD=[Aineq_constraint_generation_limit_3P_HD;constraint_generation_limit_3P_HD];
%     bineq_constraint_generation_limit_3P_HD=[bineq_constraint_generation_limit_3P_HD;0];
    
    %3P_HD_Pr
    %pit-1wan <= uit-1 - dit*(1-Pwan,shut)  t=2
    constraint_generation_limit_3P_HD_Pr=sparse(1,5*N*T);
    constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+1)=-1;%uit
    constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pitwan
    constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=1-Pishutdown_wan(i);%dit
    Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
    bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;0];
    
    %pit+1wan <= uit+1 - sit+1 *(1-Pwan,start)  t=T-1
    constraint_generation_limit_3P_HD_Pr=sparse(1,5*N*T);
    constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+T)=-1;%uit
    constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+N*T+T)=1;%pitwan
    constraint_generation_limit_3P_HD_Pr(1,(i-1)*T+2*N*T+T)=1-Pistartup_wan(i);%sit
    Aineq_constraint_generation_limit_3P_HD_Pr=[Aineq_constraint_generation_limit_3P_HD_Pr;constraint_generation_limit_3P_HD_Pr];
    bineq_constraint_generation_limit_3P_HD_Pr=[bineq_constraint_generation_limit_3P_HD_Pr;0];
    
end

Aineq_constraint_generation_limit_2P_Co=[Aineq_constraint_generation_limit,sparse(size(Aineq_constraint_generation_limit,1),3*N*T)];
bineq_constraint_generation_limit_2P_Co=bineq_constraint_generation_limit;


%% Ramping constraints
%init
Aineq_constraints_ramp_up_2P_Co=[];
bineq_constraints_ramp_up_2P_Co=[];
Aineq_constraints_ramp_up_2P_Ti=[];
bineq_constraints_ramp_up_2P_Ti=[];
Aineq_constraints_ramp_up_3P_Ti=[];
bineq_constraints_ramp_up_3P_Ti=[];
Aineq_constraints_ramp_up_3P_HD=[];
bineq_constraints_ramp_up_3P_HD=[];
Aineq_constraints_ramp_up_3P_Ti_ST=[];
bineq_constraints_ramp_up_3P_Ti_ST=[];
Aineq_constraints_ramp_up_3P_HD_Pr=[];
bineq_constraints_ramp_up_3P_HD_Pr=[];

for i=1:N
    %t=1
    %2P_Ti
    %pi1-pi0<=ui1*(Pup+Pmin)-ui0*Pmin+si1*(Pstart-Pup-Pmin)
    constraints_ramp_up_2P_Ti=sparse(1,5*N*T);
    constraints_ramp_up_2P_Ti(1,(i-1)*T+1)=-1*(Piup(i)+ThPimin(i));%uit
    constraints_ramp_up_2P_Ti(1,(i-1)*T+N*T+1)=1;%pit
    constraints_ramp_up_2P_Ti(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));%sit
    Aineq_constraints_ramp_up_2P_Ti=[Aineq_constraints_ramp_up_2P_Ti;constraints_ramp_up_2P_Ti];
    bineq_constraints_ramp_up_2P_Ti=[bineq_constraints_ramp_up_2P_Ti;Thon_off_init(i)-Ui0(i)*ThPimin(i)];
    
    %pi0-pit<=ui0*(Pdowm+Pmin)-ui1*Pmin+di1*(Pshut-Pdown-Pmin)
    constraints_ramp_up_2P_Ti=sparse(1,5*N*T);
    constraints_ramp_up_2P_Ti(1,(i-1)*T+1)=ThPimin(i);%uit
    constraints_ramp_up_2P_Ti(1,(i-1)*T+N*T+1)=-1;%pit
    constraints_ramp_up_2P_Ti(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
    Aineq_constraints_ramp_up_2P_Ti=[Aineq_constraints_ramp_up_2P_Ti;constraints_ramp_up_2P_Ti];
    bineq_constraints_ramp_up_2P_Ti=[bineq_constraints_ramp_up_2P_Ti;Ui0(i)*(Pidown(i)+ThPimin(i))-Thon_off_init(i)];
    
    %2P_Co
    %pi1-pi0<=ui0*Pup+si1*Pstart
    constraints_ramp_up_2P_Co=sparse(1,5*N*T);
    constraints_ramp_up_2P_Co(1,(i-1)*T+N*T+1)=1;%pit
    constraints_ramp_up_2P_Co(1,(i-1)*T+2*N*T+1)=-1*Pistartup(i);%sit
    Aineq_constraints_ramp_up_2P_Co=[Aineq_constraints_ramp_up_2P_Co;constraints_ramp_up_2P_Co];
    bineq_constraints_ramp_up_2P_Co=[bineq_constraints_ramp_up_2P_Co;Ui0(i)*Piup(i)+Thon_off_init(i)];
    
    %pi0-pi1<=ui1*Pdowm+di1*Pshut
    constraints_ramp_up_2P_Co=sparse(1,5*N*T);
    constraints_ramp_up_2P_Co(1,(i-1)*T+1)=-1*Pidown(i);%uit
    constraints_ramp_up_2P_Co(1,(i-1)*T+N*T+1)=-1;%pit
    constraints_ramp_up_2P_Co(1,(i-1)*T+4*N*T+1)=-1*Pishutdown(i);%dit
    Aineq_constraints_ramp_up_2P_Co=[Aineq_constraints_ramp_up_2P_Co;constraints_ramp_up_2P_Co];
    bineq_constraints_ramp_up_2P_Co=[bineq_constraints_ramp_up_2P_Co;-1*Thon_off_init(i)];
    
    %t>1
    %pit- pit-1 <=uit*(Pup+Pmin)-uit-1 *Pmin+sit*(Pstart-Pup-Pmin)
    constraints_ramp_up_2P_Ti=sparse(T-1,5*N*T);
    constraints_ramp_up_2P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*(Piup(i)+ThPimin(i))*ones(1,T-1)))];%uit
    constraints_ramp_up_2P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(ones(1,T-1)))];%pit
    constraints_ramp_up_2P_Ti(1:T-1,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=sparse(diag(-1*(Pistartup(i)-Piup(i)-ThPimin(i))*ones(1,T-1)));%sit
    Aineq_constraints_ramp_up_2P_Ti=[Aineq_constraints_ramp_up_2P_Ti;constraints_ramp_up_2P_Ti];
    bineq_constraints_ramp_up_2P_Ti=[bineq_constraints_ramp_up_2P_Ti;sparse(T-1,1)];
    
    %pit-1 -pit<=uit-1 *(Pdown+Pmin)-uit*Pmin+dit*(Pshut-Pdown-Pmin)
    constraints_ramp_up_2P_Ti=sparse(T-1,5*N*T);
    constraints_ramp_up_2P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T)=[sparse(diag(-1*(Pidown(i)+ThPimin(i))*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(ThPimin(i)*ones(1,T-1)))];%uit
    constraints_ramp_up_2P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
    constraints_ramp_up_2P_Ti(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown(i)-Pidown(i)-ThPimin(i))*ones(1,T-1)));%dit
    Aineq_constraints_ramp_up_2P_Ti=[Aineq_constraints_ramp_up_2P_Ti;constraints_ramp_up_2P_Ti];
    bineq_constraints_ramp_up_2P_Ti=[bineq_constraints_ramp_up_2P_Ti;sparse(T-1,1)];
    
    %2P_Co
    %t>1
    %pit- pit-1 <=uit-1 *Pup+sit*Pstart
    constraints_ramp_up_2P_Co=sparse(T-1,5*N*T);
    constraints_ramp_up_2P_Co(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*Piup(i)*ones(1,T-1)));%uit
    constraints_ramp_up_2P_Co(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(ones(1,T-1)))];%pit
    constraints_ramp_up_2P_Co(1:T-1,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=sparse(diag(-1*Pistartup(i)*ones(1,T-1)));%sit
    Aineq_constraints_ramp_up_2P_Co=[Aineq_constraints_ramp_up_2P_Co;constraints_ramp_up_2P_Co];
    bineq_constraints_ramp_up_2P_Co=[bineq_constraints_ramp_up_2P_Co;sparse(T-1,1)];
    
    %pit-1 -pit<=uit*Pdowm+dit*Pshut
    constraints_ramp_up_2P_Co=sparse(T-1,5*N*T);
    constraints_ramp_up_2P_Co(1:T-1,(i-1)*T+2:(i-1)*T+T)=sparse(diag(-1*Pidown(i)*ones(1,T-1)));%uit
    constraints_ramp_up_2P_Co(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
    constraints_ramp_up_2P_Co(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*Pishutdown(i)*ones(1,T-1)));%dit
    Aineq_constraints_ramp_up_2P_Co=[Aineq_constraints_ramp_up_2P_Co;constraints_ramp_up_2P_Co];
    bineq_constraints_ramp_up_2P_Co=[bineq_constraints_ramp_up_2P_Co;sparse(T-1,1)];
    
    
    %3P_HD
    %pit+1Wan - pitWan<=uit+1 *PupWan+sit+1 *(PstartWan-PupWan)  t=T-1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+T)=-1*Piup_wan(i);%uit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+T-1)=-1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+T)=1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+T)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;0];
    
    %pit-1wan  -pitWan <=uit-1 *PdownWan +dit*(PshutWan -PdownWan) t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+1)=-1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
    
    %pit-1wan  -pitWan <=uit-1 *PdownWan +dit*(PshutWan -PdownWan) t=2
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+1)=-1*Pidown_wan(i);%uit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+1)=1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+2)=-1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;0];
    
    
    if(ThTime_on_min(i)>1&(Piup(i)>Pishutdown(i)-ThPimin(i)))
        %3P_Ti
        %t=1    pi1-pi0<= ui1*Pup-di1*Pmin-di2*(Pup-Pshut+Pmin)+si1*(Pstart-Pup)
        constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti(1,(i-1)*T+1)=-1*Piup(i);%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+1)=1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i));%sit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+1)=ThPimin(i);%dit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+2)=Piup(i)-Pishutdown(i)+ThPimin(i);%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;Thon_off_init(i)];
        
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;Thon_off_init(i)];
        
        %2<=t<T   pit- pit-1 <=uit*Pup-dit*Pmin-dit+1 *(Pup-Pshut+Pmin)+sit*(Pstart-Pup)
        constraints_ramp_up_3P_Ti=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*Piup(i)*ones(1,T-2)));%uit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup(i)-Piup(i))*ones(1,T-2)));%sit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((Piup(i)-Pishutdown(i)+ThPimin(i))*ones(1,T-2)))];%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;sparse(T-2,1)];
        
        %t=T PiT - PiT-1 <=uiT*(Pup+Pmin)-uiT-1 *Pmin+siT*(Pstart-Pup-Pmin)
        constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti(1,(i-1)*T+T-1)=ThPimin(i);%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+T)=-1*(Piup(i)+ThPimin(i));%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+T-1)=-1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+T)=1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2*N*T+T)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));%sit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;0];
        
        %3P_Ti_ST
        %2<=t<T   pit- pit-1 <=oit*Pup -dit*Pmin-dit+1 *(Pup-Pshut+Pmin)+sit*Pstart
        constraints_ramp_up_3P_Ti_ST=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*Piup(i)*ones(1,T-2)));%uit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup(i))*ones(1,T-2)));%sit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((Piup(i)-Pishutdown(i)+ThPimin(i))*ones(1,T-2)))];%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;sparse(T-2,1)];
        
        %t=T PiT - PiT-1 <=oiT*(Pup+Pmin)- oiT-1 *Pmin- siT-1 *Pmin +siT*Pstart
        constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+T-1)=ThPimin(i);%oit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+T)=-1*(Piup(i)+ThPimin(i));%oit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+T-1)=-1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+T)=1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+T-1)=ThPimin(i);%sit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+T)=-1*Pistartup(i);%sit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;0];
        
        
    else
        %3P_Ti
        %Pit - Pit-1 <=uit*(Pup+Pmin)- uit-1 *Pmin+sit*(Pstart-Pup-Pmin)
        %t=1
        constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti(1,(i-1)*T+1)=-1*(Piup(i)+ThPimin(i));%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+1)=1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));%sit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;Thon_off_init(i)-Ui0(i)*ThPimin(i)];
        
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;Thon_off_init(i)-Ui0(i)*ThPimin(i)];
        
        %t>1
        constraints_ramp_up_3P_Ti=sparse(T-1,5*N*T);
        constraints_ramp_up_3P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*(Piup(i)+ThPimin(i))*ones(1,T-1)))];%uit
        constraints_ramp_up_3P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(ones(1,T-1)))];%pit
        constraints_ramp_up_3P_Ti(1:T-1,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=sparse(diag(-1*(Pistartup(i)-Piup(i)-ThPimin(i))*ones(1,T-1)));%sit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;sparse(T-1,1)];
        
        %3P_Ti_ST
        
        %t>2 Pit - Pit-1 <=oit*(Pup+Pmin)- oit-1 *Pmin- sit-1 *Pmin+sit*Pstart
        constraints_ramp_up_3P_Ti_ST=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2:(i-1)*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Piup(i)+ThPimin(i))*ones(1,T-2)))];%oit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup(i))*ones(1,T-2)))];%sit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;sparse(T-2,1)];
        
        %t=2 %Pi2 - Pi1 <=oi2*(Pup+Pmin)- ui1 *Pmin+si2*Pstart
        constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+1)=ThPimin(i);%uit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2)=-1*(Piup(i)+ThPimin(i));%oit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+1)=-1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+2)=1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+2)=-1*Pistartup(i);%sit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;0];
    end
    
    if(ThTime_on_min(i)>1&(Pidown(i)>Pistartup(i)-ThPimin(i)))
        %3P_Ti
        
        %t=2   pi1-pi2<=ui1*(Pdown+Pmin)-ui2*Pmin+di2*(Pshut-Pdown-Pmin)
        constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti(1,(i-1)*T+1)=-1*(Pidown(i)+ThPimin(i));%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2)=ThPimin(i);%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+1)=1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+2)=-1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;0];
        
        %t>2   pit-1  -pit  <= uit*Pdown+dit*Pshutn- sit-1 *(Pdown-Pstart+Pmin)-sit*(Pdown+Pmin)
        constraints_ramp_up_3P_Ti=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-1*Pidown(i)*ones(1,T-2)));%uit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag((Pidown(i)-Pistartup(i)+ThPimin(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((Pidown(i)+ThPimin(i))*ones(1,T-2)))];%sit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag(-1*Pishutdown(i)*ones(1,T-2)));%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;sparse(T-2,1)];
        
        %3P_Ti_ST
        %t=2   pi1-pi2<=ui1*(Pdown+Pmin)-oi2*Pmin-si2*Pmin+di2*(Pshut-Pdown-Pmin)
        constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+1)=-1*(Pidown(i)+ThPimin(i));%uit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2)=ThPimin(i);%oit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+1)=1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+2)=-1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+2)=ThPimin(i);%sit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;0];
        
        %t>2   pit-1  -pit  <= oit*Pdown+dit*Pshutn- sit-1 *(Pdown-Pstart+Pmin)-sit*Pmin
        constraints_ramp_up_3P_Ti_ST=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-1*Pidown(i)*ones(1,T-2)));%oit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag((Pidown(i)-Pistartup(i)+ThPimin(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((ThPimin(i))*ones(1,T-2)))];%sit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag(-1*Pishutdown(i)*ones(1,T-2)));%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;sparse(T-2,1)];
        
    else
        %3P_Ti
        %t=1   pi0-pi1<= ui0*(Pdown+Pmin)-ui1*Pmin+di1*(Pshut-Pdown-Pmin)
        constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti(1,(i-1)*T+1)=ThPimin(i);%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+1)=-1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;Ui0(i)*(Pidown(i)+ThPimin(i))-Thon_off_init(i)];
        
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;Ui0(i)*(Pidown(i)+ThPimin(i))-Thon_off_init(i)];
        
        %t>1   pit-1 -pit <= uit-1 *(Pdown+Pmin)-uit*Pmin+dit*(Pshut-Pdown-Pmin)
        constraints_ramp_up_3P_Ti=sparse(T-1,5*N*T);
        constraints_ramp_up_3P_Ti(1:T-1,(i-1)*T+1:(i-1)*T+T)=[sparse(diag(-1*(Pidown(i)+ThPimin(i))*ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(ThPimin(i)*ones(1,T-1)))];%uit
        constraints_ramp_up_3P_Ti(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
        constraints_ramp_up_3P_Ti(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown(i)-Pidown(i)-ThPimin(i))*ones(1,T-1)));%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;sparse(T-1,1)];
        
        %3P_Ti_ST
        %t=2   pi1 -pi2 <= ui1 *(Pdown+Pmin)-oi2*Pmin-si2*Pmin+di2*(Pshut-Pdown-Pmin)
        constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+1)=-1*(Pidown(i)+ThPimin(i));%uit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2)=ThPimin(i);%oit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+1)=1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+2)=-1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+2)=ThPimin(i);%sit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;0];
        
        %t>2   pit-1 -pit <= oit-1 *(Pdown+Pmin)+ sit-1 *(Pdown+Pmin)-oit*Pmin-sit*Pmin+dit*(Pshut-Pdown-Pmin)
        constraints_ramp_up_3P_Ti_ST=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2:(i-1)*T+T)=[sparse(diag(-1*(Pidown(i)+ThPimin(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ThPimin(i)*ones(1,T-2)))];%oit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+N*T+2:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pidown(i)+ThPimin(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ThPimin(i)*ones(1,T-2)))];%sit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown(i)-Pidown(i)-ThPimin(i))*ones(1,T-2)));%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;sparse(T-2,1)];
        
    end
    
    
    
    if(ThTime_on_min(i)>1&(Piup(i)>Pishutdown(i)-ThPimin(i))&ThTime_off_min(i)>=2)
        %3P_Ti
        %t=1  pi2 -pi0<=2*ui2*Pup -di1*Pmin -di2*Pmin+si1*(Pstart-Pup)+si2*(Pstart-2*Pup)
        constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2)=-2*Piup(i);%uit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+2)=1;%pit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i));%sit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+2*N*T+2)=-1*(Pistartup(i)-2*Piup(i));%sit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+1)=ThPimin(i);%dit
        constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+2)=ThPimin(i);%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;Thon_off_init(i)];
        
        
        %t>1  pit+1 - pit-1 <=2*uit+1 *Pup -dit*Pmin - dit+1 *Pmin+sit*(Pstart-Pup)+ sit+1 *(Pstart-2*Pup)
        constraints_ramp_up_3P_Ti=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*Piup(i)*ones(1,T-2)));%uit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup(i)-Piup(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup(i)-2*Piup(i))*ones(1,T-2)))];%sit
        constraints_ramp_up_3P_Ti(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ThPimin(i)*ones(1,T-2)))];%dit
        Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
        bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;sparse(T-2,1)];
        
        %3P_Ti_ST
        %t=1  pi2 -pi0<=2*oi2*Pup -di1*Pmin -di2*Pmin+si1*(Pstart-Pup)+si2*(Pstart)
        constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2)=-2*Piup(i);%oit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+2)=1;%pit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i));%sit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+2)=-1*(Pistartup(i));%sit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+4*N*T+1)=ThPimin(i);%dit
        constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+4*N*T+2)=ThPimin(i);%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;Thon_off_init(i)];
        
        %t>1  pit+1 - pit-1 <=2*oit+1 *Pup -dit*Pmin - dit+1 *Pmin+sit*(Pstart-Pup)+ sit+1 *Pstart
        constraints_ramp_up_3P_Ti_ST=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*Piup(i)*ones(1,T-2)));%oit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup(i)-Piup(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup(i))*ones(1,T-2)))];%sit
        constraints_ramp_up_3P_Ti_ST(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(ThPimin(i)*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ThPimin(i)*ones(1,T-2)))];%dit
        Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
        bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;sparse(T-2,1)];
        
    end
    
    %3P_Ti
    %pi1-pi0<=ui1*(Pup+Pmin)-ui0*Pmin+si1*(Pstart-Pup-Pmin)
%     constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
%     constraints_ramp_up_3P_Ti(1,(i-1)*T+1)=-1*(Piup(i)+ThPimin(i));%uit
%     constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+1)=1;%pit
%     constraints_ramp_up_3P_Ti(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));%sit
%     Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
%     bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;Thon_off_init(i)-Ui0(i)*ThPimin(i)];
    
    %pi0-pit<=ui0*(Pdowm+Pmin)-ui1*Pmin+di1*(Pshut-Pdown-Pmin)
%     constraints_ramp_up_3P_Ti=sparse(1,5*N*T);
%     constraints_ramp_up_3P_Ti(1,(i-1)*T+1)=ThPimin(i);%uit
%     constraints_ramp_up_3P_Ti(1,(i-1)*T+N*T+1)=-1;%pit
%     constraints_ramp_up_3P_Ti(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
%     Aineq_constraints_ramp_up_3P_Ti=[Aineq_constraints_ramp_up_3P_Ti;constraints_ramp_up_3P_Ti];
%     bineq_constraints_ramp_up_3P_Ti=[bineq_constraints_ramp_up_3P_Ti;Ui0(i)*(Pidown(i)+ThPimin(i))-Thon_off_init(i)];
    
    %3P_Ti_ST
    %pi1-pi0<=ui1*(Pup+Pmin)-ui0*Pmin+si1*(Pstart-Pup-Pmin)
%     constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
%     constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+1)=-1*(Piup(i)+ThPimin(i));%uit
%     constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+1)=1;%pit
%     constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+2*N*T+1)=-1*(Pistartup(i)-Piup(i)-ThPimin(i));%sit
%     Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
%     bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;Thon_off_init(i)-Ui0(i)*ThPimin(i)];
    
    %pi0-pit<=ui0*(Pdowm+Pmin)-ui1*Pmin+di1*(Pshut-Pdown-Pmin)
%     constraints_ramp_up_3P_Ti_ST=sparse(1,5*N*T);
%     constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+1)=ThPimin(i);%uit
%     constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+N*T+1)=-1;%pit
%     constraints_ramp_up_3P_Ti_ST(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown(i)-Pidown(i)-ThPimin(i));%dit
%     Aineq_constraints_ramp_up_3P_Ti_ST=[Aineq_constraints_ramp_up_3P_Ti_ST;constraints_ramp_up_3P_Ti_ST];
%     bineq_constraints_ramp_up_3P_Ti_ST=[bineq_constraints_ramp_up_3P_Ti_ST;Ui0(i)*(Pidown(i)+ThPimin(i))-Thon_off_init(i)];
    
    
    
    
    if(ThTime_on_min(i)>1)
        %3P_HD
    %pitWan - pit-1Wan <=uit*PupWan+sit*(PstartWan-PupWan)-dit+1*(max(PupWan-PshutWan,0))
    %t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+1)=-1*Piup_wan(i);%uit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+1)=1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+2)=max(Piup_wan(i)-Pishutdown_wan(i),0);%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;ThPit_wani0(i)];
    %t>1
    constraints_ramp_up_3P_HD=sparse(T-2,6*N*T);
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(Piup_wan(i))*ones(1,T-2)));%uit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2)));%sit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Piup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-2,1)];
    
    %pitWan-pit+1Wan<=uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(PshutWan-PdownWan)
    %t>=1
    constraints_ramp_up_3P_HD=sparse(T-1,6*N*T);
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i))*ones(1,T-1)));%uit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pidown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-1,1)];
    
    %pit+1Wan - pit-1Wan <=uit+1*2PupWan+sit*(PstartWan-PupWan)+sit+1*(PstartWan -2PupWan)
    %t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+2)=1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;ThPit_wani0(i)];
    %t>1
    constraints_ramp_up_3P_HD=sparse(T-2,6*N*T);
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)));%uit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)))];%sit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-2,1)];
    
    %pit-1Wan -pit+1Wan <= uit-1 *2*PdownWan- dit*(PshutWan-2*PdownWan)+ dit+1 *(PshutWan-PdownWan)
    %t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+2)=-1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+1)=Pishutdown_wan(i)-2*Pidown_wan(i);%dit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
    %t>1
    constraints_ramp_up_3P_HD=sparse(T-2,6*N*T);
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2)));%uit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag((Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-2,1)];
        
        
        %3P_HD_Pr
        %pitWan - pit-1Wan <=uit*PupWan+sit*(PstartWan-PupWan)-dit+1*(max(PupWan-PshutWan,0))
        %t=1
        constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=-1*Piup_wan(i);%uit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=max(Piup_wan(i)-Pishutdown_wan(i),0);%dit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
        %t>1
        constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(Piup_wan(i))*ones(1,T-2)));%uit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2)));%sit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Piup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        
        %pitWan-pit+1Wan<=uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(PshutWan-PdownWan)
        %     %t>=1
        constraints_ramp_up_3P_HD_Pr=sparse(T-1,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i))*ones(1,T-1)));%uit
        constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
        constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pidown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
        constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-1,1)];
        

        %pit+1Wan - pit-1Wan <=uit+1*2PupWan+sit*(PstartWan-PupWan)+sit+1*(PstartWan -2PupWan)
        %t=1
        constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=1;%pit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
        %t>1
        constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)));%uit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)))];%sit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        
        %pit-1Wan -pit+1Wan <=uit-1 *2*PdownWan- dit*(PshutWan-2*PdownWan)+ dit+1 *(PshutWan-PdownWan)
        %t=1
        constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=-1;%pit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=Pishutdown_wan(i)-2*Pidown_wan(i);%dit
        constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
        %t>1
        constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2)));%uit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
        constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag((Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
        Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
        bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        
        
    else
        %3P_HD
    %pitWan - pit-1Wan <=T3*(max(PupWan-PshutWan,0)-max(PstartWan-PshutWan,0))+uit*PupWan+sit*(PstartWan-PupWan)-dit+1*(max(PupWan-PshutWan,0))
    %t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+1)=-1*Piup_wan(i);%uit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+1)=1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+2)=max(Piup_wan(i)-Pishutdown_wan(i),0);%dit
    constraints_ramp_up_3P_HD(1,(i-1)*T+5*N*T+1)=-1*(max(Piup_wan(i)-Pishutdown_wan(i),0)-max(Pistartup_wan(i)-Pishutdown_wan(i),0));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;ThPit_wani0(i)];
    %t>1
    constraints_ramp_up_3P_HD=sparse(T-2,6*N*T);
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(Piup_wan(i))*ones(1,T-2)));%uit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2)));%sit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Piup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+5*N*T+2:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(max(Piup_wan(i)-Pishutdown_wan(i),0)-max(Pistartup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-2,1)];
    
    %pitWan-pit+1Wan<=T3*(max(PdownWan-PstartWan,0)-max(PshutWan-PstartWan,0))+uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(PshutWan-PdownWan)
    %t>=1
    constraints_ramp_up_3P_HD=sparse(T-1,6*N*T);
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i))*ones(1,T-1)));%uit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pidown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
    constraints_ramp_up_3P_HD(1:T-1,(i-1)*T+5*N*T+1:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(max(Pidown_wan(i)-Pistartup_wan(i),0)-max(Pishutdown_wan(i)-Pistartup_wan(i),0))*ones(1,T-1)));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-1,1)];
    
    %pit+1Wan - pit-1Wan <=T3*(PupWan-PstartWan)+uit+1*2PupWan+sit*(PstartWan-PupWan)+sit+1*(PstartWan -2PupWan)
    %t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+2)=1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
    constraints_ramp_up_3P_HD(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
    constraints_ramp_up_3P_HD(1,(i-1)*T+5*N*T+1)=-1*(Piup_wan(i)-Pistartup_wan(i));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;ThPit_wani0(i)];
    %t>1
    constraints_ramp_up_3P_HD=sparse(T-2,6*N*T);
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)));%uit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)))];%sit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+5*N*T+2:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(Piup_wan(i)-Pistartup_wan(i))*ones(1,T-2)));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-2,1)];
    
    %pit-1Wan -pit+1Wan <=T3*(PdownWan - PshutWan)+ uit-1 *2*PdownWan+ dit*(PshutWan-2*PdownWan)+ dit+1 *(PshutWan-PdownWan)
    %t=1
    constraints_ramp_up_3P_HD=sparse(1,6*N*T);
    constraints_ramp_up_3P_HD(1,(i-1)*T+N*T+2)=-1;%pit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+1)=Pishutdown_wan(i)-2*Pidown_wan(i);%dit
    constraints_ramp_up_3P_HD(1,(i-1)*T+4*N*T+2)=(Pishutdown_wan(i)-Pidown_wan(i));%dit
    constraints_ramp_up_3P_HD(1,(i-1)*T+5*N*T+1)=-1*(Pidown_wan(i)-Pishutdown_wan(i));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
    %t>1
    constraints_ramp_up_3P_HD=sparse(T-2,6*N*T);
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2)));%uit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(-1*(Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
    constraints_ramp_up_3P_HD(1:T-2,(i-1)*T+5*N*T+2:(i-1)*T+5*N*T+T-1)=sparse(diag(-1*(Pidown_wan(i)-Pishutdown_wan(i))*ones(1,T-2)));%T3
    Aineq_constraints_ramp_up_3P_HD=[Aineq_constraints_ramp_up_3P_HD;constraints_ramp_up_3P_HD];
    bineq_constraints_ramp_up_3P_HD=[bineq_constraints_ramp_up_3P_HD;sparse(T-2,1)];
        
        %3P_HD_Pr
        %pitWan - pit-1Wan <=T3*(max(PupWan-PshutWan,0)-max(PstartWan-PshutWan,0))+uit*PupWan+sit*(PstartWan-PupWan)-dit+1*(max(PupWan-PshutWan,0))
        if(max(Piup_wan(i)-Pishutdown_wan(i),0)-max(Pistartup_wan(i)-Pishutdown_wan(i),0)>=0)
            %pitWan - pit-1Wan <=sit*(min(PstartWan,PshutWan)-min(PupWan,PshutWan))+uit*PupWan-dit+1*(max(PupWan-PshutWan,0))
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=-1*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(min(Pistartup_wan(i),Pishutdown_wan(i))-min(Piup_wan(i),Pishutdown_wan(i)));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=max(Piup_wan(i)-Pishutdown_wan(i),0);%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(Piup_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(min(Pistartup_wan(i),Pishutdown_wan(i))-min(Piup_wan(i),Pishutdown_wan(i)))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Piup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
            %pitWan - pit-1Wan <=sit*(PstartWan-PupWan)+uit*PupWan-dit+1*(max(PstartWan-PshutWan,0))
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=-1*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=max(Pistartup_wan(i)-Pishutdown_wan(i),0);%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(Piup_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Pistartup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
        else
            %pitWan - pit-1Wan <=uit*PupWan+sit*(PstartWan-PupWan)-dit+1*(max(PupWan-PshutWan,0))
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=-1*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=max(Piup_wan(i)-Pishutdown_wan(i),0);%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(Piup_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Piup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
            %pitWan - pit-1Wan <=sit*(min(PstartWan,PshutWan)-min(PupWan,PshutWan))+uit*(min(PupWan,PshutWan)+max(PstartWan-PshutWan,0))-dit+1*(max(PstartWan-PshutWan,0))
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=-1*(min(Piup_wan(i),Pishutdown_wan(i))+max(Pistartup_wan(i)-Pishutdown_wan(i),0));%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(min(Pistartup_wan(i),Pishutdown_wan(i))-min(Piup_wan(i),Pishutdown_wan(i)));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=max(Pistartup_wan(i)-Pishutdown_wan(i),0);%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2:(i-1)*T+T-1)=sparse(diag(-1*(min(Piup_wan(i),Pishutdown_wan(i))+max(Pistartup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T-1)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(min(Pistartup_wan(i),Pishutdown_wan(i))-min(Piup_wan(i),Pishutdown_wan(i)))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag((max(Pistartup_wan(i)-Pishutdown_wan(i),0))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
        end
        
        %3P_HD_Pr
        %pitWan-pit+1Wan<=T3*(max(PdownWan-PstartWan,0)-max(PshutWan-PstartWan,0))+uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(PshutWan-PdownWan)
        if((max(Pidown_wan(i)-Pistartup_wan(i),0)-max(Pishutdown_wan(i)-Pistartup_wan(i),0))>=0)
            %          %pitWan- pit+1Wan <=uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(min(PshutWan,PstartWan)-min(PstartWan,PdownWan))
            %pitWan-pit+1Wan<=uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(max(PdownWan-PstartWan,0)-max(PshutWan-PstartWan,0)+PshutWan-PdownWan)
            %         %t>=1
            constraints_ramp_up_3P_HD_Pr=sparse(T-1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i))*ones(1,T-1)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pidown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(max(Pidown_wan(i)-Pistartup_wan(i),0)-max(Pishutdown_wan(i)-Pistartup_wan(i),0)+Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-1,1)];
            %
            %pitWan-pit+1Wan<=uit*PdownWan-sit*max(PshutWan-PstartWan,0)+dit+1*(PshutWan-PdownWan)
            %t>=1
            constraints_ramp_up_3P_HD_Pr=sparse(T-1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i))*ones(1,T-1)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pishutdown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-1,1)];
            
            
        else
            %pitWan-pit+1Wan<=uit*PdownWan-sit*(max(PdownWan-PstartWan,0))+dit+1*(PshutWan-PdownWan)
            %t>=1
            constraints_ramp_up_3P_HD_Pr=sparse(T-1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i))*ones(1,T-1)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pidown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-1,1)];
            
            %pitWan-pit+1Wan<=uit*(min(PstartWan,PdownWan)+max(PshutWan-PstartWan,0))-sit*(max(PshutWan-PstartWan,0))+dit+1*(min(PshutWan,PstartWan)-min(PdownWan,Pstart))
            %pitWan-pit+1Wan<=uit*(-max(PdownWan-PstartWan,0)+max(PshutWan-PstartWan,0)+PdownWan)+sit*(-max(PshutWan-PstartWan,0))+dit+1 *(max(PdownWan-PstartWan,0)-max(PshutWan-PstartWan,0)+PshutWan-PdownWan)
            constraints_ramp_up_3P_HD_Pr=sparse(T-1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+1:(i-1)*T+T-1)=sparse(diag(-1*(Pidown_wan(i)-max(Pidown_wan(i)-Pistartup_wan(i),0)+max(Pishutdown_wan(i)-Pistartup_wan(i),0))*ones(1,T-1)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-1))),sparse(T-1,1)]+[sparse(T-1,1),sparse(diag(-1*ones(1,T-1)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T-1)=sparse(diag(max(Pishutdown_wan(i)-Pistartup_wan(i),0)*ones(1,T-1)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-1,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=sparse(diag(-1*(max(Pidown_wan(i)-Pistartup_wan(i),0)-max(Pishutdown_wan(i)-Pistartup_wan(i),0)+Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-1)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-1,1)];
        end
        %3P_HD_Pr
        %pit+1Wan - pit-1Wan <=T3*(PupWan-PstartWan)+uit+1*2PupWan+sit*(PstartWan-PupWan)+sit+1*(PstartWan -2PupWan)
        if(Piup_wan(i)-Pistartup_wan(i)>=0)
            %pit+1Wan - pit-1Wan <=uit+1 *2PupWan+sit+1 *(PstartWan -2PupWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+3:(i-1)*T+2*N*T+T)=sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)));%sit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
            %pit+1Wan - pit-1Wan<=uit+1*2PupWan+sit*(PstartWan-PupWan)+sit+1*(PstartWan -2PupWan)+ dit+1(PupWan-PstartWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=-1*(Piup_wan(i)-Pistartup_wan(i));
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)))];%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T+T)=sparse(diag(-1*(Piup_wan(i)-Pistartup_wan(i))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        else
            %pit+1Wan - pit-1Wan <=uit+1*2PupWan+sit*(PstartWan-PupWan)+sit+1*(PstartWan -2PupWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+3:(i-1)*T+T)=sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T)=[sparse(diag(-1*(Pistartup_wan(i)-Piup_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)))];%sit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
            
            %pit+1Wan - pit-1Wan<=-uit*(PupWan-PstartWan)+dit+1*(PupWan-PstartWan)+ uit+1*2Pup+ sit+1*(PstartWan-2*PupWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=Piup_wan(i)-Pistartup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2)=-2*Piup_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+2)=-1*(Pistartup_wan(i)-2*Piup_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=-1*(Piup_wan(i)-Pistartup_wan(i));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2:(i-1)*T+T)=[sparse(diag((Piup_wan(i)-Pistartup_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-2*(Piup_wan(i))*ones(1,T-2)))];%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(-1*ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+3:(i-1)*T+2*N*T+T)=sparse(diag(-1*(Pistartup_wan(i)-2*Piup_wan(i))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+3:(i-1)*T+4*N*T++T)=sparse(diag(-1*(Piup_wan(i)-Pistartup_wan(i))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        end
        %3P_HD_Pr
        %pit-1Wan -pit+1Wan <=T3*(PdownWan - PshutWan)+ uit-1 *2*PdownWan+ dit*(PshutWan-2*PdownWan)+ dit+1 *(PshutWan-PdownWan)
        if(Pidown_wan(i)-Pishutdown_wan(i)>=0)
            
            %pit-1Wan -pit+1Wan <=sit*(PdownWan - PshutWan)+ uit-1 *2*PdownWan+ dit*(PshutWan-2*PdownWan)+ dit+1 *(PshutWan-PdownWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=-1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pidown_wan(i)-Pishutdown_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown_wan(i)-2*Pidown_wan(i));%dit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pidown_wan(i)-Pishutdown_wan(i))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(-1*(Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
            %pit-1Wan -pit+1Wan <= uit-1 *2*PdownWan+ dit*(PshutWan-2*PdownWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=-1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown_wan(i)-2*Pidown_wan(i));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T-1)=sparse(diag(-1*(Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        else
            %pit-1Wan -pit+1Wan <= uit-1 *2*PdownWan+ dit*(PshutWan-2*PdownWan)+ dit+1 *(PshutWan-PdownWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=-1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown_wan(i)-2*Pidown_wan(i));%dit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-2)=sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2)));%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T)=[sparse(diag(-1*(Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag(-1*(Pishutdown_wan(i)-Pidown_wan(i))*ones(1,T-2)))];%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
            
            %pit-1Wan -pit+1Wan <=sit*(PdownWan - PshutWan)-uit*(PdownWan -PshutWan) +uit-1 *2*PdownWan+ dit*(PshutWan-2*PdownWan)
            %t=1
            constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=Pidown_wan(i)-Pishutdown_wan(i);%uit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=-1;%pit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+1)=-1*(Pidown_wan(i)-Pishutdown_wan(i));%sit
            constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown_wan(i)-2*Pidown_wan(i));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;2*Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
            %t>1
            constraints_ramp_up_3P_HD_Pr=sparse(T-2,5*N*T);
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+1:(i-1)*T+T-1)=[sparse(diag(-2*(Pidown_wan(i))*ones(1,T-2))),sparse(T-2,1)]+[sparse(T-2,1),sparse(diag((Pidown_wan(i)-Pishutdown_wan(i))*ones(1,T-2)))];%uit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+N*T+1:(i-1)*T+N*T+T)=[sparse(diag(ones(1,T-2))),sparse(T-2,2)]+[sparse(T-2,2),sparse(diag(-1*ones(1,T-2)))];%pit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+2*N*T+2:(i-1)*T+2*N*T+T-1)=sparse(diag(-1*(Pidown_wan(i)-Pishutdown_wan(i))*ones(1,T-2)));%sit
            constraints_ramp_up_3P_HD_Pr(1:T-2,(i-1)*T+4*N*T+2:(i-1)*T+4*N*T+T-1)=sparse(diag(-1*(Pishutdown_wan(i)-2*Pidown_wan(i))*ones(1,T-2)));%dit
            Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
            bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;sparse(T-2,1)];
        end
        
    end
    
    %3P_HD_Pr
    %pit+1Wan - pitWan<=uit+1 *PupWan+sit+1 *(PstartWan-PupWan)  t=T-1
    constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
    constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+T)=-1*Piup_wan(i);%uit
    constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+T-1)=-1;%pit
    constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+T)=1;%pit
    constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+2*N*T+T)=-1*(Pistartup_wan(i)-Piup_wan(i));%sit
    Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
    bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;0];
    
    %pit-1wan  -pitWan <=uit-1 *PdownWan +dit*(PshutWan -PdownWan) t=1
    constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
    constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=-1;%pit
    constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+1)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
    Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
    bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;Ui0(i)*Pidown_wan(i)-ThPit_wani0(i)];
    
    %pit-1wan  -pitWan <=uit-1 *PdownWan +dit*(PshutWan -PdownWan) t=2
    % constraints_ramp_up_3P_HD_Pr=sparse(1,5*N*T);
    % constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+1)=-1*Pidown_wan(i);%uit
    % constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+1)=1;%pit
    % constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+N*T+2)=-1;%pit
    % constraints_ramp_up_3P_HD_Pr(1,(i-1)*T+4*N*T+2)=-1*(Pishutdown_wan(i)-Pidown_wan(i));%dit
    % Aineq_constraints_ramp_up_3P_HD_Pr=[Aineq_constraints_ramp_up_3P_HD_Pr;constraints_ramp_up_3P_HD_Pr];
    % bineq_constraints_ramp_up_3P_HD_Pr=[bineq_constraints_ramp_up_3P_HD_Pr;0];
    
end


%% Power balance constraints

pit=[];
for t=1:N
    pit=sparse([pit,eye(T)]);
end
Aeq_constraint_power_balance=[sparse(T,N*T),pit];
beq_constraint_power_balance=Dt;

uit=[];
pitWan=[];
for i=1:N
    uit=[uit,eye(T)*ThPimin(i)];
    pitWan=[pitWan,eye(T)*(ThPimax(i) - ThPimin(i))];
end
Aeq_constraint_power_balance_project=[uit,pitWan];
beq_constraint_power_balance_project=Dt;

Aeq_constraint_power_balance_2P_Co=[Aeq_constraint_power_balance,sparse(size(Aeq_constraint_power_balance,1),3*N*T)];
beq_constraint_power_balance_2P_Co=beq_constraint_power_balance;

Aeq_constraint_power_balance_2P_Ti=[Aeq_constraint_power_balance,sparse(size(Aeq_constraint_power_balance,1),3*N*T)];
beq_constraint_power_balance_2P_Ti=beq_constraint_power_balance;

Aeq_constraint_power_balance_3P_Ti=[Aeq_constraint_power_balance,sparse(size(Aeq_constraint_power_balance,1),3*N*T)];
beq_constraint_power_balance_3P_Ti=beq_constraint_power_balance;

Aeq_constraint_power_balance_3P_Ti_ST=[Aeq_constraint_power_balance,sparse(size(Aeq_constraint_power_balance,1),3*N*T)];
beq_constraint_power_balance_3P_Ti_ST=beq_constraint_power_balance;

Aeq_constraint_power_balance_3P_HD=[Aeq_constraint_power_balance_project,sparse(size(Aeq_constraint_power_balance_project,1),4*N*T)];
beq_constraint_power_balance_3P_HD=beq_constraint_power_balance_project;

Aeq_constraint_power_balance_3P_HD_Pr=[Aeq_constraint_power_balance_project,sparse(size(Aeq_constraint_power_balance_project,1),3*N*T)];
beq_constraint_power_balance_3P_HD_Pr=beq_constraint_power_balance_project;
%% system spinning reserve requirement
Aineq_constraints_spinning_reserve_requirment=[];
bineq_constraints_spinning_reserve_requirment=-Dt-Spin;
uit=[];
sit=[];
constraints_spinning_reserve_requirment=sparse(1,5*N*T);
for i=1:N
    uit=[uit,-1*eye(T)*ThPimax(i)];
    sit=[sit,[sparse(1,T);sparse(T-1,1),-1*eye(T-1)*ThPimax(i)]];
end
Aineq_constraints_spinning_reserve_requirment=[uit,sparse(T,4*N*T)];

Aineq_constraints_spinning_reserve_requirment_2P_Co=Aineq_constraints_spinning_reserve_requirment;
bineq_constraints_spinning_reserve_requirment_2P_Co=bineq_constraints_spinning_reserve_requirment;

Aineq_constraints_spinning_reserve_requirment_2P_Ti=Aineq_constraints_spinning_reserve_requirment;
bineq_constraints_spinning_reserve_requirment_2P_Ti=bineq_constraints_spinning_reserve_requirment;

Aineq_constraints_spinning_reserve_requirment_3P_Ti=Aineq_constraints_spinning_reserve_requirment;
bineq_constraints_spinning_reserve_requirment_3P_Ti=bineq_constraints_spinning_reserve_requirment;

Aineq_constraints_spinning_reserve_requirment_3P_Ti_ST=[uit,sparse(T,N*T),sit,sparse(T,2*N*T)];
bineq_constraints_spinning_reserve_requirment_3P_Ti_ST=bineq_constraints_spinning_reserve_requirment;

Aineq_constraints_spinning_reserve_requirment_3P_HD=[Aineq_constraints_spinning_reserve_requirment,sparse(T,N*T)];
bineq_constraints_spinning_reserve_requirment_3P_HD=bineq_constraints_spinning_reserve_requirment;

Aineq_constraints_spinning_reserve_requirment_3P_HD_Pr=Aineq_constraints_spinning_reserve_requirment;
bineq_constraints_spinning_reserve_requirment_3P_HD_Pr=bineq_constraints_spinning_reserve_requirment;
%% CET
if(CET==1)
    variable_CET={'uit/oit','B',N*T;'Pit/PitWan','C',N*T;'sit','B',N*T;'Sit/SitWan','C',N*T;'dit','B',N*T;'T3','B',N*T;'E_b','C',1;'E_s','C',1;'u_b','B',1;'u_s','B',1;};
    variable_num_3P_HD=6*N*T+4;
    variable_num_3P_HD_Pr=5*N*T+4;
    E0=dataUC.E0;
    a=dataUC.a;
    b=dataUC.b;
    c=dataUC.c;
    price_buy=dataUC.price_buy;
    price_sell=dataUC.price_sell;
    emmission_buy_max=dataUC.emmission_buy_max;
    emmission_sell_max=dataUC.emmission_sell_max;
    
     %3P_HD
    
    Aineq_emission_buy_up=sparse(1,variable_num_3P_HD);
    Aineq_emission_buy_up(1,6*N*T+1)=1;
    Aineq_emission_buy_up(1,6*N*T+3)=-1*emmission_buy_max;
    bineq_emission_buy_up=0;
    
    Aineq_emission_sell_up=sparse(1,variable_num_3P_HD);
    Aineq_emission_sell_up(1,6*N*T+2)=1;
    Aineq_emission_sell_up(1,6*N*T+4)=-1*emmission_sell_max;
    bineq_emission_sell_up=0;
    
    Aineq_sell_buy_constraint=sparse(1,variable_num_3P_HD);
    Aineq_sell_buy_constraint(1,6*N*T+3)=1;
    Aineq_sell_buy_constraint(1,6*N*T+4)=1;
    bineq_sell_buy_constraint=1;
    
    Aineq_trade_constraint_3P_HD=[Aineq_emission_sell_up;Aineq_emission_buy_up;Aineq_sell_buy_constraint];
    bineq_trade_constraint_3P_HD=[bineq_emission_sell_up;bineq_emission_buy_up;bineq_sell_buy_constraint];
    
    Q_emission=sparse(variable_num_3P_HD,variable_num_3P_HD);
    l_emission=sparse(1,variable_num_3P_HD);
    r_emission=E0;
    
    Pwan_emission=sparse(N*T,variable_num_3P_HD);
    Pwan_emission(1:N*T,N*T+1:2*N*T)=diag(reshape(repmat(c,1,T)',N*T,1));
    Q_emission=sparse(variable_num_3P_HD,variable_num_3P_HD);
    Q_emission(N*T+1:2*N*T,:)=Pwan_emission;
    
    l_emission=sparse(variable_num_3P_HD,1);
    l_emission(1:N*T,1)=reshape(repmat(a,1,T)',N*T,1);%uit
    l_emission(N*T+1:2*N*T,1)=reshape(repmat(b,1,T)',N*T,1);%pit
    l_emission(6*N*T+1,1)=-1;%E_b
    l_emission(6*N*T+2,1)=1;%E_s
   
    Q_3P_HD=Q_emission;
    l_3P_HD=l_emission;
    r_3P_HD=r_emission;
    
    %3P_HD_Pr

    Aineq_emission_buy_up=sparse(1,variable_num_3P_HD_Pr);
    Aineq_emission_buy_up(1,5*N*T+1)=1;
    Aineq_emission_buy_up(1,5*N*T+3)=-1*emmission_buy_max;
    bineq_emission_buy_up=0;
    
    Aineq_emission_sell_up=sparse(1,variable_num_3P_HD_Pr);
    Aineq_emission_sell_up(1,5*N*T+2)=1;
    Aineq_emission_sell_up(1,5*N*T+4)=-1*emmission_sell_max;
    bineq_emission_sell_up=0;
    
    Aineq_sell_buy_constraint=sparse(1,variable_num_3P_HD_Pr);
    Aineq_sell_buy_constraint(1,5*N*T+3)=1;
    Aineq_sell_buy_constraint(1,5*N*T+4)=1;
    bineq_sell_buy_constraint=1;
    
    Aineq_trade_constraint_3P_HD_Pr=[Aineq_emission_sell_up;Aineq_emission_buy_up;Aineq_sell_buy_constraint];
    bineq_trade_constraint_3P_HD_Pr=[bineq_emission_sell_up;bineq_emission_buy_up;bineq_sell_buy_constraint];
    
    Q_emission=sparse(variable_num_3P_HD_Pr,variable_num_3P_HD_Pr);
    l_emission=sparse(1,variable_num_3P_HD_Pr);
    r_emission=E0;
    
    Pwan_emission=sparse(N*T,variable_num_3P_HD_Pr);
    Pwan_emission(1:N*T,N*T+1:2*N*T)=diag(reshape(repmat(c,1,T)',N*T,1));
    Q_emission=sparse(variable_num_3P_HD_Pr,variable_num_3P_HD_Pr);
    Q_emission(N*T+1:2*N*T,:)=Pwan_emission;
    
    l_emission=sparse(variable_num_3P_HD_Pr,1);
    l_emission(1:N*T,1)=reshape(repmat(a,1,T)',N*T,1);%uit
    l_emission(N*T+1:2*N*T,1)=reshape(repmat(b,1,T)',N*T,1);%pit
    l_emission(5*N*T+1,1)=-1;%E_b
    l_emission(5*N*T+2,1)=1;%E_s
    
    
    Q_3P_HD_Pr=Q_emission;
    l_3P_HD_Pr=l_emission;
    r_3P_HD_Pr=r_emission;
end





%% objective
%二次项
Pwan_coefficient=sparse(N*T,5*N*T);
Pwan_coefficient(1:N*T,N*T+1:2*N*T)=2*diag(reshape(repmat(Gama,1,T)',N*T,1));
H=sparse(5*N*T,5*N*T);
H(N*T+1:2*N*T,:)=Pwan_coefficient;

Pwan_coefficient_wan=sparse(N*T,6*N*T);
Pwan_coefficient_wan(1:N*T,N*T+1:2*N*T)=2*diag(reshape(repmat(Gama_wan,1,T)',N*T,1));
H_wan=sparse(6*N*T,6*N*T);
H_wan(N*T+1:2*N*T,:)=Pwan_coefficient_wan;

Pwan_coefficient_wan_Pr=sparse(N*T,5*N*T);
Pwan_coefficient_wan_Pr(1:N*T,N*T+1:2*N*T)=2*diag(reshape(repmat(Gama_wan,1,T)',N*T,1));
H_wan_Pr=sparse(5*N*T,5*N*T);
H_wan_Pr(N*T+1:2*N*T,:)=Pwan_coefficient_wan_Pr;

%一次项
f=sparse(5*N*T,1);
f(1:N*T,1)=reshape(repmat(Alpha,1,T)',N*T,1);%uit
f(N*T+1:2*N*T,1)=reshape(repmat(Beta,1,T)',N*T,1);%pit
f(3*N*T+1:4*N*T,1)=ones(N*T,1);%Sit

f_3P_Hi_ST=sparse(5*N*T,1);
f_3P_Hi_ST(1:N*T,1)=reshape(repmat(Alpha,1,T)',N*T,1);%uit
f_3P_Hi_ST(N*T+1:2*N*T,1)=reshape(repmat(Beta,1,T)',N*T,1);%pit
f_3P_Hi_ST(3*N*T+1:4*N*T,1)=ones(N*T,1);%Sit_wan
for i=1:N
    f_3P_Hi_ST((i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T,1)=sparse([0;ones(T-1,1)*Alpha(i)]+ThHot_cost_start(i)*ones(T,1));%sit
end

f_3P_H=sparse(6*N*T,1);
f_3P_H(1:N*T,1)=reshape(repmat(Alpha_wan,1,T)',N*T,1);%uit
f_3P_H(N*T+1:2*N*T,1)=reshape(repmat(Beta_wan,1,T)',N*T,1);%pit
f_3P_H(2*N*T+1:3*N*T,1)=reshape(repmat(ThHot_cost_start,1,T)',N*T,1);%sit
f_3P_H(3*N*T+1:4*N*T,1)=ones(N*T,1);%SitWan

f_3P_Pr=sparse(5*N*T,1);
f_3P_Pr(1:N*T,1)=reshape(repmat(Alpha_wan,1,T)',N*T,1);%uit
f_3P_Pr(N*T+1:2*N*T,1)=reshape(repmat(Beta_wan,1,T)',N*T,1);%pit
f_3P_Pr(2*N*T+1:3*N*T,1)=reshape(repmat(ThHot_cost_start,1,T)',N*T,1);%sit
f_3P_Pr(3*N*T+1:4*N*T,1)=ones(N*T,1);%SitWan
%不等式约束
Aineq_2P_Co=[Aineq_constraint_min_upordown_time_2P_Co;Aineq_constraint_start_cost_2P_Co;Aineq_constraints_spinning_reserve_requirment_2P_Co;Aineq_constraint_generation_limit_2P_Co;Aineq_constraints_ramp_up_2P_Co];%];%];%];
bineq_2P_Co=[bineq_constraint_min_upordown_time_2P_Co;bineq_constraint_start_cost_2P_Co;bineq_constraints_spinning_reserve_requirment_2P_Co;bineq_constraint_generation_limit_2P_Co;bineq_constraints_ramp_up_2P_Co];%];%];%];

Aineq_2P_Ti=[Aineq_constraint_min_upordown_time_2P_Ti;Aineq_constraint_start_cost_2P_Ti;Aineq_constraints_spinning_reserve_requirment_2P_Ti;Aineq_constraint_generation_limit_2P_Ti;Aineq_constraints_ramp_up_2P_Ti];%];%];%];
bineq_2P_Ti=[bineq_constraint_min_upordown_time_2P_Ti;bineq_constraint_start_cost_2P_Ti;bineq_constraints_spinning_reserve_requirment_2P_Ti;bineq_constraint_generation_limit_2P_Ti;bineq_constraints_ramp_up_2P_Ti];%];%];%];

Aineq_3P_Ti=[Aineq_constraint_min_upordown_time_3P_Ti;Aineq_constraint_start_cost_3P_Ti;Aineq_constraints_spinning_reserve_requirment_3P_Ti;Aineq_constraint_generation_limit_3P_Ti;Aineq_constraints_ramp_up_3P_Ti];%];%];%];%];
bineq_3P_Ti=[bineq_constraint_min_upordown_time_3P_Ti;bineq_constraint_start_cost_3P_Ti;bineq_constraints_spinning_reserve_requirment_3P_Ti;bineq_constraint_generation_limit_3P_Ti;bineq_constraints_ramp_up_3P_Ti];%];%];%];%];

Aineq_3P_Ti_ST=[Aineq_constraint_min_upordown_time_3P_Ti_ST;Aineq_constraint_start_cost_3P_Ti_ST;Aineq_constraints_spinning_reserve_requirment_3P_Ti_ST;Aineq_constraint_generation_limit_3P_Ti_ST;Aineq_constraints_ramp_up_3P_Ti_ST];%];%];%];%];
bineq_3P_Ti_ST=[bineq_constraint_min_upordown_time_3P_Ti_ST;bineq_constraint_start_cost_3P_Ti_ST;bineq_constraints_spinning_reserve_requirment_3P_Ti_ST;bineq_constraint_generation_limit_3P_Ti_ST;bineq_constraints_ramp_up_3P_Ti_ST];%];%];%];%];

Aineq_3P_HD=[Aineq_constraint_min_upordown_time_3P_HD;Aineq_constraint_start_cost_3P_HD;Aineq_constraints_spinning_reserve_requirment_3P_HD;Aineq_constraint_state_3P_HD;Aineq_constraint_generation_limit_3P_HD;Aineq_constraints_ramp_up_3P_HD];%];%];%];%];
bineq_3P_HD=[bineq_constraint_min_upordown_time_3P_HD;bineq_constraint_start_cost_3P_HD;bineq_constraints_spinning_reserve_requirment_3P_HD;bineq_constraint_state_3P_HD;bineq_constraint_generation_limit_3P_HD;bineq_constraints_ramp_up_3P_HD];%];%];%];%];

Aineq_3P_HD_Pr=[Aineq_constraint_min_upordown_time_3P_HD_Pr;Aineq_constraint_start_cost_3P_HD_Pr;Aineq_constraints_spinning_reserve_requirment_3P_HD_Pr;Aineq_constraint_generation_limit_3P_HD_Pr;Aineq_constraints_ramp_up_3P_HD_Pr];%];%];%];%];
bineq_3P_HD_Pr=[bineq_constraint_min_upordown_time_3P_HD_Pr;bineq_constraint_start_cost_3P_HD_Pr;bineq_constraints_spinning_reserve_requirment_3P_HD_Pr;bineq_constraint_generation_limit_3P_HD_Pr;bineq_constraints_ramp_up_3P_HD_Pr];%];%];%];%];
%等式约束
Aeq_2P_Co=[Aeq_constraint_power_balance_2P_Co;Aeq_constraints_initial_statues_2P_Co;Aeq_constraint_state_2P_Co];%];%];
beq_2P_Co=[beq_constraint_power_balance_2P_Co;beq_constraints_initial_statues_2P_Co;beq_constraint_state_2P_Co];%];%];

Aeq_2P_Ti=[Aeq_constraint_power_balance_2P_Ti;Aeq_constraints_initial_statues_2P_Ti;Aeq_constraint_state_2P_Ti];%];%];
beq_2P_Ti=[beq_constraint_power_balance_2P_Ti;beq_constraints_initial_statues_2P_Ti;beq_constraint_state_2P_Ti];%];%];

Aeq_3P_Ti=[Aeq_constraint_power_balance_3P_Ti;Aeq_constraints_initial_statues_3P_Ti;Aeq_constraint_state_3P_Ti];%];%];
beq_3P_Ti=[beq_constraint_power_balance_3P_Ti;beq_constraints_initial_statues_3P_Ti;beq_constraint_state_3P_Ti];%];%];

Aeq_3P_Ti_ST=[Aeq_constraint_power_balance_3P_Ti_ST;Aeq_constraints_initial_statues_3P_Ti_ST;Aeq_constraint_state_3P_Ti_ST];%];%];
beq_3P_Ti_ST=[beq_constraint_power_balance_3P_Ti_ST;beq_constraints_initial_statues_3P_Ti_ST;beq_constraint_state_3P_Ti_ST];%];%];

Aeq_3P_HD=[Aeq_constraint_power_balance_3P_HD;Aeq_constraints_initial_statues_3P_HD;Aeq_constraint_state_3P_HD];%];%];
beq_3P_HD=[beq_constraint_power_balance_3P_HD;beq_constraints_initial_statues_3P_HD;beq_constraint_state_3P_HD];%];%];

Aeq_3P_HD_Pr=[Aeq_constraint_power_balance_3P_HD_Pr;Aeq_constraints_initial_statues_3P_HD_Pr;Aeq_constraint_state_3P_HD_Pr];%];%];;
beq_3P_HD_Pr=[beq_constraint_power_balance_3P_HD_Pr;beq_constraints_initial_statues_3P_HD_Pr;beq_constraint_state_3P_HD_Pr];%];%];;
%变量上下界
lb=sparse(5*N*T,1);
ub=inf*ones(5*N*T,1);
ub(1:N*T,1)=ones(N*T,1);%uit
ub(N*T+1:2*N*T,1)=sparse(reshape(repmat(ThPimax,1,T)',N*T,1));%pit
ub(2*N*T+1:3*N*T,1)=ones(N*T,1);%sit
ub(4*N*T+1:5*N*T,1)=ones(N*T,1);%dit

lb_3P_HD=sparse(6*N*T,1);
ub_3P_HD=inf*ones(6*N*T,1);
ub_3P_HD(1:N*T,1)=ones(N*T,1);%uit
ub_3P_HD(N*T+1:2*N*T,1)=ones(N*T,1);%pit
ub_3P_HD(2*N*T+1:3*N*T,1)=ones(N*T,1);%sit
ub_3P_HD(4*N*T+1:5*N*T,1)=ones(N*T,1);%dit
ub_3P_HD(5*N*T+1:6*N*T,1)=ones(N*T,1);%T3

lb_3P_HD_Pr=sparse(5*N*T,1);
ub_3P_HD_Pr=inf*ones(5*N*T,1);
ub_3P_HD_Pr(1:N*T,1)=ones(N*T,1);%uit
ub_3P_HD_Pr(N*T+1:2*N*T,1)=ones(N*T,1);%pit
ub_3P_HD_Pr(2*N*T+1:3*N*T,1)=ones(N*T,1);%sit
ub_3P_HD_Pr(4*N*T+1:5*N*T,1)=ones(N*T,1);%dit

%变量类型
ctype='';
ctype(1,1:5*N*T)='C';
ctype(1,1:N*T)='B';
ctype(1,2*N*T+1:3*N*T)='B';
ctype(1,4*N*T+1:5*N*T)='B';

ctype_3P_HD='';
ctype_3P_HD(1,1:6*N*T)='C';
ctype_3P_HD(1,1:N*T)='B';
ctype_3P_HD(1,2*N*T+1:3*N*T)='B';
ctype_3P_HD(1,4*N*T+1:5*N*T)='B';
ctype_3P_HD(1,5*N*T+1:6*N*T)='B';

ctype_3P_HD_Pr='';
ctype_3P_HD_Pr(1,1:5*N*T)='C';
ctype_3P_HD_Pr(1,1:N*T)='B';
ctype_3P_HD_Pr(1,2*N*T+1:3*N*T)='B';
ctype_3P_HD_Pr(1,4*N*T+1:5*N*T)='B';
%% construct six models


if(CET==0)
    %2P_Co
    model.M_2P_Co.H=H;
    model.M_2P_Co.f=f;
    model.M_2P_Co.Aeq=Aeq_2P_Co;
    model.M_2P_Co.beq=beq_2P_Co;
    model.M_2P_Co.Aineq=Aineq_2P_Co;
    model.M_2P_Co.bineq=bineq_2P_Co;
    model.M_2P_Co.ctype=ctype;
    model.M_2P_Co.lb=lb;
    model.M_2P_Co.ub=ub;
    
    model.M_2P_Ti.H=H;
    model.M_2P_Ti.f=f;
    model.M_2P_Ti.Aeq=Aeq_2P_Ti;
    model.M_2P_Ti.beq=beq_2P_Ti;
    model.M_2P_Ti.Aineq=Aineq_2P_Ti;
    model.M_2P_Ti.bineq=bineq_2P_Ti;
    model.M_2P_Ti.ctype=ctype;
    model.M_2P_Ti.lb=lb;
    model.M_2P_Ti.ub=ub;
    
    model.M_3P_Ti.H=H;
    model.M_3P_Ti.f=f;
    model.M_3P_Ti.Aeq=Aeq_3P_Ti;
    model.M_3P_Ti.beq=beq_3P_Ti;
    model.M_3P_Ti.Aineq=Aineq_3P_Ti;
    model.M_3P_Ti.bineq=bineq_3P_Ti;
    model.M_3P_Ti.ctype=ctype;
    model.M_3P_Ti.lb=lb;
    model.M_3P_Ti.ub=ub;
    
    model.M_3P_Ti_ST.H=H;
    model.M_3P_Ti_ST.f=f_3P_Hi_ST;
    model.M_3P_Ti_ST.Aeq=Aeq_3P_Ti_ST;
    model.M_3P_Ti_ST.beq=beq_3P_Ti_ST;
    model.M_3P_Ti_ST.Aineq=Aineq_3P_Ti_ST;
    model.M_3P_Ti_ST.bineq=bineq_3P_Ti_ST;
    model.M_3P_Ti_ST.ctype=ctype;
    model.M_3P_Ti_ST.lb=lb;
    model.M_3P_Ti_ST.ub=ub;
    
    model.M_3P_HD.H=H_wan;
    model.M_3P_HD.f=f_3P_H;
    model.M_3P_HD.Aeq=Aeq_3P_HD;
    model.M_3P_HD.beq=beq_3P_HD	;
    model.M_3P_HD.Aineq=Aineq_3P_HD;
    model.M_3P_HD.bineq=bineq_3P_HD;
    model.M_3P_HD.ctype=ctype_3P_HD;
    model.M_3P_HD.lb=lb_3P_HD;
    model.M_3P_HD.ub=ub_3P_HD;
    
    model.M_3P_HD_Pr.H=H_wan_Pr;
    model.M_3P_HD_Pr.f=f_3P_Pr;
    model.M_3P_HD_Pr.Aeq=Aeq_3P_HD_Pr;
    model.M_3P_HD_Pr.beq=beq_3P_HD_Pr;
    model.M_3P_HD_Pr.Aineq=Aineq_3P_HD_Pr;
    model.M_3P_HD_Pr.bineq=bineq_3P_HD_Pr;
    model.M_3P_HD_Pr.ctype=ctype_3P_HD_Pr;
    model.M_3P_HD_Pr.lb=lb_3P_HD_Pr;
    model.M_3P_HD_Pr.ub=ub_3P_HD_Pr;
else
    model.M_2P_Co.H=[H,sparse(size(H,1),4);sparse(4,size([H,sparse(size(H,1),4)],2))];
    model.M_2P_Co.f=[f;price_buy;-1*price_sell;0;0];
    model.M_2P_Co.Aeq=[Aeq_2P_Co,sparse(size(Aeq_3P_Ti_ST,1),4)];
    model.M_2P_Co.beq=beq_2P_Co;
    model.M_2P_Co.Aineq=[Aineq_2P_Co,sparse(size(Aineq_2P_Co,1),4);Aineq_trade_constraint_3P_HD_Pr];
    model.M_2P_Co.bineq=[bineq_2P_Co;bineq_trade_constraint_3P_HD_Pr];
    model.M_2P_Co.ctype=strcat(ctype,'CCBB');
    model.M_2P_Co.lb=[lb;sparse(4,1)];
    model.M_2P_Co.ub=[ub;inf;inf;1;1];
    model.M_2P_Co.Q=Q_emission;
    model.M_2P_Co.l=l_emission;
    model.M_2P_Co.r=r_emission;
    
    model.M_2P_Ti.H=[H,sparse(size(H,1),4);sparse(4,size([H,sparse(size(H,1),4)],2))];
    model.M_2P_Ti.f=[f;price_buy;-1*price_sell;0;0];
    model.M_2P_Ti.Aeq=[Aeq_2P_Ti,sparse(size(Aeq_2P_Ti,1),4)];
    model.M_2P_Ti.beq=beq_2P_Ti;
    model.M_2P_Ti.Aineq=[Aineq_2P_Ti,sparse(size(Aineq_2P_Ti,1),4);Aineq_trade_constraint_3P_HD_Pr];
    model.M_2P_Ti.bineq=[bineq_2P_Ti;bineq_trade_constraint_3P_HD_Pr];
    model.M_2P_Ti.ctype=strcat(ctype,'CCBB');
    model.M_2P_Ti.lb=[lb;sparse(4,1)];
    model.M_2P_Ti.ub=[ub;inf;inf;1;1];
    model.M_2P_Ti.Q=Q_emission;
    model.M_2P_Ti.l=l_emission;
    model.M_2P_Ti.r=r_emission;
    
    model.M_3P_Ti.H=[H,sparse(size(H,1),4);sparse(4,size([H,sparse(size(H,1),4)],2))];
    model.M_3P_Ti.f=[f;price_buy;-1*price_sell;0;0];
    model.M_3P_Ti.Aeq=[Aeq_3P_Ti,sparse(size(Aeq_3P_Ti,1),4)];
    model.M_3P_Ti.beq=beq_3P_Ti;
    model.M_3P_Ti.Aineq=[Aineq_3P_Ti,sparse(size(Aineq_3P_Ti,1),4);Aineq_trade_constraint_3P_HD_Pr];
    model.M_3P_Ti.bineq=[bineq_3P_Ti;bineq_trade_constraint_3P_HD_Pr];
    model.M_3P_Ti.ctype=strcat(ctype,'CCBB');
    model.M_3P_Ti.lb=[lb;sparse(4,1)];
    model.M_3P_Ti.ub=[ub;inf;inf;1;1];
    model.M_3P_Ti.Q=Q_emission;
    model.M_3P_Ti.l=l_emission;
    model.M_3P_Ti.r=r_emission;
    
    model.M_3P_Ti_ST.H=[H,sparse(size(H,1),4);sparse(4,size([H,sparse(size(H,1),4)],2))];
    model.M_3P_Ti_ST.f=[f_3P_Hi_ST;price_buy;-1*price_sell;0;0];
    model.M_3P_Ti_ST.Aeq=[Aeq_3P_Ti_ST,sparse(size(Aeq_3P_Ti_ST,1),4)];
    model.M_3P_Ti_ST.beq=beq_3P_Ti_ST;
    model.M_3P_Ti_ST.Aineq=[Aineq_3P_Ti_ST,sparse(size(Aineq_3P_Ti_ST,1),4);Aineq_trade_constraint_3P_HD_Pr];
    model.M_3P_Ti_ST.bineq=[bineq_3P_Ti_ST;bineq_trade_constraint_3P_HD_Pr];
    model.M_3P_Ti_ST.ctype=strcat(ctype,'CCBB');
    model.M_3P_Ti_ST.lb=[lb;sparse(4,1)];
    model.M_3P_Ti_ST.ub=[ub;inf;inf;1;1];   
    for i=1:N
        l_emission((i-1)*T+2*N*T+1:(i-1)*T+2*N*T+T,1)=sparse([0;ones(T-1,1)*Alpha(i)]+ThHot_cost_start(i)*ones(T,1));%sit
    end
    model.M_3P_Ti_ST.Q=Q_emission;
    model.M_3P_Ti_ST.l=l_emission;
    model.M_3P_Ti_ST.r=r_emission;
    
    model.M_3P_HD.H=[H_wan,sparse(size(H_wan,1),4);sparse(4,size([H_wan,sparse(size(H_wan,1),4)],2))];
    model.M_3P_HD.f=[f_3P_H;price_buy;-1*price_sell;0;0];
    model.M_3P_HD.Aeq=[Aeq_3P_HD,sparse(size(Aeq_3P_HD,1),4)];
    model.M_3P_HD.beq=beq_3P_HD	;
    model.M_3P_HD.Aineq=[Aineq_3P_HD,sparse(size(Aineq_3P_HD,1),4);Aineq_trade_constraint_3P_HD];
    model.M_3P_HD.bineq=[bineq_3P_HD;bineq_trade_constraint_3P_HD];   
    model.M_3P_HD.ctype=strcat(ctype_3P_HD,'CCBB');
    model.M_3P_HD.lb=[lb_3P_HD;sparse(4,1)];
    model.M_3P_HD.ub=[ub_3P_HD;inf;inf;1;1];
    model.M_3P_HD.Q=Q_3P_HD;
    model.M_3P_HD.l=l_3P_HD;
    model.M_3P_HD.r=r_3P_HD;    
    
    model.M_3P_HD_Pr.H=[H_wan_Pr,sparse(size(H_wan_Pr,1),4);sparse(4,size([H_wan_Pr,sparse(size(H_wan_Pr,1),4)],2))];
    model.M_3P_HD_Pr.f=[f_3P_Pr;price_buy;-1*price_sell;0;0];
    model.M_3P_HD_Pr.Aeq=[Aeq_3P_HD_Pr,sparse(size(Aeq_3P_HD_Pr,1),4)];
    model.M_3P_HD_Pr.beq=beq_3P_HD_Pr;
    model.M_3P_HD_Pr.Aineq=[Aineq_3P_HD_Pr,sparse(size(Aineq_3P_HD_Pr,1),4);Aineq_trade_constraint_3P_HD_Pr];
    model.M_3P_HD_Pr.bineq=[bineq_3P_HD_Pr;bineq_trade_constraint_3P_HD_Pr];
    model.M_3P_HD_Pr.ctype=strcat(ctype_3P_HD_Pr,'CCBB');
    model.M_3P_HD_Pr.lb=[lb_3P_HD_Pr;sparse(4,1)];
    model.M_3P_HD_Pr.ub=[ub_3P_HD_Pr;inf;inf;1;1];
    model.M_3P_HD_Pr.Q=Q_3P_HD_Pr;
    model.M_3P_HD_Pr.l=l_3P_HD_Pr;
    model.M_3P_HD_Pr.r=r_3P_HD_Pr;
end

model.N=N;
model.T=T;
model.ThHot_cost_start=ThHot_cost_start;
model.ThPimin=ThPimin;
model.ThPimax=ThPimax;
end