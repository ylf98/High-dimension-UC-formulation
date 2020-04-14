%function Result=mian(FileName)
function model=constraints_produce(dataUC)

 dataUC=readdataUC('UC_AF/10_std.mod',0);
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
NodeNum=dataUC.N;

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];

%Zit=[];


%目标函数需要的变量--start--
Alpha_wan = Alpha + Beta .* ThPimin + Gama .* ThPimin .* ThPimin;           %Alpha弯弯--N*1矩阵
Beta_wan = (ThPimax - ThPimin) .* (Beta + 2 * Gama .* ThPimin);             %Beta弯弯--N*1矩阵
Gama_wan = Gama .* (ThPimax - ThPimin) .* (ThPimax - ThPimin);              %Gama弯弯--N*1矩阵

% 目标函数需要的变量--end--

                                   % 电机组爬坡约束需要的变量--start--
Piup_wan = Piup ./ (ThPimax - ThPimin);                                     %Piup弯弯--N*1矩阵
Pidown_wan = Pidown ./ (ThPimax - ThPimin);                                 %Pidown弯弯--N*1矩阵
Pistartup_wan = (Pistartup - ThPimin) ./ (ThPimax - ThPimin);               %Pistartup弯弯--N*1矩阵
Pishutdown_wan = (Pishutdown - ThPimin) ./ (ThPimax - ThPimin);             %Pishutdown弯弯--N*1矩阵
                                   % 电机组爬坡约束需要的变量--end--
Ui0 = full(spones(Thon_off_init));                                   %火电机组机组初始开/停机状态
ThPit_wani0 = (Thon_off_init - Ui0 .* ThPimin) ./ (ThPimax - ThPimin);      %火电Pit弯弯初始值
                                  %电机组最小开/停机约束需要的变量--start--
Ui = max(0,min(ones(N,1) * T,Ui0 .* (ThTime_on_min - ThTime_on_off_init)));                    %--N*1矩阵
Li = max(0,min(ones(N,1) * T,(ones(N,1) - Ui0) .*  (ThTime_off_min + ThTime_on_off_init)));     %--N*1矩阵

                            
                                              %--约束条件---
%等式约束

%电机组初始状态
Uit_state=[];
Aeq_UC_init=[];
beq_UC_init=[];

for i=1:N
    if(Ui(i) + Li(i) >= 1)
        for t=1:Ui(i)+Li(i)
            Uit_state=sparse(N,T);
            Uit_state(i,t)=1;
        
            Uit=reshape(Uit_state',1,N*T);
            PitWan=sparse(1,N*T);
            Sit=sparse(1,N*T);
            SitWan=sparse(1,N*T);
        
            Aeq_UC_init=[Aeq_UC_init;Uit,PitWan,Sit,SitWan];
            beq_UC_init=[beq_UC_init;Ui0(i)];
        end
    end
end
Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


%功率平衡约束
Aeq_Power_balance=[];
for i=1:N
   Uit=[Uit,eye(T)*ThPimin(i)];
end
for i=1:N
    PitWan=[PitWan,eye(T)*(ThPimax(i) - ThPimin(i))];
end
Sit=sparse(T,N*T);
SitWan=sparse(T,N*T);

Aeq_Power_balance=[Uit,PitWan,Sit,SitWan];
beq_Power_balance=Dt;

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];



%直流潮流
% Node_mark=sparse(1,NodeNum);
% for i=1:NodeNum
%     if(all(ismember(i,dataUC.busLoad.bus_PDQR))&&all(ismember(i,Bus_G)))
%         Node_mark(i)=3;
%     end
%     if(~all(ismember(i,dataUC.busLoad.bus_PDQR))&&all(ismember(i,Bus_G)))
%         Node_mark(i)=2;
%     end
%     if(all(ismember(i,dataUC.busLoad.bus_PDQR))&&~all(ismember(i,Bus_G)))
%         Node_mark(i)=1;
%     end
%     if(~all(ismember(i,dataUC.busLoad.bus_PDQR))&&~all(ismember(i,Bus_G)))
%         Node_mark(i)=0;
%     end
% end
% 
% 
Aeq_DCPowerFlow=[];
beq_DCPowerFlow=[];
% for i=1:NodeNum
%     for j=1:T
%     if( Node_mark(i)==3)
%         Uit_state=sparse(N,T);
%         Uit_state(find(i==Bus_G),j)=ThPimin(find(i==Bus_G))*1;
%         PitWan_state=sparse(N,T);
%         PitWan_state(find(i==Bus_G),j)=(ThPimax(find(i==Bus_G)) - ThPimin(find(i==Bus_G)))*1;
%         seita_state=sparse(NodeNum,T);
%        seita_state(:,j)=-1*ones(1,NodeNum).*B(i,:);
%        
%         
%         
%         Uit=reshape(Uit_state',1,N*T);
%         PitWan=reshape(PitWan_state',1,N*T);
%         Sit=sparse(1,N*T);
%         SitWan=sparse(1,N*T);
%         Seita=reshape((seita_state)',1,NodeNum*T);
%         
%         Aeq_DCPowerFlow=[Aeq_DCPowerFlow;Uit,PitWan,Sit,SitWan,Seita];
%         beq_DCPowerFlow=[beq_DCPowerFlow;bus_P(j,find(i==dataUC.busLoad.bus_PDQR))];
%     end
%     if(Node_mark(i)==2)
%        Uit_state=sparse(N,T);
%         Uit_state(find(i==Bus_G),j)=ThPimin(find(i==Bus_G))*1;
%         PitWan_state=sparse(N,T);
%         PitWan_state(find(i==Bus_G),j)=(ThPimax(find(i==Bus_G)) - ThPimin(find(i==Bus_G)))*1;
%         seita_state=sparse(NodeNum,T);
%         seita_state(:,j)=-1*ones(1,NodeNum).*B(i,:);
%         
%         
%         Uit=reshape(Uit_state',1,N*T);
%         PitWan=reshape(PitWan_state',1,N*T);
%         Sit=sparse(1,N*T);
%         SitWan=sparse(1,N*T);
%         Seita=reshape((seita_state)',1,NodeNum*T);
%         
%         Aeq_DCPowerFlow=[Aeq_DCPowerFlow;Uit,PitWan,Sit,SitWan,Seita];
%         beq_DCPowerFlow=[beq_DCPowerFlow;0];
%     end
%     if(Node_mark(i)==1)
%         
%         seita_state=sparse(NodeNum,T);
%        seita_state(:,j)=-1*ones(1,NodeNum).*B(i,:);
% 
%         
%         
%         Uit=sparse(1,N*T);
%         PitWan=sparse(1,N*T);
%         Sit=sparse(1,N*T);
%         SitWan=sparse(1,N*T);
%         Seita=reshape((seita_state)',1,NodeNum*T);
%         
%         Aeq_DCPowerFlow=[Aeq_DCPowerFlow;Uit,PitWan,Sit,SitWan,Seita];
%         beq_DCPowerFlow=[beq_DCPowerFlow;bus_P(j,find(i==dataUC.busLoad.bus_PDQR))];
%     end
%     if(Node_mark(i)==0)
%         seita_state=sparse(NodeNum,T);
%         seita_state(:,j)=-1*ones(1,NodeNum).*B(i,:);
%         
%         Uit=sparse(1,N*T);
%         PitWan=sparse(1,N*T);
%         Sit=sparse(1,N*T);
%         SitWan=sparse(1,N*T);
%         Seita=reshape((seita_state)',1,NodeNum*T);
%         
%         Aeq_DCPowerFlow=[Aeq_DCPowerFlow;Uit,PitWan,Sit,SitWan,Seita];
%         beq_DCPowerFlow=[beq_DCPowerFlow;0];
%     end
%     end
% end

%不等式约束

%启动状态约束 有问题
Aineq_state_init=[];

for i=1:N
    Uit_state=sparse(N,T);
    Sit_state=sparse(N,T);
    
    Uit_state(i,1)=1;
    Sit_state(i,1)=-1;
    
    Uit=reshape(Uit_state',1,N*T);
    PitWan=sparse(1,N*T);
    Sit=reshape(Sit_state',1,N*T);
    SitWan=sparse(1,N*T);

    
    Aineq_state_init=[Aineq_state_init;Uit,PitWan,Sit,SitWan];
end

bineq_state_init=Ui0;

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


Aineq_state=Aineq_state_init;
for i=1:N
    for t=2:T
        Uit_state=sparse(N,T);
        Sit_state=sparse(N,T);
        
        Uit_state(i,t-1)=-1;
        Uit_state(i,t)=1;
        Sit_state(i,t)=-1;
        
        Uit=reshape(Uit_state',1,N*T);
        PitWan=sparse(1,N*T);
        Sit=reshape(Sit_state',1,N*T);
        SitWan=sparse(1,N*T);

        Aineq_state=[Aineq_state;Uit,PitWan,Sit,SitWan];
    end
end
bineq_state=[bineq_state_init;sparse(N*(T-1),1)];

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


%Sit弯弯约束
 bineq_startup_cost=[];
 Aineq_startup_cost=[];
for i = 1:N
    for t = 1:T
            Uit_state=sparse(N,T);
            Sit_state=sparse(N,T);
            SitWan_state=sparse(N,T);
            
            Uit_state(i,max(t - ThTime_off_min(i) - ThCold_time_start(i) - 1,1):t - 1)=1;
            Sit_state(i,t)=1;
            SitWan_state(i,t)=-1;
            
            Uit=-1*(ThCold_cost_start(i) - ThHot_cost_start(i))*reshape(Uit_state',1,N*T);
            PitWan=sparse(1,N*T);
            Sit=(ThCold_cost_start(i) - ThHot_cost_start(i))*reshape(Sit_state',1,N*T);
            SitWan=reshape(SitWan_state',1,N*T);

            
            Aineq_startup_cost=[Aineq_startup_cost;Uit,PitWan,Sit,SitWan];
            
        if t - ThTime_off_min(i) - ThCold_time_start(i) - 1 <= 0 && max(0,-ThTime_on_off_init(i)) < ...
            abs(t - ThTime_off_min(i) - ThCold_time_start(i) - 1) + 1 

            bineq_startup_cost=[bineq_startup_cost;(ThCold_cost_start(i) - ThHot_cost_start(i))*1];
        else
            bineq_startup_cost=[bineq_startup_cost;(ThCold_cost_start(i) - ThHot_cost_start(i))*0];
            
        end
    end
end

for i = 1:N
    for t = 1:T
        SitWan_state=sparse(N,T);
        
        SitWan_state(i,t)=-1;
        
        Uit=sparse(1,N*T);
        PitWan=sparse(1,N*T);
        Sit=sparse(1,N*T);
        SitWan=reshape(SitWan_state',1,N*T);

        
        Aineq_startup_cost=[Aineq_startup_cost;Uit,PitWan,Sit,SitWan];
        bineq_startup_cost=[bineq_startup_cost;0];
    end
end

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


%电机组运行约束
% 0 <= ThPit_wan
% Aineq_Pwan=[zeros(N*T,N*T),-1*eye(N*T),zeros(N*T,2*N*T)];
% bineq_Pwan=zeros(N*T,1);
%ThPit_wan <= Uit
Aineq_Pwan=[];
bineq_Pwan=[];

Uit=-1*eye(N*T);
PitWan=eye(N*T);
Sit=sparse(N*T,N*T);
SitWan=sparse(N*T,N*T);


Aineq_Pwan=[Aineq_Pwan;Uit,PitWan,Sit,SitWan];
bineq_Pwan=[bineq_Pwan;sparse(N*T,1)];

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


%直流潮流
Aineq_DCPowerFlow=[];
bineq_DCPowerFlow=[];
% for t=1:T
%     for i=1:numel(branchX)
%         Seita_state=sparse(NodeNum,T);
%         Seita_state(branchI(i),t)=1*branchX(i);
%         Seita_state(branchJ(i),t)=-1*branchX(i);
%         
%         Uit=sparse(1,N*T);
%         PitWan=sparse(1,N*T);
%         Sit=sparse(1,N*T);
%         SitWan=sparse(1,N*T);
%         Seita=reshape(Seita_state',1,NodeNum*T);
%         
%         Aineq_DCPowerFlow=[Aineq_DCPowerFlow;Uit,PitWan,Sit,SitWan,Seita];
%     end
% end
% bineq_DCPowerFlow=[bineq_DCPowerFlow;reshape(repmat(branchP,1,T)',length(branchP)*T,1)];
% 
% for t=1:T
%     for i=1:numel(branchX)
%         Seita_state=sparse(NodeNum,T);
%         Seita_state(branchI(i),t)=-1*branchX(i);
%         Seita_state(branchJ(i),t)=1*branchX(i);
%         
%         Uit=sparse(1,N*T);
%         PitWan=sparse(1,N*T);
%         Sit=sparse(1,N*T);
%         SitWan=sparse(1,N*T);
%         Seita=reshape(Seita_state',1,NodeNum*T);
%         
%         Aineq_DCPowerFlow=[Aineq_DCPowerFlow;Uit,PitWan,Sit,SitWan,Seita];
%     end
% end
% bineq_DCPowerFlow=[bineq_DCPowerFlow;reshape(repmat(branchP,1,T)',length(branchP)*T,1)];

%电机组上坡约束
%ThPit_wan(:,1) - ThPit_wani0 <= Uit(:,1) .* Piup_wan + Sit(:,1) .* (Pistartup_wan - Piup_wan)
Aineq_ramp_up=[];
for i=1:N
    Ui1_state=sparse(N,T);
    Si1_state=sparse(N,T);
    Pi1Wan=sparse(N,T);
    
    Ui1_state(i,1)=-1;
    Si1_state(i,1)=-1;
    Pi1Wan(i,1)=1;
    
    Uit=Piup_wan(i)*reshape(Ui1_state',1,N*T);
    PitWan=reshape(Pi1Wan',1,N*T);
    Sit=(Pistartup_wan(i) - Piup_wan(i))*reshape(Si1_state',1,N*T);
    SitWan=sparse(1,N*T);

    
    Aineq_ramp_up=[Aineq_ramp_up;Uit,PitWan,Sit,SitWan];
end
bineq_ramp_up=ThPit_wani0;


%ThPit_wan(:,2:T) - ThPit_wan(:,1:T-1) <= Uit(:,2:T) .* (Piup_wan * ones(1,T-1)) + Sit(:,2:T) .* ((Pistartup_wan - Piup_wan) * ones(1,T-1))
Uit_state=[];
Uit=[];
for i=1:N
    for t=2:T
        Uit_state=sparse(N,T);
        Uit_state(i,t)=-1;
        Uit=[Uit;Piup_wan(i)*reshape(Uit_state',1,N*T)];
        
    end
end

PitWan=[];
PitWan_part=[];
for i=1:N
    for t=2:T
        PitWan_part=sparse(N,T);
        PitWan_part(i,t-1)=-1;
        PitWan_part(i,t)=1;
        PitWan=[PitWan;reshape(PitWan_part',1,N*T)];
    end
end

Sit_state=[];
Sit=[];
for i=1:N
    for t=2:T
        Sit_state=sparse(N,T);
        Sit_state(i,t)=-1;
        Sit=[Sit;(Pistartup_wan(i) - Piup_wan(i))*reshape(Sit_state',1,N*T)];
    end
end
SitWan=sparse(N*(T-1),N*T);


Aineq_ramp_up=[Aineq_ramp_up;Uit,PitWan,Sit,SitWan];
bineq_ramp_up=[bineq_ramp_up;sparse(N*(T-1),1)];

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


%电机组下坡约束
%ThPit_wani0 - ThPit_wan(:,1) <= Ui0 .* Pishutdown_wan + (Sit(:,1) - Uit(:,1)) .* 
Aineq_ramp_down=[];
for i=1:N
    Ui1_state=sparse(N,T);
    Si1_state=sparse(N,T);
    Pi1Wan=sparse(N,T);
    
    Ui1_state(i,1)=1;
    Si1_state(i,1)=-1;
    Pi1Wan(i,1)=-1;
    
    

    Aineq_ramp_down=[Aineq_ramp_down;(Pishutdown_wan(i) - Pidown_wan(i))*reshape(Ui1_state',1,N*T),reshape(Pi1Wan',1,N*T),...
        (Pishutdown_wan(i) - Pidown_wan(i))*reshape(Si1_state',1,N*T),sparse(1,N*T)];
end
bineq_ramp_down=Pishutdown_wan.*Ui0-ThPit_wani0;
%ThPit_wan(:,1:T-1) - ThPit_wan(:,2:T) <= Uit(:,1:T-1) .* (Pishutdown_wan * ones(1,T-1)) + (Sit(:,2:T) - Uit(:,2:T)) .* ((Pishutdown_wan - Pidown_wan) * ones(1,T-1))
Uit_state=[];
Uit_part=[];
for i=1:N
    for t=2:T
        Uit_state=sparse(N,T);
        Uit_state(i,t-1)=Pishutdown_wan(i)*-1;
        Uit_state(i,t)=(Pishutdown_wan(i) - Pidown_wan(i))*1;
        Uit_part=[Uit_part;reshape(Uit_state',1,N*T)];
        
    end
end

PitWan=[];
PitWan_part=[];
for i=1:N
    for t=2:T
        PitWan=sparse(N,T);
        PitWan(i,t-1)=1;
        PitWan(i,t)=-1;
        PitWan_part=[PitWan_part;reshape(PitWan',1,N*T)];
    end
end

Sit_state=[];
Sit_part=[];
for i=1:N
    for t=2:T
        Sit_state=sparse(N,T);
        Sit_state(i,t)=(Pishutdown_wan(i) - Pidown_wan(i))*-1;
        Sit_part=[Sit_part;reshape(Sit_state',1,N*T)];
    end
end

 Aineq_ramp_down=[Aineq_ramp_down;Uit_part,PitWan_part,...
     Sit_part,sparse(N*(T-1),N*T)];
 
 
bineq_ramp_down=[bineq_ramp_down;sparse(N*(T-1),1)];

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


%电机组最小开机/停机时间 
Aineq_upORdowm=[];
bineq_upORdowm=[];
for i = 1:N
    for t = Ui(i) + 1:T
       Uit_state=sparse(N,T);
       Sit_state=sparse(N,T);
       
       Uit_state(i,t)=-1;
       Sit_state(i,max(0,t - ThTime_on_min(i)) + 1:t)=1;
       
        Uit=reshape(Uit_state',1,N*T);
        PitWan=sparse(1,N*T);
        Sit=reshape(Sit_state',1,N*T);
        SitWan=sparse(1,N*T);


       Aineq_upORdowm=[Aineq_upORdowm;Uit,PitWan,Sit,SitWan];
       bineq_upORdowm=[bineq_upORdowm;0];
    end
end

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];

for i = 1:N
    for t = Li(i) + 1:T
        if t - ThTime_off_min(i) <= 0
           Sit_state=sparse(N,T);
           Sit_state(i,max(0,t - ThTime_off_min(i)) + 1:t)=1;
           
           Uit=sparse(1,N*T);
           PitWan=sparse(1,N*T);
           Sit=reshape(Sit_state',1,N*T);
           SitWan=sparse(1,N*T);

            
           Aineq_upORdowm=[Aineq_upORdowm;Uit,PitWan,Sit,SitWan];
           bineq_upORdowm=[bineq_upORdowm;1 - Ui0(i)];
        else
           Uit_state=sparse(N,T);
           Sit_state=sparse(N,T);
           Uit_state(i,max(0,t - ThTime_off_min(i)))=1;
           Sit_state(i,max(0,t - ThTime_off_min(i)) + 1:t)=1;
           
           Uit=reshape(Uit_state',1,N*T);
           PitWan=sparse(1,N*T);
           Sit=reshape(Sit_state',1,N*T);
           SitWan=sparse(1,N*T);

           Aineq_upORdowm=[Aineq_upORdowm;Uit,PitWan,Sit,SitWan];
           bineq_upORdowm=[bineq_upORdowm;1];
        end
    end
end

Uit=[];
PitWan=[];
Sit=[];
SitWan=[];


% 旋转热备用约束
Uit_state=[];
for i=1:N
   Uit_state=[Uit_state,-1*eye(T)*ThPimax(i)];
end

Uit=Uit_state;
PitWan=sparse(T,N*T);
Sit=sparse(T,N*T);
SitWan=sparse(T,N*T);


Aineq_spinning_reserve_require=[Uit,PitWan,Sit,SitWan];
bineq_spinning_reserve_require=-1*(Dt + Spin);



%目标函数约束--start--
%火电Pit弯弯约束
% Aineq_SitWan=[zeros(N*T,N*T),zeros(N*T,N*T),-1*eye(N*T),zeros(N*T,N*T)];
% bineq_SitWan=[zeros(N*T,1)];


%MILP
% Aineq_MILP_approximate=[];
% bineq_MILP_approximate=[];
% for l=0:L
%     a=l/L;
%     Uit=repmat(Alpha_wan-a^2*Gama_wan,T,1).*eye(N*T);
%     Sit=zeros(N*T,N*T);
%     SitWan=zeros(N*T,N*T);
%     PitWan=repmat(2*a*Gama_wan+Beta_wan,T,1).*eye(N*T);
%     Zit=eye(N*T)*-1;
%     Aineq_MILP_approximate=[Aineq_MILP_approximate;Uit_state,Sit_state,Sitwan_state,Pitwan_state,Zit_state];
%     bineq_MILP_approximate=[bineq_MILP_approximate;zeros(N*T,1)]
% end

%目标函数约束--end--

%目标函数--start--

Pwan_coefficient=[sparse(N*T,N*T),2*diag(reshape(repmat(Gama_wan,1,T)',N*T,1)),sparse(N*T,N*T),sparse(N*T,N*T)];
H=[sparse(N*T,(4*N*T));Pwan_coefficient;sparse(N*T,(4*N*T));sparse(N*T,(4*N*T))];
f=[reshape(repmat(Alpha_wan,1,T)',N*T,1);reshape(repmat(Beta_wan,1,T)',N*T,1);reshape(repmat(ThHot_cost_start,1,T)',N*T,1);ones(N*T,1)];
Aineq=[Aineq_state;Aineq_startup_cost;Aineq_upORdowm;Aineq_spinning_reserve_require;Aineq_DCPowerFlow;Aineq_Pwan;Aineq_ramp_up;Aineq_ramp_down];%
bineq=[bineq_state;bineq_startup_cost;bineq_upORdowm;bineq_spinning_reserve_require;bineq_DCPowerFlow;bineq_Pwan;bineq_ramp_up;bineq_ramp_down];%;
Aeq=[Aeq_UC_init;Aeq_Power_balance;Aeq_DCPowerFlow];
beq=[beq_UC_init;beq_Power_balance;beq_DCPowerFlow]; 
lb=[sparse(N*T,1);sparse(N*T,1);sparse(N*T,1);sparse(N*T,1)];
ub=[ones(N*T,1);ones(N*T,1);ones(N*T,1);inf*ones(N*T,1)];
ctype=[];
for i=1:(4*N*T)
    if(i<=N*T)
        ctype=strcat(ctype,'B');
    end
    if(i>N*T&&i<=2*N*T)
        ctype=strcat(ctype,'C');
    end
    if(i>2*N*T&&i<=3*N*T)
        ctype=strcat(ctype,'B');
    end
    if(i>3*N*T&&i<=4*N*T)
        ctype=strcat(ctype,'C');
    end
    if(i>4*N*T&&i<=4*N*T+NodeNum*T)
        ctype=strcat(ctype,'C');
    end
end


model.H=H;
model.f=f;
model.Aineq=Aineq;
model.bineq=bineq;
model.Aeq=Aeq;
model.beq=beq;
model.lb=lb;
model.ub=ub;
model.ctype=ctype;

%调用接口cplexmiqp
options = cplexoptimset;
options.Display='on';
% options.mip.tolerances.mipgap=5e-3;
[x,fval,exitflag,output]=cplexmiqp(H,f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
fprintf('\n solution status=%s\n',output.cplexstatusstring);
fprintf('Solution value = %f \n', fval);
end





