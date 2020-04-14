function [dataUC]=expand( dataUC,times_T,times_Uit )
%expand 扩展机组数据
%   times_T为扩展时间，24时间段的倍数，1指24时段，2指24*2时段……
%   times_Uit为扩展机组个数，1指当前数据文件中机组个数的一倍,即N*1，2指N*2……

if(times_T>1)                           %扩展时间段
%     PD_ini = [dataUC.PD;dataUC.PD(1)];
%     PD=[];
%     for i=2:(dataUC.T+1)
%             a= PD_ini(i-1);
%             b= PD_ini(i);
%             deltaPD = (b-a)/times_T;
%         for j=1:times_T
%             PD = [PD;a+(j-1)*deltaPD];
%         end
%     end
%     dataUC.PD=PD;
%     
%     spin_ini= [dataUC.spin;dataUC.spin(1)];
%     spin=[];
%     for i=2:dataUC.T+1
%             a= spin_ini(i-1);
%             b= spin_ini(i);
%             deltaspin = (b-a)/times_T;
%         for j=1:times_T
%             spin = [spin;a+(j-1)*deltaspin];
%         end
%     end
%     dataUC.spin=spin;
%     dataUC.T=dataUC.T*times_T;
dataUC.PD=repmat(dataUC.PD,times_T,1);
dataUC.spin=repmat(dataUC.spin,times_T,1);
dataUC.T=dataUC.T*times_T;
end

if(times_Uit>1)                          %扩展机组相关数据
    
    dataUC.PD=dataUC.PD*times_Uit;
    dataUC.spin=dataUC.spin*times_Uit;
    dataUC.N=dataUC.N*times_Uit;
    dataUC.p_rampdown=reshape(repmat(dataUC.p_rampdown,1,times_Uit)',dataUC.N,1);
    dataUC.p_rampup=reshape(repmat(dataUC.p_rampup,1,times_Uit)',dataUC.N,1);
    dataUC.alpha=reshape(repmat(dataUC.alpha,1,times_Uit)',dataUC.N,1);
    dataUC.beta=reshape(repmat(dataUC.beta,1,times_Uit)',dataUC.N,1);
    dataUC.gamma=reshape(repmat(dataUC.gamma,1,times_Uit)',dataUC.N,1);
    dataUC.p_low=reshape(repmat(dataUC.p_low,1,times_Uit)',dataUC.N,1);
    dataUC.p_up=reshape(repmat(dataUC.p_up,1,times_Uit)',dataUC.N,1);
    dataUC.time_on_off_ini=reshape(repmat(dataUC.time_on_off_ini,1,times_Uit)',dataUC.N,1);
    dataUC.time_min_on=reshape(repmat(dataUC.time_min_on,1,times_Uit)',dataUC.N,1);
    dataUC.time_min_off=reshape(repmat(dataUC.time_min_off,1,times_Uit)',dataUC.N,1);
    dataUC.p_initial=reshape(repmat(dataUC.p_initial,1,times_Uit)',dataUC.N,1);
    dataUC.u0=reshape(repmat(dataUC.u0,1,times_Uit)',dataUC.N,1);
    dataUC.p_startup=reshape(repmat(dataUC.p_startup,1,times_Uit)',dataUC.N,1);
    dataUC.p_shutdown=reshape(repmat(dataUC.p_shutdown,1,times_Uit)',dataUC.N,1);
    dataUC.Cold_hour=reshape(repmat(dataUC.Cold_hour,1,times_Uit)',dataUC.N,1);
    dataUC.Cold_cost=reshape(repmat(dataUC.Cold_cost,1,times_Uit)',dataUC.N,1);
    dataUC.Hot_cost=reshape(repmat(dataUC.Hot_cost,1,times_Uit)',dataUC.N,1);
    
%     dataUC.a=reshape(repmat(dataUC.a,1,times_Uit)',dataUC.N,1);
%     dataUC.b=reshape(repmat(dataUC.b,1,times_Uit)',dataUC.N,1);
%     dataUC.c=reshape(repmat(dataUC.c,1,times_Uit)',dataUC.N,1);
%     dataUC.E0=dataUC.E0*times_Uit;
%     dataUC.emmission_buy_max= dataUC.emmission_buy_max*times_Uit;
%     dataUC.emmission_sell_max= dataUC.emmission_sell_max*times_Uit;
end    

end

