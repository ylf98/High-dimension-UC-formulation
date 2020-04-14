clear all;

TimeLimit=3600;
L =4;
index=1;
boolean_MILP=0;%1表示目标函数线性化，0表示目标函数不线性化
GAP=5e-3;
GAP_Str='0005';%如果GAP改变，此处也要变。
data=1;%数据集，1表示10-1080，2表示10-100
access=1;%1表示求解MIP问题，0表示不求解MIP问题
access_relax=1;%1表示求解连续松弛，0表示不求解连续松弛
times_T=1;%时间拓展，1表示24，2表示2*24，...
times_Uit=1;%机组拓展，1表示1*N,2表示2*N,...
solver='CPLEX';%求解器选择 "CPLEX","CP"
CET=0;%是否求解碳排放，1表示是，0表示否


if(CET==0)
    FileFoder=fullfile('.\UC_AF\');
    if(mod(data,3)==1)
        dirOutput=dir(fullfile(FileFoder,'c*.mod'));
    else
        dirOutput=dir(fullfile(FileFoder,'*w.mod'));
    end
else
    FileFoder=fullfile('.\CET_data_1\');
    dirOutput=dir(fullfile(FileFoder,'c*.mod'));
end

compare_fileName={dirOutput.name}';
FileFold='UC_AF\';
File_num=size(compare_fileName,1);
% File_num=24;
FileRemark='';

if(boolean_MILP~=0)
    FileRemark=strcat('_',GAP_Str,'_MILP_');
else
    FileRemark=strcat('_',GAP_Str,'_MIQP_');
end

if(CET==1)
    FileRemark=strcat(FileRemark,'CET_');
end

if(strcmp(solver,'CP'))
    FileRemark=strcat(FileRemark,'CP_');
end
if(strcmp(solver,'CPLEX'))
    FileRemark=strcat(FileRemark,'CPLEX_');
end

if(access==1)
    if(CET==0)
        fid_CPLEX_MIP_2P_Co=fopen(strcat('result_',FileFold,'2P_Co_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid_CPLEX_MIP_2P_Ti=fopen(strcat('result_',FileFold,'2P_Ti_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid_CPLEX_MIP_3P_Ti=fopen(strcat('result_',FileFold,'3P_Ti_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid_CPLEX_MIP_3P_Ti_ST=fopen(strcat('result_',FileFold,'3P_Ti_ST_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid_CPLEX_MIP_3P_HD=fopen(strcat('result_',FileFold,'3P_HD_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid_CPLEX_MIP_3P_HD_Pr=fopen(strcat('result_',FileFold,'3P_HD_Pr_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        
        fid.MIP_2P_Co=fid_CPLEX_MIP_2P_Co;
        fid.MIP_2P_Ti=fid_CPLEX_MIP_2P_Ti;
        fid.MIP_3P_Ti=fid_CPLEX_MIP_3P_Ti;
        fid.MIP_3P_Ti_ST=fid_CPLEX_MIP_3P_Ti_ST;
        fid.MIP_3P_HD=fid_CPLEX_MIP_3P_HD;
        fid.MIP_3P_HD_Pr=fid_CPLEX_MIP_3P_HD_Pr;
    else
        fid_CPLEX_MIP_3P_HD_Pr=fopen(strcat('result_',FileFold,'3P_HD_Pr_',num2str(index),'_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid.MIP_3P_HD_Pr=fid_CPLEX_MIP_3P_HD_Pr;
    end
end
if(access_relax==1)
    if(boolean_MILP==1)
        FileRemark=strcat(FileRemark,'lP','_');
    else
        FileRemark=strcat(FileRemark,'QP','_');
    end
    if(CET==0)
    fid_CPLEX_Relax_2P_Co=fopen(strcat('result_',FileFold,'2P_Co_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
    fid_CPLEX_Relax_2P_Ti=fopen(strcat('result_',FileFold,'2P_Ti_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
    fid_CPLEX_Relax_3P_Ti=fopen(strcat('result_',FileFold,'3P_Ti_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
    fid_CPLEX_Relax_3P_Ti_ST=fopen(strcat('result_',FileFold,'3P_Ti_ST_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
    fid_CPLEX_Relax_3P_HD=fopen(strcat('result_',FileFold,'3P_HD_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
    fid_CPLEX_Relax_3P_HD_Pr=fopen(strcat('result_',FileFold,'3P_HD_Pr_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
    
    fid.Relax_2P_Co=fid_CPLEX_Relax_2P_Co;
    fid.Relax_2P_Ti=fid_CPLEX_Relax_2P_Ti;
    fid.Relax_3P_Ti=fid_CPLEX_Relax_3P_Ti;
    fid.Relax_3P_Ti_ST=fid_CPLEX_Relax_3P_Ti_ST;
    fid.Relax_3P_HD=fid_CPLEX_Relax_3P_HD;
    fid.Relax_3P_HD_Pr=fid_CPLEX_Relax_3P_HD_Pr;
    else
            fid_CPLEX_Relax_3P_HD_Pr=fopen(strcat('result_',FileFold,'3P_HD_Pr_Relax',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
             fid.Relax_3P_HD_Pr=fid_CPLEX_Relax_3P_HD_Pr;
    end
end

while(index<=File_num)
    FileName=compare_fileName{index};
    FileName=strcat(FileName,'');
    FileType='.mod';
    FilePath=strcat(FileFoder,FileName);
    dataUC=readdataUC(FilePath,CET);      %读取文件数据
    [dataUC]=expand( dataUC,times_T,times_Uit);
    %     dataUC=readdataUC('UC_AF/5_std.mod');
    [model]=model_produce(dataUC,CET);
    orig_model=model;
    if(boolean_MILP~=0)
        [ model.M_2P_Co	 ] = obj_linear_general_tangent( model.M_2P_Co,L,model.ThPimin,model.ThPimax,model.N,model.T,CET);
        [ model.M_2P_Ti ] = obj_linear_general_tangent( model.M_2P_Ti,L,model.ThPimin,model.ThPimax,model.N,model.T,CET);
        [ model.M_3P_HD ] = obj_linear_general_tangent_project( model.M_3P_HD,L,CET);
        [ model.M_3P_HD_Pr ] = obj_linear_general_tangent_project( model.M_3P_HD_Pr,L,CET);
        [ model.M_3P_Ti ] = obj_linear_general_tangent( model.M_3P_Ti,L,model.ThPimin,model.ThPimax,model.N,model.T,CET);
        [ model.M_3P_Ti_ST ] = obj_linear_general_tangent( model.M_3P_Ti_ST,L,model.ThPimin,model.ThPimax,model.N,model.T,CET);
    end
    
    [result]=solve(model,TimeLimit,GAP,access,access_relax,data,boolean_MILP,times_T,times_Uit,index,fid,solver,CET);
   
    index=index+1;
    if(CET==1)
        fid_CPLEX_MIP_3P_HD_Pr=fopen(strcat('result_',FileFold,'3P_HD_Pr_',num2str(index),'_',FileRemark,'expandtime_',num2str(times_T),'_expandUnit_',num2str(times_Uit),'_',num2str(data),'.mod'),'a+');
        fid.MIP_3P_HD_Pr=fid_CPLEX_MIP_3P_HD_Pr;
    end
end

function [result]=solve(model,TimeLimit,GAP,access,access_relax,data,boolean_MILP,times_T,times_Uit,index,fid,solver,CET)
if(strcmp(solver,'CPLEX'))
    success_num=0;
    if(access==1)
       
        cplex_class=Cplex('MIP');
        fid_MIP_2P_Co=fid.MIP_2P_Co;
        fid_MIP_2P_Ti=fid.MIP_2P_Ti;
        fid_MIP_3P_Ti=fid.MIP_3P_Ti;
        fid_MIP_3P_Ti_ST=fid.MIP_3P_Ti_ST;
        fid_MIP_3P_HD=fid.MIP_3P_HD;
        fid_MIP_3P_HD_Pr=fid.MIP_3P_HD_Pr;
        clc;
        if(CET==1)
            strlog=strcat('log\log_CET',num2str(data));
        else
            strlog=strcat('log\log_',num2str(data));
        end
        
        if(boolean_MILP==1)
            diary(strcat(strlog,'_MILP_expand_',num2str(model.N),'_time_',num2str(times_T),'_unit_',num2str(times_Uit),'_',num2str(index),'.txt'));
        else
            diary(strcat(strlog,'_MIQP_expand_',num2str(model.N),'_time_',num2str(times_T),'_unit_',num2str(times_Uit),'_',num2str(index),'.txt'));
        end
        
        bool_M_2P_Co=0;
        bool_M_3P_HD_Pr=0;
        bool_M_2P_Ti=0;
        bool_M_3P_Ti=0;
        bool_M_3P_Ti_ST=0;
        bool_M_3P_HD=0;
        
        while(success_num<6)
            diary on;
            try
                
                if(bool_M_2P_Co==0)
                    if(CET==1)
                        [result.M_2P_Co.MIP]=CplexSolve('MIP','2P_Co',model.M_2P_Co,GAP,TimeLimit,model.N,boolean_MILP,CET,fid_MIP_2P_Co);
                    else
                        [result.M_2P_Co.MIP]=CplexSolve('MIP','2P_Co',model.M_2P_Co,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_class);
                        WriteMIPResultToTxt( fid_MIP_2P_Co,result.M_2P_Co.MIP.time,result.M_2P_Co.MIP.Obj,result.M_2P_Co.MIP.gap,result.M_2P_Co.MIP.node,result.M_2P_Co.MIP.iteration,result.M_2P_Co.MIP.cons,result.M_2P_Co.MIP.nonzs,result.M_2P_Co.MIP.VarsI,result.M_2P_Co.MIP.VarsC,index,model.N );
                        
                    end
                    bool_M_2P_Co=1;
                end
                
                if(bool_M_3P_HD_Pr==0)
                    if(CET==1)
                        [result.M_3P_HD_Pr.MIP]=CplexSolve('MIP','3P_HD_Pr',model.M_3P_HD_Pr,GAP,TimeLimit,model.N,boolean_MILP,CET,fid_MIP_3P_HD_Pr);
                    else
                        [result.M_3P_HD_Pr.MIP]=CplexSolve('MIP','3P_HD_Pr',model.M_3P_HD_Pr,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_class);
                        WriteMIPResultToTxt( fid_MIP_3P_HD_Pr,result.M_3P_HD_Pr.MIP.time,result.M_3P_HD_Pr.MIP.Obj,result.M_3P_HD_Pr.MIP.gap,result.M_3P_HD_Pr.MIP.node,result.M_3P_HD_Pr.MIP.iteration,result.M_3P_HD_Pr.MIP.cons,result.M_3P_HD_Pr.MIP.nonzs,result.M_3P_HD_Pr.MIP.VarsI,result.M_3P_HD_Pr.MIP.VarsC,index,model.N );
                        
                    end
                    bool_M_3P_HD_Pr=1;
                end
                
                if(bool_M_2P_Ti==0)
                    if(CET==1)
                        [result.M_2P_Ti.MIP]=CplexSolve('MIP','2P_Ti',model.M_2P_Ti,GAP,TimeLimit,model.N,boolean_MILP,CET,fid_MIP_2P_Ti);
                    else
                        [result.M_2P_Ti.MIP]=CplexSolve('MIP','2P_Ti',model.M_2P_Ti,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_class);
                        WriteMIPResultToTxt( fid_MIP_2P_Ti,result.M_2P_Ti.MIP.time,result.M_2P_Ti.MIP.Obj,result.M_2P_Ti.MIP.gap,result.M_2P_Ti.MIP.node,result.M_2P_Ti.MIP.iteration,result.M_2P_Ti.MIP.cons,result.M_2P_Ti.MIP.nonzs,result.M_2P_Ti.MIP.VarsI,result.M_2P_Ti.MIP.VarsC,index,model.N );
                        
                    end
                    bool_M_2P_Ti=1;
                end
                
                if(bool_M_3P_Ti==0)
                    if(CET==1)
                        [result.M_3P_Ti.MIP]=CplexSolve('MIP','3P_Ti',model.M_3P_Ti,GAP,TimeLimit,model.N,boolean_MILP,CET,fid_MIP_3P_Ti);
                    else
                        [result.M_3P_Ti.MIP]=CplexSolve('MIP','3P_Ti',model.M_3P_Ti,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_class);
                        WriteMIPResultToTxt( fid_MIP_3P_Ti,result.M_3P_Ti.MIP.time,result.M_3P_Ti.MIP.Obj,result.M_3P_Ti.MIP.gap,result.M_3P_Ti.MIP.node,result.M_3P_Ti.MIP.iteration,result.M_3P_Ti.MIP.cons,result.M_3P_Ti.MIP.nonzs,result.M_3P_Ti.MIP.VarsI,result.M_3P_Ti.MIP.VarsC,index,model.N );
                        
                    end
                    bool_M_3P_Ti=1;
                end
                
                if(bool_M_3P_Ti_ST==0)
                    if(CET==1)
                        [result.M_3P_Ti_ST.MIP]=CplexSolve('MIP','3P_Ti_ST',model.M_3P_Ti_ST,GAP,TimeLimit,model.N,boolean_MILP,CET,fid_MIP_3P_Ti_ST);
                    else
                        [result.M_3P_Ti_ST.MIP]=CplexSolve('MIP','3P_Ti_ST',model.M_3P_Ti_ST,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_class);
                        WriteMIPResultToTxt( fid_MIP_3P_Ti_ST,result.M_3P_Ti_ST.MIP.time,result.M_3P_Ti_ST.MIP.Obj,result.M_3P_Ti_ST.MIP.gap,result.M_3P_Ti_ST.MIP.node,result.M_3P_Ti_ST.MIP.iteration,result.M_3P_Ti_ST.MIP.cons,result.M_3P_Ti_ST.MIP.nonzs,result.M_3P_Ti_ST.MIP.VarsI,result.M_3P_Ti_ST.MIP.VarsC,index,model.N );
                        
                    end
                    bool_M_3P_Ti_ST=1;
                end
                
                if(bool_M_3P_HD==0)
                    if(CET==1)
                        [result.M_3P_HD.MIP]=CplexSolve('MIP','3P_HD',model.M_3P_HD,GAP,TimeLimit,model.N,boolean_MILP,CET,fid_MIP_3P_HD);
                    else
                        [result.M_3P_HD.MIP]=CplexSolve('MIP','3P_HD',model.M_3P_HD,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_class);
                        WriteMIPResultToTxt( fid_MIP_3P_HD,result.M_3P_HD.MIP.time,result.M_3P_HD.MIP.Obj,result.M_3P_HD.MIP.gap,result.M_3P_HD.MIP.node,result.M_3P_HD.MIP.iteration,result.M_3P_HD.MIP.cons,result.M_3P_HD.MIP.nonzs,result.M_3P_HD.MIP.VarsI,result.M_3P_HD.MIP.VarsC,index,model.N );
                        
                    end
                    bool_M_3P_HD=1;
                end
            catch
                diary off;
            end
            
            success_num=bool_M_2P_Co+bool_M_3P_HD_Pr+bool_M_2P_Ti+bool_M_3P_Ti+bool_M_3P_Ti_ST+bool_M_3P_HD;
        end
        diary off;
    end
    
    if(access_relax==1)
        cplex_relax = Cplex('QP');
        bool_M_2P_Co=1;
        bool_M_3P_HD_Pr=0;
        bool_M_2P_Ti=1;
        bool_M_3P_Ti=1;
        bool_M_3P_Ti_ST=1;
        bool_M_3P_HD=1;
        
        if(CET==0)
        fid_CPLEX_Relax_2P_Co=fid.Relax_2P_Co;
        fid_CPLEX_Relax_2P_Ti=fid.Relax_2P_Ti;
        fid_CPLEX_Relax_3P_Ti=fid.Relax_3P_Ti;
        fid_CPLEX_Relax_3P_Ti_ST=fid.Relax_3P_Ti_ST;
        fid_CPLEX_Relax_3P_HD=fid.Relax_3P_HD;
        fid_CPLEX_Relax_3P_HD_Pr=fid.Relax_3P_HD_Pr;
        else
                    fid_CPLEX_Relax_3P_HD_Pr=fid.Relax_3P_HD_Pr;
        end
        
        if(bool_M_2P_Co==0)
            [result.M_2P_Co.QP]=CplexSolve('QP','',model.M_2P_Co,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_relax);
            result.M_2P_Co.QP.I_LP_u=size(find((result.M_2P_Co.QP.x(1:model.N*model.T,1)-fix(result.M_2P_Co.QP.x(1:model.N*model.T,1)))==0),1)/(model.N*model.T);
            result.M_2P_Co.QP.LP_all=size(find((result.M_2P_Co.QP.x(find(model.M_2P_Co.ctype=='B'))-fix(result.M_2P_Co.QP.x(find(model.M_2P_Co.ctype=='B'))))==0),1)/size(result.M_2P_Co.QP.x(find(model.M_2P_Co.ctype=='B')),1);
            WriteRELAXResultToTxt( fid_CPLEX_Relax_2P_Co,result.M_2P_Co.QP.relaxObj,result.M_2P_Co.QP.I_LP_u,result.M_2P_Co.QP.LP_all,index,model.N );
            bool_M_2P_Co=1;
            
        end
        if(bool_M_3P_HD_Pr==0)
            [result.M_3P_HD_Pr.QP]=CplexSolve('QP','',model.M_3P_HD_Pr,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_relax);
            result.M_3P_HD_Pr.QP.I_LP_u=size(find((result.M_3P_HD_Pr.QP.x(1:model.N*model.T,1)-fix(result.M_3P_HD_Pr.QP.x(1:model.N*model.T,1)))==0),1)/(model.N*model.T);
            result.M_3P_HD_Pr.QP.LP_all=size(find((result.M_3P_HD_Pr.QP.x(find(model.M_3P_HD_Pr.ctype=='B'))-fix(result.M_3P_HD_Pr.QP.x(find(model.M_3P_HD_Pr.ctype=='B'))))==0),1)/size(result.M_3P_HD_Pr.QP.x(find(model.M_3P_HD_Pr.ctype=='B')),1);
            WriteRELAXResultToTxt( fid_CPLEX_Relax_3P_HD_Pr,result.M_3P_HD_Pr.QP.relaxObj,result.M_3P_HD_Pr.QP.I_LP_u,result.M_3P_HD_Pr.QP.LP_all,index,model.N );
            bool_M_3P_HD_Pr=1;
        end
        if(bool_M_2P_Ti==0)
            [result.M_2P_Ti.QP]=CplexSolve('QP','',model.M_2P_Ti,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_relax);
            result.M_2P_Ti.QP.I_LP_u=size(find((result.M_2P_Ti.QP.x(1:model.N*model.T,1)-fix(result.M_2P_Ti.QP.x(1:model.N*model.T,1)))==0),1)/(model.N*model.T);
            result.M_2P_Ti.QP.LP_all=size(find((result.M_2P_Ti.QP.x(find(model.M_2P_Ti.ctype=='B'))-fix(result.M_2P_Ti.QP.x(find(model.M_2P_Ti.ctype=='B'))))==0),1)/size(result.M_2P_Ti.QP.x(find(model.M_2P_Ti.ctype=='B')),1);
            WriteRELAXResultToTxt( fid_CPLEX_Relax_2P_Ti,result.M_2P_Ti.QP.relaxObj,result.M_2P_Ti.QP.I_LP_u,result.M_2P_Ti.QP.LP_all,index,model.N );
            bool_M_2P_Ti=1;
        end
        if(bool_M_3P_Ti==0)
            [result.M_3P_Ti.QP]=CplexSolve('QP','',model.M_3P_Ti,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_relax);
            
            result.M_3P_Ti.QP.I_LP_u=size(find((result.M_3P_Ti.QP.x(1:model.N*model.T,1)-fix(result.M_3P_Ti.QP.x(1:model.N*model.T,1)))==0),1)/(model.N*model.T);
            result.M_3P_Ti.QP.LP_all=size(find((result.M_3P_Ti.QP.x(find(model.M_3P_Ti.ctype=='B'))-fix(result.M_3P_Ti.QP.x(find(model.M_3P_Ti.ctype=='B'))))==0),1)/size(result.M_3P_Ti.QP.x(find(model.M_3P_Ti.ctype=='B')),1);
            WriteRELAXResultToTxt( fid_CPLEX_Relax_3P_Ti,result.M_3P_Ti.QP.relaxObj,result.M_3P_Ti.QP.I_LP_u,result.M_3P_Ti.QP.LP_all,index,model.N );
            bool_M_3P_Ti=1;
        end
        if(bool_M_3P_Ti_ST==0)
            [result.M_3P_Ti_ST.QP]=CplexSolve('QP','',model.M_3P_Ti_ST,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_relax);
            result.M_3P_Ti_ST.QP.I_LP_u=size(find((result.M_3P_Ti_ST.QP.x(1:model.N*model.T,1)-fix(result.M_3P_Ti_ST.QP.x(1:model.N*model.T,1)))==0),1)/(model.N*model.T);
            result.M_3P_Ti_ST.QP.LP_all=size(find((result.M_3P_Ti_ST.QP.x(find(model.M_3P_Ti_ST.ctype=='B'))-fix(result.M_3P_Ti_ST.QP.x(find(model.M_3P_Ti_ST.ctype=='B'))))==0),1)/size(result.M_3P_Ti_ST.QP.x(find(model.M_3P_Ti_ST.ctype=='B')),1);
            WriteRELAXResultToTxt( fid_CPLEX_Relax_3P_Ti_ST,result.M_3P_Ti_ST.QP.relaxObj,result.M_3P_Ti_ST.QP.I_LP_u,result.M_3P_Ti_ST.QP.LP_all,index,model.N );
            bool_M_3P_Ti_ST=1;
        end
        if(bool_M_3P_HD==0)
            [result.M_3P_HD.QP]=CplexSolve('QP','',model.M_3P_HD,GAP,TimeLimit,model.N,boolean_MILP,CET,[],cplex_relax);
            result.M_3P_HD.QP.I_LP_u=size(find((result.M_3P_HD.QP.x(1:model.N*model.T,1)-fix(result.M_3P_HD.QP.x(1:model.N*model.T,1)))==0),1)/(model.N*model.T);
            result.M_3P_HD.QP.LP_all=size(find((result.M_3P_HD.QP.x(find(model.M_3P_HD.ctype=='B'))-fix(result.M_3P_HD.QP.x(find(model.M_3P_HD.ctype=='B'))))==0),1)/size(result.M_3P_HD.QP.x(find(model.M_3P_HD.ctype=='B')),1);
            WriteRELAXResultToTxt( fid_CPLEX_Relax_3P_HD,result.M_3P_HD.QP.relaxObj,result.M_3P_HD.QP.I_LP_u,result.M_3P_HD.QP.LP_all,index,model.N );
            bool_M_3P_HD=1;
        end
        
        
        
        
        
        
    end
end
if(strcmp(solver,'CP'))%未完成
    bool_M_2P_Co=1;
    bool_M_3P_HD_Pr=0;
    bool_M_2P_Ti=1;
    bool_M_3P_Ti=1;
    bool_M_3P_Ti_ST=1;
    bool_M_3P_HD=1;
    
    N=model.N;
    T=model.T;
    if(bool_M_2P_Co==0)
        fid_MIP_2P_Co=fid.MIP_2P_Co;
        [x,result,time]=function_CP(model,fid_MIP_2P_Co);
        bool_M_2P_Co=1;
    end
    
    if(bool_M_3P_HD_Pr==0)
        fid_MIP_3P_HD_Pr=fid.MIP_3P_HD_Pr;
        [ x,result,time]=function_CP(model.M_3P_HD_Pr,fid_MIP_3P_HD_Pr,N,T);
        bool_M_3P_HD_Pr=1;
    end
    
    if(bool_M_2P_Ti==0)
        fid_MIP_2P_Ti=fid.MIP_2P_Ti;
        [result.M_2P_Ti.MIP]=function_CP(model,fid_MIP_2P_Ti);
        bool_M_2P_Ti=1;
    end
    
    if(bool_M_3P_Ti==0)
        fid_MIP_3P_Ti=fid.MIP_3P_Ti;
        [result.M_3P_Ti.MIP]=function_CP(model,fid_MIP_3P_Ti);
        bool_M_3P_Ti=1;
    end
    
    if(bool_M_3P_Ti_ST==0)
        fid_MIP_3P_Ti_ST=fid.MIP_3P_Ti_ST;
        [result.M_3P_Ti_ST.MIP]=function_CP(model,fid_MIP_3P_Ti_ST);
        bool_M_3P_Ti_ST=1;
    end
    
    if(bool_M_3P_HD==0)
        fid_MIP_3P_HD=fid.MIP_3P_HD;
        [result.M_3P_HD.MIP]=function_CP(model,fid_MIP_3P_HD);
        bool_M_3P_HD=1;
    end
end
end

function [result]=CplexSolve(problem,model_type,model,GAP,TimeLimit,N,boolean_MILP,CET,fid_CET,cplex_clas)
H=model.H;
f=model.f;
Aineq=model.Aineq;
bineq=model.bineq;
Aeq=model.Aeq;
beq=model.beq;
ctype=model.ctype;
lb=model.lb;
ub=model.ub;

if(CET==1)
    Q=model.Q;
    l=model.l;
    r=model.r;
end

if(strcmp(problem,'MIP'))
    disp(N);
    disp(model_type);
    
    options = cplexoptimset;
    options.Display='off';
    %options.mip.limits.nodes=0;
    options.mip.tolerances.mipgap=GAP;
    options.timelimit=TimeLimit;
    
    if(CET==0)
        if(boolean_MILP==0)
            [x,fval,exitflag,output]=cplexmiqp(H,f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
        else
            [x,fval,exitflag,output]=cplexmilp(f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype,[],options);
        end
    end
    %     result.obj=fval;
    %     result.x=x;
    %     result.time=output.time;
    result.iteration=output.iterations;
    disp(fval);
    try
%          pause(30);
        cplex_class = Cplex('MIP');
       
        disp(cplex_class.Param.mip.tolerances.mipgap.Cur);
        cplex_class.Model.Q = H;
        cplex_class.Model.obj = f;
        cplex_class.Model.lb = lb;
        cplex_class.Model.ub = ub;
        cplex_class.Model.A = [Aineq;Aeq];
        cplex_class.Model.lhs = [-Inf.*ones(size(bineq,1),1);beq];
        cplex_class.Model.rhs = [bineq;beq];
        cplex_class.Model.ctype = ctype;
        
        cplex_class.Param.mip.tolerances.mipgap.Cur = GAP;%
        cplex_class.Param.timelimit.Cur=TimeLimit;
        cplex_class.Model.sense = 'minimize';
        %         cplex_class.DisplayFunc= 'on';
        %         cplex_class.Param.mip.limits.cutpasses.Cur=-1;
        if(CET==1)
            cplex_class.addQCs(l, Q,'L',r);
            cplex_class.InfoCallback.func = @mipex4cb;%调用callback
            cplex_class.InfoCallback.data.fid=fid_CET;
            cplex_class.InfoCallback.data.Q=Q;
            cplex_class.InfoCallback.data.r=r;
            cplex_class.InfoCallback.data.l=l;
            cplex_class.DisplayFunc= 'off'; %不显示Cplex求解过程
        else
            cplex_class.Param.mip.tolerances.mipgap.Cur = GAP;
        end
        cplex_class.solve();
        result.x = cplex_class.Solution.x;
        result.time=cplex_class.Solution.time;
        result.Obj=cplex_class.Solution.objval;
        result.gap=cplex_class.Solution.miprelgap;
        %     result.node=cplex_class.Solution.nodeint; %
        result.node=cplex_class.Solution.nodecnt;
        result.cons=size([Aineq;Aeq],1);
        result.nonzs=size(find(f~=0),1)+size(find(H~=0),1)+size(find([Aineq;Aeq]~=0),1)+find(find([bineq;beq]~=0),1);%
        result.VarsI=size(find(ctype=='B'),2);
        result.VarsC=size(find(ctype=='C'),2);
    catch
        disp('CPLEX类');
    end
end

if(strcmp(problem,'QP'))
    
    
   if(CET==0)
    try
             cplex_class = Cplex('QP');
%              pause(30);
        disp(cplex_class.DisplayFunc);
        cplex_class.Model.Q = H;
        cplex_class.Model.obj = f;
        cplex_class.Model.lb = lb;
        cplex_class.Model.ub = ub;
        cplex_class.Model.A = [Aineq;Aeq];
        cplex_class.Model.lhs = [-Inf.*ones(size(bineq,1),1);beq];
        cplex_class.Model.rhs = [bineq;beq];
        
        cplex_class.Param.timelimit.Cur=TimeLimit;
        cplex_class.DisplayFunc= 'off';
        ctype(find(ctype=='B'))='C';
        cplex_class.Param.timelimit.Cur=TimeLimit;
        
        cplex_class.Model.sense = 'minimize';
        
        
        cplex_class.Model.ctype = ctype;
        cplex_class.Param.barrier.convergetol.Cur=1e-12;
        cplex_class.solve();
        result.x = cplex_class.Solution.x;
        result.relaxObj=cplex_class.Solution.objval;
    catch
        disp('QP出错');
        
        options = cplexoptimset;
        options.Display='off';
        options.timelimit=TimeLimit;
        if(boolean_MILP==0)
            [x,fval,exitflag,output]=cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,options);
        else
            [x,fval,exitflag,output]=cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,[],options);
        end
        fprintf('\n solution status=%s\n',output.cplexstatusstring);
        result.relaxObj=fval;
        result.x=x;
        result.relaxtime=output.time;
        
    end
   else
            [ model ] = obj_linear( model,4 );
             model.ctype(find(model.ctype=='B'))='C';
            cplex_solver = Cplex('Master problem for HTC');
            cplex_solver.Model.sense = 'minimize';
            cplex_solver.Model.Q = [];
            cplex_solver.Model.obj = model.f;
            cplex_solver.Model.lb = model.lb;
            cplex_solver.Model.ub = model.ub;
            cplex_solver.Model.A = [model.Aineq;model.Aeq];
            cplex_solver.Model.lhs = [-Inf.*ones(size(model.bineq,1),1);model.beq];
            cplex_solver.Model.rhs = [model.bineq;model.beq];
            cplex_solver.Model.ctype = model.ctype;
            cplex_solver.Param.timelimit.Cur=TimeLimit;
            cplex_solver.addQCs(model.l, model.Q,'L',model.r);
%             cplex_solver.InfoCallback.func = @mipex4cb;%调用callback
%             cplex_solver.InfoCallback.data.fid=fid_CPLEX;
%             cplex_solver.InfoCallback.data.Q=model.Q;
%             cplex_solver.InfoCallback.data.r=model.r;
%             cplex_solver.InfoCallback.data.l=model.l;
%             cplex_solver.DisplayFunc= 'off'; %不显示Cplex求解过程
            cplex_solver.solve();
            if(cplex_solver.Solution.status==108)
                fprintf(fid_CPLEX,'超时');
            else
                result_cplex=cplex_solver.Solution.objval;
                gap_cplex=cplex_solver.Solution.miprelgap;
                relaxa_solution=cplex_solver.Solution.bestobjval;
                x=cplex_solver.Solution.x;
                time=cplex_solver.Solution.time;
%                 result_Obj_cplex=0.5*x(1:H_num,1)'*H*x(1:H_num,1)+f'*x(1:H_num);
%                 CET_cplex=x(1:H_num,1)'*Q*x(1:H_num,1)+l'*x(1:H_num,1)-r;
%                 fprintf(fid_CPLEX,'%d    %f    %f     %f\n',0,time,result_Obj_cplex,CET_cplex);
            end
   end
end

%Cplex.Param.barrier.convergetol
%Cplex.Param.simplex.tolerances.feasibility
end

function stop = mipex4cb(info,data)
if info.IncObj<1e73
    gap = info.MipGap * 100.0;
    info.Time
    info.BestObj
    info.IncObj
    if(~isempty(info.IncX))
        x=info.IncX;
        CET_cplex_linear=x'*data.Q*x+data.l'*x-data.r;
        fprintf(data.fid,'%f    %f    %f     %f\n',info.Time,info.BestObj,info.IncObj,CET_cplex_linear);
    end
end
stop = false;
end
