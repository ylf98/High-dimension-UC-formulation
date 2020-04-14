FileFoder=fullfile('.\log\');
data=2;
MILP_bool=0;
if(MILP_bool==0)
    FileRemark=strcat('log_',num2str(data),'_MIQP_expand_');
    getEliminatedRowAndCol='MIQP Presolve eliminated';
    getReduce_nonzeros_cons_Vars='Reduced MIQP has';
else
    FileRemark=strcat('log_',num2str(data),'_MILP_expand_');
    getEliminatedRowAndCol='MIP Presolve eliminated';
    getReduce_nonzeros_cons_Vars='Reduced MIP has';
end
dirOutput=dir(fullfile(FileFoder,strcat(FileRemark,'*.txt')));
compare_fileName={dirOutput.name}';
FileFold='UC_AF\';
File_num=size(compare_fileName,1);

fid_CPLEX_MIP_2P_Co=fopen(strcat('result_',FileFold,'2P_Co_',FileRemark,'.mod'),'a+');
fid_CPLEX_MIP_2P_Ti=fopen(strcat('result_',FileFold,'2P_Ti_',FileRemark,'.mod'),'a+');
fid_CPLEX_MIP_3P_Ti=fopen(strcat('result_',FileFold,'3P_Ti_',FileRemark,'.mod'),'a+');
fid_CPLEX_MIP_3P_Ti_ST=fopen(strcat('result_',FileFold,'3P_Ti_ST_',FileRemark,'.mod'),'a+');
fid_CPLEX_MIP_3P_HD=fopen(strcat('result_',FileFold,'3P_HD_',FileRemark,'.mod'),'a+');
fid_CPLEX_MIP_3P_HD_Pr=fopen(strcat('result_',FileFold,'3P_HD_Pr_',FileRemark,'.mod'),'a+');
index=1;
mark=0;
while(index<=File_num)
    FileName=compare_fileName{index};
    FileName=strcat(FileName,'');
    FileType='.txt';
    FilePath=strcat(FileFoder,FileName);
    fidin=fopen(FilePath); %打开文件
    
    num_Gomory_ractional_cuts=0;
    num_Zero_half_cuts=0;
    num_Mixed_integer_rounding_cuts=0;
    num_Flow_cuts=0;
    num_Implied_bound_cuts=0;
    num_Cover_cuts=0;
    tline=fgetl(fidin); %读取一行
    N=str2double(tline);
    while(~feof(fidin)) %judging whether the end of file
        tline=fgetl(fidin); %读取一行
        if ~isempty(tline) %判断是不是空行
            switch (tline)
                case '2P_Co'
                    mark=1;
                case '2P_Ti'
                    mark=2;
                case '3P_Ti'
                    mark=3;
                case '3P_Ti_ST'
                    mark=4;
                case '3P_HD'
                    mark=5;
                case '3P_HD_Pr'
                    mark=6;
            end
            
            if(~isempty(strfind(tline,getEliminatedRowAndCol)))
                if(MILP_bool==0)
                    str_location=strfind(tline,'MIQP Presolve');
                    row_location=strfind(tline,'rows');
                    eliminated_rows=str2double(tline(str_location+length('MIQP Presolve eliminated'):row_location-1));
                    and_location=strfind(tline,'and');
                    columns_location=strfind(tline,'columns');
                    eliminated_columns=str2double(tline(and_location+length('and '):columns_location-1));
                else
                    str_location=strfind(tline,'MIP Presolve');
                    row_location=strfind(tline,'rows');
                    eliminated_rows=str2double(tline(str_location+length('MIP Presolve eliminated'):row_location-1));
                    and_location=strfind(tline,'and');
                    columns_location=strfind(tline,'columns');
                    eliminated_columns=str2double(tline(and_location+length('and '):columns_location-1));
                end
            end
            
            if(~isempty(strfind(tline,getReduce_nonzeros_cons_Vars)))
                if(~isempty(strfind(tline,'binaries')))
                    binaries_location=strfind(tline,'binaries');
                    reduce_Vars_I=str2double(tline(length(getReduce_nonzeros_cons_Vars)+1:binaries_location-1));
                    reduce_Vars_C=reduce_Vars-reduce_Vars_I;
                else
                    row_location=strfind(tline,'rows');
                    reduce_cons=str2double(tline(length(getReduce_nonzeros_cons_Vars)+1:row_location-1));
                    columns_location=strfind(tline,'columns');
                    reduce_Vars=str2double(tline(row_location+length('rows, '):columns_location-1));
                    and_location=strfind(tline,'and');
                    nonzeros_location=strfind(tline,'nonzeros');
                    reduce_nonzeros=str2double(tline(and_location+length('and'):nonzeros_location-1));
                end
                
            end
            
            if(~isempty(strfind(tline,'   Node  Left     Objective  IInf  Best Integer    Best Bound    ItCnt     Gap')))
                fgetl(fidin); %读取一行
                data=fscanf(fidin,'%f',[1,3]); %读取一行
                while(isempty(data))
                    fgetl(fidin); %读取一行
                    data=fscanf(fidin,'%f',[1,3]); %读取一行
                end
                root_obj=data(3);
            end
            
            if(~isempty(strfind(tline,'Cover cuts applied:')))
                cuts_location=strfind(tline,':');
                num_Cover_cuts=str2double(tline(cuts_location+length(':'):end));
            end
            
            if(~isempty(strfind(tline,'Implied bound cuts applied:')))
                cuts_location=strfind(tline,':');
                num_Implied_bound_cuts=str2double(tline(cuts_location+length(':'):end));
            end
            
            if(~isempty(strfind(tline,'Flow cuts applied:')))
                cuts_location=strfind(tline,':');
                num_Flow_cuts=str2double(tline(cuts_location+length(':'):end));
            end
            
            if(~isempty(strfind(tline,'Mixed integer rounding cuts applied:')))
                cuts_location=strfind(tline,':');
                num_Mixed_integer_rounding_cuts=str2double(tline(cuts_location+length(':'):end));
            end
            
            if(~isempty(strfind(tline,'Zero-half cuts applied:')))
                cuts_location=strfind(tline,':');
                num_Zero_half_cuts=str2double(tline(cuts_location+length(':'):end));
            end
            
            if(~isempty(strfind(tline,'Gomory fractional cuts applied:')))
                cuts_location=strfind(tline,':');
                num_Gomory_ractional_cuts=str2double(tline(cuts_location+length(':'):end));
            end
            
            if(~isempty(strfind(tline,'Total (root+branch&cut)')))
                cuts=num_Gomory_ractional_cuts+num_Zero_half_cuts+num_Mixed_integer_rounding_cuts+num_Flow_cuts+num_Implied_bound_cuts+num_Cover_cuts;
                switch (mark)
                    case 1
                        writeLogtxt(fid_CPLEX_MIP_2P_Co,index,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
                    case 2
                        writeLogtxt(fid_CPLEX_MIP_2P_Ti,index,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
                    case 3
                        writeLogtxt(fid_CPLEX_MIP_3P_Ti,index,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
                    case 4
                        writeLogtxt(fid_CPLEX_MIP_3P_Ti_ST,index,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
                    case 5
                        writeLogtxt(fid_CPLEX_MIP_3P_HD,index,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
                    case 6
                        writeLogtxt(fid_CPLEX_MIP_3P_HD_Pr,index,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
                end
                num_Gomory_ractional_cuts=0;
                num_Zero_half_cuts=0;
                num_Mixed_integer_rounding_cuts=0;
                num_Flow_cuts=0;
                num_Implied_bound_cuts=0;
                num_Cover_cuts=0;
                root_obj=0;
                cuts=0;
                reduce_cons=0;
                reduce_nonzeros=0;
                reduce_Vars=0;
                reduce_Vars_I=0;
                reduce_Vars_C=0;
            end
            
        end
    end
    index=index+1;
end

function writeLogtxt(fid,i,N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C)
if(i==1)
    fprintf(fid,'Unit_Number\t\t\t Cuts\t\t\t rootObj\t\t\t reduce_cons\t\t\t reduce_nonzeros\t\t\t reduce_Vars\t\t\t reduce_Vars_I\t\t\t  reduce_Vars_C\n');
    fprintf(fid,'%d\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t\n',N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
else
    fprintf(fid,'%d\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t\n',N,cuts,root_obj,reduce_cons,reduce_nonzeros,reduce_Vars,reduce_Vars_I,reduce_Vars_C);
end
end
