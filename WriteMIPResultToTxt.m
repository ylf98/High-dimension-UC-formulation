function WriteMIPResultToTxt( fid,time,Obj,Gap,nodes,iterations,cons,nonzs,Vars_I,Vars_C,i,N )
%WRITERESULTTOTXT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if(i==1)
    fprintf(fid,'Unit_Number\t\t\t Time\t\t\t Obj\t\t\t Gap\t\t\t nodes\t\t\t iterations\t\t\t cons\t\t\t nonzs\t\t\t Vars_I\t\t\t Vars_C\n');
    fprintf(fid,'%d\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t  %f\t\t  %f\n',N,time,Obj,Gap,nodes,iterations,cons,nonzs,Vars_I,Vars_C);
else
    fprintf(fid,'%d\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t %f\t\t  %f\t\t  %f\n',N,time,Obj,Gap,nodes,iterations,cons,nonzs,Vars_I,Vars_C);
end
end

