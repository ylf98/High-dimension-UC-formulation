function WriteRELAXResultToTxt( fid,Obj,I_LP_u,LP_all,i,N )
%WRITERESULTTOTXT �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if(i==1)
    fprintf(fid,'Unit_Number\t\t\t Obj\t\t\t I_LP_u\t\t\t LP_all\t\t\t \n');
    fprintf(fid,'%d\t\t %f\t\t %f\t\t %f\t\t\n',N,Obj,I_LP_u,LP_all);
else
   fprintf(fid,'%d\t\t %f\t\t %f\t\t %f\t\t\n',N,Obj,I_LP_u,LP_all);
end
end

