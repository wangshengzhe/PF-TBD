function [z_rdb]=computRealTMul2(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z,d,b,r,Dd,Lb,fai,TE)
%获取目标二的单独量测
[NumTarget, timesensor]=size(TE);%目标个数和仿真时间长度
h_rdb=zeros(Nr,Nd,Nb,timesensor);
for k=1:timesensor
if TE(:,k)==0%没有雷达有探测信息
   for i=1:Nr % i corresponding to r
        for j=1:Nd % j corresponding to d
            for m=1:Nb % m corresponding to b
                n_I=delta_n*randn;
                n_Q=delta_n*randn;
                z_rdb_I(i,j,m,k)=n_I;
                z_rdb_Q(i,j,m,k)=n_Q;
%                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);%amplitude
                z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;% power
            end
        end
    end
else %有雷达有信息
    for i=1:Nr % i corresponding to r
        for j=1:Nd % j corresponding to d
            for m=1:Nb % m corresponding to b
                n_I=delta_n*randn;
                n_Q=delta_n*randn;
                for p=2:NumTarget
                    if TE(p,k)==1
                        %每一个目标能量值的加和 
                        h_rdb(i,j,m,k)=h_rdb(i,j,m,k)+exp(-0.5*Lr*(r(i)-z((p-1)*3+1,k))^2/Dr-0.5*Ld*(d(j)-z((p-1)*3+2,k))^2/Dd-0.5*Lb*(b(m)-z((p-1)*3+3,k))^2/Db);
                    end
                end
                z_rdb_I(i,j,m,k)=Amp*h_rdb(i,j,m,k)*cos(fai)+n_I;
                z_rdb_Q(i,j,m,k)=Amp*h_rdb(i,j,m,k)*sin(fai)+n_Q;
%                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);%amplitude
                z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;
            end
        end
    end



     
end   
     
end 
end