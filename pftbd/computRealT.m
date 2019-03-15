function [h_rdb, z_rdb_I, z_rdb_Q, z_rdb]=computRealT(interval,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z,d,b,r,Dd,Lb,fai,steps)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
for k=1:interval(1)-1%1:5
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
end
for k=interval(1):interval(end)%6:20
    for i=1:Nr
        for j=1:Nd
            for m=1:Nb
                n_I=delta_n*randn;
                n_Q=delta_n*randn;
                h_rdb(i,j,m,k)=exp(-0.5*Lr*(r(i)-z(1,k))^2/Dr-0.5*Ld*(d(j)-z(2,k))^2/Dd-0.5*Lb*(b(m)-z(3,k))^2/Db);
                z_rdb_I(i,j,m,k)=Amp*h_rdb(i,j,m,k)*cos(fai)+n_I;
                z_rdb_Q(i,j,m,k)=Amp*h_rdb(i,j,m,k)*sin(fai)+n_Q;
%                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);%amplitude
                z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;
            end
        end
    end
end
for k=interval(end)+1:steps
    for i=1:Nr % i corresponding to r
        for j=1:Nd % j corresponding to d
            for m=1:Nb % m corresponding to b
                n_I=delta_n*randn;
                n_Q=delta_n*randn;
                z_rdb_I(i,j,m,k)=n_I;
                z_rdb_Q(i,j,m,k)=n_Q;
%                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);
                z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;
            end
        end
    end
end


end

