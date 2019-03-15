function q=computweight(xp,LR1,ST,RT,E,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb,r,d,b,Nb,Ld)
%先把粒子坐标矫正到雷达自己的坐标系，再计算权重
 q=1;
 xp(1)=xp(1)-LR1(1);
 xp(3)=xp(3)-LR1(2);
 %时间推移
xp(1)=xp(1)+xp(2)*(ST-RT);
xp(3)=xp(3)+xp(4)*(ST-RT);
% ------------------------ update the weight q(k+1,i)
        if E==0
            q=1;
        else
            zp(1)=sqrt((xp(1))^2+(xp(3))^2);
            zp(2)=(xp(1)*xp(2)+xp(3)*xp(4))/zp(1);
            zp(3)=atan(xp(1)/xp(3));
            if zp(1)<r(1)||zp(1)>r(end)||zp(2)<d(1)||zp(2)>d(end)||zp(3)<b(1)||zp(3)>b(end) % ensure the particle is valid
                q=1; 
            else
                p=3; % local likelihood parameter, very important,
                ik=ceil((zp(1)-r(1))/Dr); % the ceil bin of the target signal, represent Nr
                jk=ceil((zp(2)-d(1))/Dd);
                mk=ceil((zp(3)-b(1))/Db);
                
                for ii=max(1,ik-p):min(Nr,ik+p)
                    for jj=max(1,jk-p):min(Nd,jk+p)
                        for mm=max(1,mk-p):min(Nb,mk+p)
%                             h_rdb_p(k+1,i)=Amp*exp(-0.5*Lr*(r(ii)-zp(1,k+1,i))^2/Dr-0.5*Ld*(d(jj)-zp(2,k+1,i))^2/Dd-0.5*Lb*(b(mm)-zp(3,k+1,i))^2/Db); % note! in here, fai=0, h_rdb_p(k+1,i) is abs
%                             z_rdb_p(k+1,i)=h_rdb_p(k+1,i)*z_rdb(ii,jj,mm,k);% z_rdb_p(ii,jj,mm,k) take abs 
%                             q(k+1,i)=q(k+1,i)*exp(-h_rdb_p(k+1,i)^2/2/delta_n^2)*besseli(0,z_rdb_p(k+1,i)/delta_n^2);% note maybe besellk is correct, I do not know
                              h_rdb_p=Amp^2*exp(-Lr*(r(ii)-zp(1))^2/Dr-Ld*(d(jj)-zp(2))^2/Dd-Lb*(b(mm)-zp(3))^2/Db)+2*delta_n^2; 
                              q=q*(2*delta_n^2/h_rdb_p)*exp((0.5/delta_n^2-1/h_rdb_p)*z_rdb(ii,jj,mm));
                        end
                    end
                end
            end
        end
end

