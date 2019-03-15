function [ E1,axp1,axe1,Pb1]=PFsingle1111(E1,N,k,F_cv,axp1,G,delta_v_p1,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,LR1,ST,RT11,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11,r,d,b,Nb,Ld,Pb1,axe1)
% for k=1:steps_temp1-1
% [pb,newaxp1,fE1,axe]=PDFM(xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp1,E1,z_rdb11(:,:,:,k+1),LR1,ST(k+1),RT11(k+1),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
% Pb1(1,k+1)=pb;  %k+1时刻的目标出现概率     
% axe1(:,k+1)=axe;%k+1时刻的目标位置
% xp1=newaxp1;
% E1=fE1;
%end
for i=1:N
 % -------------------- update existing variable E(k+1,i)
        if E1(k,i)==0    
            if rand<0.1
                E1(k+1,i)=1;
            else
                E1(k+1,i)=0;
            end
        else
            if rand<0.1
                E1(k+1,i)=0;
            else
                E1(k+1,i)=1;
            end
        end
        % ---------------------- update particle xp(:,k+1,i)
        if E1(k,i)==1&&E1(k+1,i)==1   % paricle maintain 
            axp1(:,k+1,i)=F_cv*axp1(:,k,i)+G*delta_v_p1*randn(2,1);
        elseif E1(k,i)==0&&E1(k+1,i)==1 % new particle, note: we can choose a better zone for new particles 
            axp1(1,k+1,i)=xmin+(xmax-xmin)*rand;
            axp1(2,k+1,i)=vxmin+(vxmax-vxmin)*rand;
            axp1(3,k+1,i)=ymin+(ymax-ymin)*rand;
            axp1(4,k+1,i)=vymin+(vymax-vymin)*rand;
        elseif E1(k+1,i)==0
            axp1(:,k+1,i)=0;
        end
 xp1(:,k+1,i)=axp1(:,k+1,i);
 qq1(k+1,i)=computweight(xp1(:,k+1,i),LR1,ST(k+1),RT11(k+1),E1(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11(:,:,:,k+1),r,d,b,Nb,Ld);

    end
    

 
    
        q_sum1=sum(qq1(k+1,:));
    for i=1:N
        qq1(k+1,i)=qq1(k+1,i)/q_sum1;
    end
        

    xp_temp1=xp1;
    axp_temp1=axp1;
    diffr=axp1(:,k+1,:)-xp1(:,k+1,:);
    
    EE_temp1=E1(k+1,:);
    for i=1:N
        uu=rand;
        qq_sum1=0;
        for j=1:N
            qq_sum1=qq_sum1+qq1(k+1,j);
            if qq_sum1>=uu
               xp1(:,k+1,i)=xp_temp1(:,k+1,j);
                axp1(:,k+1,i)=axp_temp1(:,k+1,j);
                E1(k+1,i)=EE_temp1(j);
                break;
            end
        end
%         q(k+1,i)=1;
    end
    %xp为重采样完粒子群
   % wwww=1;

    
    for i=1:N
        Pb1(1,k+1)=Pb1(1,k+1)+E1(k+1,i);
%          xe(:,k+1)=xe(:,k+1)+E(k+1,i)*xp(:,k+1,i);
         axe1(:,k+1)=axe1(:,k+1)+E1(k+1,i)*axp1(:,k+1,i);
    end
    if Pb1(1,k+1)~=0;
%        xe(:,k+1)=xe(:,k+1)/Pb(1,k+1);
       axe1(:,k+1)=axe1(:,k+1)/Pb1(1,k+1);
    end
    
    Pb1(1,k+1)=Pb1(1,k+1)/N;
end