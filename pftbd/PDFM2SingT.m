function [pb,newaxp1,fE1,axe1]=PDFM2SingT(yuzhi,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp1,E1,z_rdb11,LR1,ST,RT,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1)

%本程序为单个跟踪粒子群跟踪单个目标


%N为粒子群个数，假设为不变的
%--------------------以下为对一个目标的跟踪与检测
%xp为上一时刻粒子群,axp为本时刻更新的粒子群
axp1=xp1;
newE1=E1;
for i=1:N
        % -------------------- update existing variable E(k+1,i)
        if E1(i)==0    
            if rand<0.1
                newE1(i)=1;
            else
                newE1(i)=0;
            end
        else
            if rand<0.1
                newE1(i)=0;
            else
                newE1(i)=1;
            end
        end
        % ---------------------- update particle xp(:,k+1,i)
        if E1(i)==1&&newE1(i)==1   % paricle maintain 
            axp1(:,i)=F_cv*xp1(:,i)+G*delta_v_p1*randn(2,1);
        elseif E1(i)==0&&newE1(i)==1 % new particle, note: we can choose a better zone for new particles 
            axp1(1,i)=xmin+(xmax-xmin)*rand;
            axp1(2,i)=vxmin+(vxmax-vxmin)*rand;
            axp1(3,i)=ymin+(ymax-ymin)*rand;
            axp1(4,i)=vymin+(vymax-vymin)*rand;
        elseif newE1(i)==0
            axp1(:,i)=0;
        end
%以上为粒子状态更新
 
 %以下为粒子权值求解
 qq1(i)=computweight(axp1(:,i),LR1,ST,RT,newE1(i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11,r,d,b,Nb,Ld);

end
        
        q_sum1=sum(qq1);
    for i=1:N
        newqq1(i)=qq1(i)/q_sum1;
    end
    
    
%     newqq1=sort(newqql);
%以下为系统重采样
    xp_temp1=axp1;
    axp_temp1=axp1;
%     diffr=axp1(:,k+1,:)-xp1(:,k+1,:);
    
    EE_temp1=newE1;
%     [newxp1 fE1]=resample(
    for i=1:N
        uu=rand;
        qq_sum1=0;
        for j=1:N
            qq_sum1=qq_sum1+newqq1(j);
            if qq_sum1>=uu
                        
                newxp1(:,i)=xp_temp1(:,j);
                newaxp1(:,i)=axp_temp1(:,j);
                fE1(i)=EE_temp1(j);
                break;
            end
        end
%         q(k+1,i)=1;
    end
    %xp为重采样完粒子群
    wwww=1;
%     figure(1)
    
%     plot(xp1(1,k+1,:),xp1(3,k+1,:));
   pb=0; 
   axe1=zeros(4,1);
    for i=1:N
        pb=pb+fE1(i);
%         Pb1(1,k+1)=Pb1(1,k+1)+E1(k+1,i);
%          xe(:,k+1)=xe(:,k+1)+E(k+1,i)*xp(:,k+1,i);
         axe1=axe1+fE1(i)*newaxp1(:,i);
    end
    if pb~=0;
%        xe(:,k+1)=xe(:,k+1)/Pb(1,k+1);
       axe1=axe1/pb;
    end
    pb=pb/N;