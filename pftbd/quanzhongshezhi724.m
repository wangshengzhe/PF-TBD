%两个目标
%粒子群发现目标一后，给新的粒子群提供总量测，并设置目标一（估计值）周围的粒子权重为1，
%观察粒子群跟踪目标二的跟踪效果
%对比量测相减法
close all;
%   -------------simulate the trajectory of target----------------
steps=30; % measurement length
T=1;% scan period T=1s
MC=1;
for round=1:MC
    round
x=zeros(4,steps); % true state variable x=[x,vx,y,vy];
z=zeros(3,steps); % measurement variable z=[range,doppler,bearing];
time1_begin=5; % first turn start
time1_end=20; % first turn end
% angle1=-90; % first turn angle 90angle=90*pi/180=1.57radian
% time2_begin=50; % second turn start
% time2_end=55; % second turn end
% angle2=180; % second turn angle 180angle=180*pi/180=3.14radian
LR1=[0 -100];
DRbegin11=0;
DRbegin12=0;
LR2=LR1;
DRbegin2=0.1;
LR3=[0,-50];
DRbegin3=0.2;
LR4=LR3;
DRbegin4=0.1;
LR5=[0,0];
DRbegin5=0.2;
LR6=[0,0];
DRbegin6=0.1;
LR7=[0,50];
DRbegin7=0.1;
LR9=[0,-150];
DRbegin9=0.1;
F_cv=[1 T 0 0; % transition matrix of cv model
    0 1 0 0;
    0 0 1 T;
    0 0 0 1];

G=[T^2/2 0; % process noise transition matrix, assuming the noise distribution is the same for x and y axis
    T 0
    0 T^2/2
    0 T];
delta_v=0.01; %process noise standard deviation here, equal to 0;
Qv=G*delta_v^2*G';% process noise covariance

t_appear=[6 6];% target appearing time
t_disappear=[20 20]; % target disappearing time
%interval=t_appear:t_disappear;% exist time of target
Oinix1=[300,0.3,-30,0.2]; % initial state %统一坐标系下的目标一初始位置
Oinix2=[300,0.3,-80,0.2]; % initial state %统一坐标系下的目标二初始位置
%有问题
[RT11 ,intervalN11, lengthRR11, detT11 ,x11, ST]=Predata(steps,Oinix1,t_appear(1),t_disappear(1),T,DRbegin11,LR1,F_cv,G,delta_v);
[RT12 ,intervalN12, lengthRR12, detT11 ,x12, ST]=Predata(steps,Oinix2,t_appear(2),t_disappear(2),T,DRbegin12,LR1,F_cv,G,delta_v);



ax1=x11;
ax1(1,:)=x11(1,:)+LR1(1);%恢复坐标，x11雷达自身坐标系的轨迹。ax1统一坐标系下的目标轨迹
ax1(3,:)=x11(3,:)+LR1(2);

ax2=x12;
ax2(1,:)=x12(1,:)+LR1(1);
ax2(3,:)=x12(3,:)+LR1(2);


z11=z;
compuz11=computT(x11,intervalN11);%在雷达自己的坐标系下求三维的量测
z11(:,1:intervalN11(lengthRR11(2)))=compuz11;%z11(:,1:21)=compuz11

z12=z;
compuz12=computT(x12,intervalN12);
z12(:,1:intervalN12(lengthRR12(2)))=compuz12;%z12(:,1:21)=compuz12

z1=[z11;z12];%雷达自己坐标系下的量测

TE=zeros(2,steps);
TE(1,intervalN11(1):intervalN11(lengthRR11(2)))=1;
TE(2,intervalN12(1):intervalN12(lengthRR12(2)))=1;





interval11=intervalN11;
interval12=intervalN12;
  

% -------------  simulate amplitude data ------------
%计算真的量测

Amp=sqrt(0.3); % signal amplitude of target; 4(3DB);8(6DB);16(9DB);32(12DB),assuming the amplitude is a constant
delta_n=1; % noise standard deviation N=N_real+j*N_imaginary, std(N_real)=delta_n; 
SNR=10*log10(Amp^2/(2*delta_n^2)); % note: it is only peak SNR dB 

Dr=2; % the ceil scale of each parameter, note: the more broarder of range of the parameter ( x,vx,y,vy,r,d,b); the worse of the initial accuracy
Dd=0.5; % for jump and sensitive doppler, the Dd need bigger enough
Db=0.04; % a smaller Db and Dr , increase the accuracy of the trajectory 

r=280:Dr:330; % the range of each parameter
d=-0.5:Dd:2.5;
b=0.5:Db:2.5; 

b3=-2.5:Db:-0.5; 
b4=-2.5:Db:-0.5; 
b5=-2.5:Db:-0.5; 
b6=-2.5:Db:-0.5; 
b7=-2.5:Db:-0.5; 
b9=0.5:Db:2.5; 
Nr=length(r); % the number of each ceil bin ; in here: 51*7*26
Nd=length(d);
Nb=length(b);
QQ=1;

Lr=1/Dr*QQ;
Ld=1/Dd*QQ;
Lb=1/Db*QQ; 

h_rdb=zeros(Nr,Nd,Nb,steps);% signal contribution to the ceil 
z_rdb_I=zeros(Nr,Nd,Nb,steps); % real part of z (complex amplitude)
z_rdb_Q=zeros(Nr,Nd,Nb,steps); % imaginary part of z
z_rdb=zeros(Nr,Nd,Nb,steps); % amplitude of z 
fai=0; % unknown phase (0,2pi), for simplicity, let fai=0;


% [z_rdb1]=computRealTMul(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1,d,b,r,Dd,Lb,fai,TE);
%获取总的量测
[z_rdb11,h_rdb]=computRealTMulcai(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1,d,b,r,Dd,Lb,fai,TE);
%surf(z_rdb11)
%[z_rdb11]=computRealT(interval1,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z,d,b,r,Dd,Lb,fai,steps)
%目标二的单独量测
%[z_rdb12]=computRealTMul2(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1,d,b,r,Dd,Lb,fai,TE);

% ------------------------ SIR based TBD ----------------

N=5000; % paricle number
xmin=250;    xmax=340; % local state space 
ymin=-100;   ymax=-10;
vxmin=-2;   vxmax=2.5;
vymin=-2;   vymax=2;

xp1=zeros(4,steps,N); %粒子群1
zp1=zeros(3,steps,N); %目标一的三维量测矩阵
xp2=zeros(4,steps,N); %粒子群2
zp2=zeros(3,steps,N); %目标二的三维量测矩阵

h_rdb_p1=zeros(steps,N); % blurring prediction
z_rdb_p1=zeros(steps,N); % amplitude prediction

h_rdb_p2=zeros(steps,N);
z_rdb_p2=zeros(steps,N);


xe1=zeros(4,steps); % state estimate
axe1=zeros(4,steps);
xe2=zeros(4,steps); 
axe2=zeros(4,steps);

Pb1=zeros(1,steps); % existing probability
Pb2=zeros(1,steps);

E1=zeros(steps,N); % existing variable
E2=zeros(steps,N);
qq1=ones(steps,N); % weight 
qq2=ones(steps,N);

EE1=[0.9 0.1 % E11,E10,E01,E00
    0.1 0.9];
delta_v_p1=1.5; % process noise std, CV model, using the same F_cv and G % it is a very important parameter for high maneuvering target
for i=1:N % initialization
    xp1(1,1,i)=xmin+(xmax-xmin)*rand;
    xp1(2,1,i)=vxmin+(vxmax-vxmin)*rand;
    xp1(3,1,i)=ymin+(ymax-ymin)*rand;
    xp1(4,1,i)=vymin+(vymax-vymin)*rand;
    
    xp2(1,1,i)=xmin+(xmax-xmin)*rand;
    xp2(2,1,i)=vxmin+(vxmax-vxmin)*rand;
    xp2(3,1,i)=ymin+(ymax-ymin)*rand;
    xp2(4,1,i)=vymin+(vymax-vymin)*rand;
    
    % -------examine if it is located in the valid space of r,d,b ,if not, resample---------
    qq(1,i)=1/N;
    if rand<0.1
        E1(1,i)=1;
    else
        E1(1,i)=0;
    end
    xe1(:,1)=xe1(:,1)+E1(1,i)*xp1(:,1,i); % initial estimation
    axe1(:,1)=axe1(:,1)+E1(1,i)*xp1(:,1,i);
    Pb1(1,1)=Pb1(1,1)+E1(1,i); % initial probability
    
     if rand<0.1
        E2(1,i)=1;
    else
        E2(1,i)=0;
    end
    xe2(:,1)=xe2(:,1)+E2(1,i)*xp2(:,1,i); % initial estimation
    axe2(:,1)=axe2(:,1)+E2(1,i)*xp2(:,1,i);
    Pb2(1,1)=Pb2(1,1)+E2(1,i); % initial probability
end

if Pb1(1,1)~=0;
   xe1(:,1)=xe1(:,1)/Pb1(1,1);
end
Pb1(1,1)=Pb1(1,1)/N;
axp1=xp1;

if Pb2(1,1)~=0;
   xe2(:,1)=xe2(:,1)/Pb2(1,1);
end
Pb2(1,1)=Pb2(1,1)/N;
axp2=xp2;

% --------------------------------- test ----------------------------------
steps_temp1=steps;%测试开始，采样次数
% --------------------------------- test ----------------------------------


for k=1:steps_temp1-1
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
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        q_sum1=sum(qq1(k+1,:));    
    for i=1:N
        qq1(k+1,i)=qq1(k+1,i)/q_sum1;
    end
    xp_temp1=xp1;
    axp_temp1=axp1;
    diffr=axp1(:,k+1,:)-xp1(:,k+1,:);
    EE_temp1=E1(k+1,:);
    for i=1:N%目标一重采样
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
%   q(k+1,i)=1;
    end
    %xp为重采样完粒子群
    wwww=1;
     for i=1:N
        Pb1(1,k+1)=Pb1(1,k+1)+E1(k+1,i);
         axe1(:,k+1)=axe1(:,k+1)+E1(k+1,i)*axp1(:,k+1,i);
    end
    if Pb1(1,k+1)~=0;
       axe1(:,k+1)=axe1(:,k+1)/Pb1(1,k+1);
    end 
    Pb1(1,k+1)=Pb1(1,k+1)/N;    

 %目标二
  for i=1:N
% -------------------- update existing variable E(k+1,i)
        if E2(k,i)==0    
            if rand<0.1
                E2(k+1,i)=1;
            else
                E2(k+1,i)=0;
            end
        else
            if rand<0.1
                E2(k+1,i)=0;
            else
                E2(k+1,i)=1;
            end
        end
% ---------------------- update particle xp(:,k+1,i)
        if E2(k,i)==1&&E2(k+1,i)==1   % paricle maintain 
            axp2(:,k+1,i)=F_cv*axp2(:,k,i)+G*delta_v_p1*randn(2,1);
        elseif E2(k,i)==0&&E2(k+1,i)==1 % new particle, note: we can choose a better zone for new particles 
            axp2(1,k+1,i)=xmin+(xmax-xmin)*rand;
            axp2(2,k+1,i)=vxmin+(vxmax-vxmin)*rand;
            axp2(3,k+1,i)=ymin+(ymax-ymin)*rand;
            axp2(4,k+1,i)=vymin+(vymax-vymin)*rand;
        elseif E2(k+1,i)==0
            axp2(:,k+1,i)=0;
        end
 xp2(:,k+1,i)=axp2(:,k+1,i);
 
 xpt=[304,0.3,-45,0.2];
qqt(1,k+1)=computweight( xpt,LR1,ST(k+1),RT11(k+1),E2(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11(:,:,:,k+1),r,d,b,Nb,Ld);

 xpt1=[304,0.3,-55,0.2];
qqt(2,k+1)=computweight( xpt1,LR1,ST(k+1),RT11(k+1),E2(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11(:,:,:,k+1),r,d,b,Nb,Ld);

xpt2=[304,0.3,-80,0.2];
qqt(3,k+1)=computweight( xpt2,LR1,ST(k+1),RT11(k+1),E2(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11(:,:,:,k+1),r,d,b,Nb,Ld);

qq2(k+1,i)=computweight(xp2(:,k+1,i),LR1,ST(k+1),RT11(k+1),E2(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb11(:,:,:,k+1),r,d,b,Nb,Ld);
 %以下为粒子权值进一步处理过程，将与已探测的目标接近（距离小于10）的粒子权值记为1，即按噪声处理
            if (((axe1(1,k+1)-xp2(1,k+1,i))^2+(axe1(3,k+1)-xp2(3,k+1,i))^2)<=20^2)
                qq2(k+1,i)=1;
            end
  
     
  end 
   
  
   
  
  %==========================计算完权重====================================
q_sum2=sum(qq2(k+1,:));
       for i=1:N
        qq2(k+1,i)=qq2(k+1,i)/q_sum2;
    end     
    xp_temp2=xp2;
    axp_temp2=axp2;
    diffr2=axp2(:,k+1,:)-xp2(:,k+1,:);
    EE_temp2=E2(k+1,:);  
        for i=1:N
        uu=rand;
        qq_sum2=0;
        for j=1:N
            qq_sum2=qq_sum2+qq2(k+1,j);
            if qq_sum2>=uu
               xp2(:,k+1,i)=xp_temp2(:,k+1,j);
                axp2(:,k+1,i)=axp_temp2(:,k+1,j);
                E2(k+1,i)=EE_temp2(j);
                break;
            end
        end
%   q(k+1,i)=1;
    end 
  
    for i=1:N
        Pb2(1,k+1)=Pb2(1,k+1)+E2(k+1,i);
         axe2(:,k+1)=axe2(:,k+1)+E2(k+1,i)*axp2(:,k+1,i);
    end
    if Pb2(1,k+1)~=0;
%        xe(:,k+1)=xe(:,k+1)/Pb(1,k+1);
       axe2(:,k+1)=axe2(:,k+1)/Pb2(1,k+1);
    end 
    Pb2(1,k+1)=Pb2(1,k+1)/N;   
    
    
    
    
    
    
    
end
%------------------------------------------test_end---------------------------------------------------------


zm1=zeros(3,steps);
for k=2:steps_temp1
    zm1(1,k)=sqrt((xe1(1,k))^2+(xe1(3,k))^2);
    zm1(2,k)=(xe1(1,k)*xe1(2,k)+xe1(3,k)*xe1(4,k))/zm1(1,k);
    zm1(3,k)=atan(xe1(1,k)/xe1(3,k));
end
PT1=zeros(1,steps);
tt1=0;

for k=1:steps
    if Pb1(k)>0.6;
        PT1(k)=1;
        tt1=tt1+1;
        Fxe1(:,tt1)=axe1(:,k);
    end
   
    
end

APb1(round,:)=Pb1;
APb2(round,:)=Pb2;
lengthINter1=size(intervalN11);
for i=1:lengthINter1(2)
%     RM(i)=sqrt((xe(1,intervalN(i))-x(1,intervalN(i)))^2+(xe(3,intervalN(i))-x(3,intervalN(i)))^2);
    aRM1(round,i)=sqrt((axe1(1,intervalN11(i))-ax1(1,intervalN11(i)))^2+(axe1(3,intervalN11(i))-ax1(3,intervalN11(i)))^2);
end
%目标二
lengthINter2=size(intervalN12);
for i=1:lengthINter2(2)
%     RM(i)=sqrt((xe(1,intervalN(i))-x(1,intervalN(i)))^2+(xe(3,intervalN(i))-x(3,intervalN(i)))^2);
    aRM2(round,i)=sqrt((axe2(1,intervalN12(i))-ax2(1,intervalN12(i)))^2+(axe2(3,intervalN12(i))-ax2(3,intervalN12(i)))^2);
end







end
%%%%%%%%%%%目标一
for i=1:steps
    avePb1(i)=mean(APb1(:,i));%检测概率
end

for i=1:lengthINter1(2)
    aveRM1(i)=mean(aRM1(:,i));%距离
end
 %-------目标二 
for i=1:steps
    avePb2(i)=mean(APb2(:,i));%检测概率
end

for i=1:lengthINter2(2)
    aveRM2(i)=mean(aRM2(:,i));%距离
end

figure(1);
plot(1:steps_temp1,avePb1(1:steps_temp1),'ro-');
title('目标一检测概率')
figure(4);
plot(1:steps_temp1,avePb2(1:steps_temp1),'ko-');
title('目标二检测概率')
figure(2);
 plot(aveRM1,'k-');
 hold on
 plot(aveRM2,'b-');
 legend('目标一误差','目标二误差')
mean(aveRM1)
figure(3);
plot(axe1(1,intervalN11),axe1(3,intervalN11),'b*-',ax1(1,intervalN11),ax1(3,intervalN11),'ro-',ax2(1,intervalN12),ax2(3,intervalN12),'ko-');
hold on 
plot(axe2(1,intervalN12),axe2(3,intervalN12),'r*-',ax1(1,intervalN11),ax1(3,intervalN11),'ro-',ax2(1,intervalN12),ax2(3,intervalN12),'ko-');
title('设定权重目标跟踪效果图')
