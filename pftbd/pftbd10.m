 %8个雷达量测，目标出现在5-20；目标初始位置300 -80 速度vx：0.5

clear all;
close all;
%   -------------simulate the trajectory of target----------------
steps=30; % measurement length
T=1;% scan period T=1s
MC=1;
for round=1:MC
    tic;%%%新增-hua%%%
x=zeros(4,steps); % true state variable x=[x,vx,y,vy];4*30
z=zeros(3,steps); % measurement variable z=[range,doppler,bearing];3*30
time1_begin=5; % first turn start
time1_end=20; % first turn end

LR1=[0 -100];%雷达坐标
DRbegin1=0.2;
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

t_appear=5;% target appearing time
t_disappear=20; % target disappearing time
%interval=t_appear:t_disappear;% exist time of target
Oinix=[300,0.5,-80,0]; % initial state %同一坐标系下的目标初始位置

[RT1 ,intervalN1, lengthRR1, detT1 ,x1, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin1,LR1,F_cv,G,delta_v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x1是相对雷达校正（以雷达为原点）后的轨迹，且向后推了0.2s的状态
%真实轨迹且向后推了0.2s的状态
ax=x1;
ax(1,:)=x1(1,:)+LR1(1);
ax(3,:)=x1(3,:)+LR1(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1=z;
compuz1=computT(x1,intervalN1);%雷达一自己为坐标系原点获取三维量测compuz1(:,6:20)，注意输入是x1不是ax，x1是起始时刻是6.2的目标轨迹
z1(:,1:intervalN1(lengthRR1(2)))=compuz1;%z1(:,1:20)=compuz1(1:20)


[RT2 ,intervalN2, lengthRR2, detT2 ,x2, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin2,LR2,F_cv,G,delta_v);
z2=z;
compuz2=computT(x2,intervalN2);
z2(:,1:intervalN2(lengthRR2(2)))=compuz2;


[RT3 ,intervalN3, lengthRR3, detT3 ,x3, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin3,LR3,F_cv,G,delta_v);
z3=z;
compuz3=computT(x3,intervalN3);
z3(:,1:intervalN3(lengthRR3(2)))=compuz3;

[RT4 ,intervalN4, lengthRR4, detT4 ,x4, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin4,LR4,F_cv,G,delta_v);
z4=z;
compuz4=computT(x4,intervalN4);
z4(:,1:intervalN4(lengthRR4(2)))=compuz4;

[RT5 ,intervalN5, lengthRR5, detT5 ,x5, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin5,LR5,F_cv,G,delta_v);
z5=z;
compuz5=computT(x5,intervalN5);
z5(:,1:intervalN5(lengthRR5(2)))=compuz5;

[RT6 ,intervalN6, lengthRR6, detT6 ,x6, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin6,LR6,F_cv,G,delta_v);
z6=z;
compuz6=computT(x6,intervalN6);
z6(:,1:intervalN6(lengthRR6(2)))=compuz6;

[RT7 ,intervalN7, lengthRR7, detT7 ,x7, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin7,LR7,F_cv,G,delta_v);
z7=z;
compuz7=computT(x7,intervalN7);
z7(:,1:intervalN7(lengthRR7(2)))=compuz7;

[RT9 ,intervalN9, lengthRR9, detT9 ,x9, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin9,LR9,F_cv,G,delta_v);
z9=z;
compuz9=computT(x9,intervalN9);
z9(:,1:intervalN9(lengthRR9(2)))=compuz9;


interval=intervalN1;%[6 7 ... 20]
% ----------------trajectory end and test --------------------

%计算真的量测

Amp=sqrt(16); % signal amplitude of target; 4(3DB);8(6DB);16(9DB);32(12DB),assuming the amplitude is a constant
delta_n=1; % noise standard deviation N=N_real+j*N_imaginary, std(N_real)=delta_n; 
SNR=10*log10(Amp^2/(2*delta_n^2)); % note: it is only peak SNR dB 

Dr=2; % the ceil scale of each parameter, note: the more broarder of range of the parameter ( x,vx,y,vy,r,d,b); the worse of the initial accuracy
Dd=0.5; % for jump and sensitive doppler, the Dd need bigger enough
Db=0.04; % a smaller Db and Dr , increase the accuracy of the trajectory 

r=280:Dr:350; % 280：2：350the range of each parameter
d=-0.5:Dd:2.5;%-0.5：0.5：2.5
b=0.5:Db:2.5; %0.5：0.04：2.5

b3=-2.5:Db:-0.5;
b4=-2.5:Db:-0.5; 

b5=-2.5:Db:-0.5; 
b6=-2.5:Db:-0.5; 

b7=-2.5:Db:-0.5; 
b9=0.5:Db:2.5; 
Nr=length(r);% 36   the number of each ceil bin ; in here: 51*41*21
Nd=length(d);%7
Nb=length(b);%51
QQ=1;

Lr=1/Dr*QQ;
Ld=1/Dd*QQ;
Lb=1/Db*QQ; 

h_rdb=zeros(Nr,Nd,Nb,steps);% signal contribution to the ceil 
z_rdb_I=zeros(Nr,Nd,Nb,steps); % real part of z (complex amplitude)
z_rdb_Q=zeros(Nr,Nd,Nb,steps); % imaginary part of z
z_rdb=zeros(Nr,Nd,Nb,steps); % amplitude of z 
fai=0; % unknown phase (0,2pi), for simplicity, let fai=0;
%每帧量测1:5噪声，6：20 目标存在，21：30 噪声 
[h_rdb1, z_rdb_I1, z_rdb_Q1, z_rdb1]=computRealT(intervalN1,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1,d,b,r,Dd,Lb,fai,steps);
[h_rdb2 ,z_rdb_I2,z_rdb_Q2, z_rdb2]=computRealT(intervalN2,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z2,d,b,r,Dd,Lb,fai,steps);
[h_rdb3, z_rdb_I3, z_rdb_Q3, z_rdb3]=computRealT(intervalN3,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z3,d,b3,r,Dd,Lb,fai,steps);
[h_rdb4, z_rdb_I4 ,z_rdb_Q4 ,z_rdb4]=computRealT(intervalN4,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z4,d,b4,r,Dd,Lb,fai,steps);
[h_rdb5, z_rdb_I5 ,z_rdb_Q5 ,z_rdb5]=computRealT(intervalN5,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z5,d,b5,r,Dd,Lb,fai,steps);
[h_rdb6, z_rdb_I6 ,z_rdb_Q6 ,z_rdb6]=computRealT(intervalN6,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z6,d,b6,r,Dd,Lb,fai,steps);
[h_rdb7, z_rdb_I7 ,z_rdb_Q7 ,z_rdb7]=computRealT(intervalN7,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z7,d,b7,r,Dd,Lb,fai,steps);
[h_rdb9, z_rdb_I9 ,z_rdb_Q9 ,z_rdb9]=computRealT(intervalN9,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z9,d,b9,r,Dd,Lb,fai,steps);




% -------------- simulate data end and test-----------------
% ---------test output range-doppler-azimuth-amplitude image ----------
%
% ------------------------ SIR based TBD ----------------

N=1000; % paricle number
% xmin=250;    xmax=340; % local state space 
% ymin=-130;    ymax=-70;
% vxmin=-0.2;   vxmax=0.2;
% vymin=-0.2;   vymax=0.2;
xmin=250;    xmax=340; % local state space 
ymin=-120;    ymax=-60;
vxmin=-2;   vxmax=2.5;
vymin=-2;   vymax=2;

xp=zeros(4,steps,N); % state particle
zp=zeros(3,steps,N); % range,doppler,bearing particle
h_rdb_p=zeros(steps,N); % blurring prediction
z_rdb_p=zeros(steps,N); % amplitude prediction
xe=zeros(4,steps); % state estimate
axe=zeros(4,steps);

Pb=zeros(1,steps); % existing probability
E=zeros(steps,N); % existing variable
q=ones(steps,N); % weight 

EE=[0.9 0.1 % E11,E10,E01,E00
    0.1 0.9];
delta_v_p=1.5; % process noise std, CV model, using the same F_cv and G % it is a very important parameter for high maneuvering target
for i=1:N % initialization
    xp(1,1,i)=xmin+(xmax-xmin)*rand;
    xp(2,1,i)=vxmin+(vxmax-vxmin)*rand;
    xp(3,1,i)=ymin+(ymax-ymin)*rand;
    xp(4,1,i)=vymin+(vymax-vymin)*rand;
    % -------examine if it is located in the valid space of r,d,b ,if not, resample---------
    q(1,i)=1/N;
    if rand<0.1
        E(1,i)=1;
    else
        E(1,i)=0;
    end
    xe(:,1)=xe(:,1)+E(1,i)*xp(:,1,i); % initial estimation
    axe(:,1)=axe(:,1)+E(1,i)*xp(:,1,i);
    Pb(1,1)=Pb(1,1)+E(1,i); % initial probability
end
if Pb(1,1)~=0;
   xe(:,1)=xe(:,1)/Pb(1,1);
end
Pb(1,1)=Pb(1,1)/N;
axp=xp;
% ------- test ----------------
steps_temp=steps;%测试开始，采样次数
% ------- test ----------------

%粒子状态的重新计算，坐标的加减及时间的推移




for k=1:steps_temp-1
%     for i=1:N
%  xp(1,k,i)=xp(1,k,i)-LR(1);
%  xp(3,k,i)=xp(3,k,i)-LR(2);
%  %以上为粒子相对雷达坐标的变换
%  %以下为计算粒子时间推移过程
% %  xp(1,k,i)=xp(1,k,i)+xp(2,k,i)*(ST(k)-RT(k));
% %  xp(3,k,i)=xp(3,k,i)+xp(4,k,i)*(ST(k)-RT(k));
%  
%     end

    
    for i=1:N
        % -------------------- update existing variable E(k+1,i)
        if E(k,i)==0    
            if rand<0.1
                E(k+1,i)=1;
            else
                E(k+1,i)=0;
            end
        else
            if rand<0.1
                E(k+1,i)=0;
            else
                E(k+1,i)=1;
            end
        end
        % ---------------------- update particle xp(:,k+1,i)
        if E(k,i)==1&&E(k+1,i)==1   % paricle maintain 
            axp(:,k+1,i)=F_cv*axp(:,k,i)+G*delta_v_p*randn(2,1);
        elseif E(k,i)==0&&E(k+1,i)==1 % new particle, note: we can choose a better zone for new particles 
            axp(1,k+1,i)=xmin+(xmax-xmin)*rand;
            axp(2,k+1,i)=vxmin+(vxmax-vxmin)*rand;
            axp(3,k+1,i)=ymin+(ymax-ymin)*rand;
            axp(4,k+1,i)=vymin+(vymax-vymin)*rand;
        elseif E(k+1,i)==0
            axp(:,k+1,i)=0;
        end
 %每一个粒子状态进行坐标转换 
%  if k==6&&E(k+1,i)==1
%      ch=1;
%  end
%对每个粒子针对每个雷达求权值
 %粒子位置推移
 
 xp(:,k+1,i)=axp(:,k+1,i);
 q1(k+1,i)=computweight(xp(:,k+1,i),LR1,ST(k+1),RT1(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb1(:,:,:,k+1),r,d,b,Nb,Ld);
 q2(k+1,i)=computweight(xp(:,k+1,i),LR2,ST(k+1),RT2(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb2(:,:,:,k+1),r,d,b,Nb,Ld);
 q3(k+1,i)=computweight(xp(:,k+1,i),LR3,ST(k+1),RT3(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb3(:,:,:,k+1),r,d,b3,Nb,Ld);
 q4(k+1,i)=computweight(xp(:,k+1,i),LR4,ST(k+1),RT4(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb4(:,:,:,k+1),r,d,b4,Nb,Ld);
 q5(k+1,i)=computweight(xp(:,k+1,i),LR5,ST(k+1),RT5(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb5(:,:,:,k+1),r,d,b5,Nb,Ld);
 q6(k+1,i)=computweight(xp(:,k+1,i),LR6,ST(k+1),RT6(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb6(:,:,:,k+1),r,d,b6,Nb,Ld);
 q7(k+1,i)=computweight(xp(:,k+1,i),LR7,ST(k+1),RT7(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb7(:,:,:,k+1),r,d,b7,Nb,Ld);
 q9(k+1,i)=computweight(xp(:,k+1,i),LR9,ST(k+1),RT9(k+1),E(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb9(:,:,:,k+1),r,d,b9,Nb,Ld);
%  xp(1,k+1,i)=xp(1,k+1,i)-LR1(1);
%  xp(3,k+1,i)=xp(3,k+1,i)-LR1(2);
%  %时间推移
%  xp(1,k+1,i)=xp(1,k+1,i)+xp(2,k+1,i)*(ST(k)-RT2(k));
% xp(3,k+1,i)=xp(3,k+1,i)+xp(4,k+1,i)*(ST(k)-RT2(k));
% %     
%         
%         
%         
%         % ------------------------ update the weight q(k+1,i)
%         if E(k+1,i)==0
%             q(k+1,i)=1;
%         else
%             zp(1,k+1,i)=sqrt((xp(1,k+1,i))^2+(xp(3,k+1,i))^2);
%             zp(2,k+1,i)=(xp(1,k+1,i)*xp(2,k+1,i)+xp(3,k+1,i)*xp(4,k+1,i))/zp(1,k+1,i);
%             zp(3,k+1,i)=atan(xp(1,k+1,i)/xp(3,k+1,i));
%             if zp(1,k+1,i)<r(1)||zp(1,k+1,i)>r(end)||zp(2,k+1,i)<d(1)||zp(2,k+1,i)>d(end)||zp(3,k+1,i)<b(1)||zp(3,k+1,i)>b(end) % ensure the particle is valid
%                 q(k+1,i)=1; 
%             else
%                 p=3; % local likelihood parameter, very important,
%                 ik=ceil((zp(1,k+1,i)-r(1))/Dr); % the ceil bin of the target signal, represent Nr
%                 jk=ceil((zp(2,k+1,i)-d(1))/Dd);
%                 mk=ceil((zp(3,k+1,i)-b(1))/Db);
%                 
%                 for ii=max(1,ik-p):min(Nr,ik+p)
%                     for jj=max(1,jk-p):min(Nd,jk+p)
%                         for mm=max(1,mk-p):min(Nb,mk+p)
% %                             h_rdb_p(k+1,i)=Amp*exp(-0.5*Lr*(r(ii)-zp(1,k+1,i))^2/Dr-0.5*Ld*(d(jj)-zp(2,k+1,i))^2/Dd-0.5*Lb*(b(mm)-zp(3,k+1,i))^2/Db); % note! in here, fai=0, h_rdb_p(k+1,i) is abs
% %                             z_rdb_p(k+1,i)=h_rdb_p(k+1,i)*z_rdb(ii,jj,mm,k);% z_rdb_p(ii,jj,mm,k) take abs 
% %                             q(k+1,i)=q(k+1,i)*exp(-h_rdb_p(k+1,i)^2/2/delta_n^2)*besseli(0,z_rdb_p(k+1,i)/delta_n^2);% note maybe besellk is correct, I do not know
%                               h_rdb_p(k+1,i)=Amp^2*exp(-Lr*(r(ii)-zp(1,k+1,i))^2/Dr-Ld*(d(jj)-zp(2,k+1,i))^2/Dd-Lb*(b(mm)-zp(3,k+1,i))^2/Db)+2*delta_n^2; 
%                               q(k+1,i)=q(k+1,i)*(2*delta_n^2/h_rdb_p(k+1,i))*exp((0.5/delta_n^2-1/h_rdb_p(k+1,i))*z_rdb(ii,jj,mm,k+1));
%                         end
%                     end
%                 end
%             end
%         end
    end   %遍历粒子结束
    
    q_sum1=sum(q1(k+1,:));
   %R1权重归一化
    for i=1:N
        q1(k+1,i)=q1(k+1,i)/q_sum1;
    end
    
        q_sum2=sum(q2(k+1,:));
    for i=1:N
        q2(k+1,i)=q2(k+1,i)/q_sum2;
    end
     q_sum3=sum(q3(k+1,:));
    for i=1:N
        q3(k+1,i)=q3(k+1,i)/q_sum3;
    end
       q_sum4=sum(q4(k+1,:));
    for i=1:N
        q4(k+1,i)=q4(k+1,i)/q_sum4;
    end
    
           q_sum5=sum(q5(k+1,:));
    for i=1:N
        q5(k+1,i)=q5(k+1,i)/q_sum5;
    end
    
             q_sum6=sum(q6(k+1,:));
    for i=1:N
        q6(k+1,i)=q6(k+1,i)/q_sum6;
    end
    
    q_sum7=sum(q7(k+1,:));
    for i=1:N
        q7(k+1,i)=q7(k+1,i)/q_sum7;
    end
    
       
             q_sum9=sum(q9(k+1,:));
    for i=1:N
        q9(k+1,i)=q9(k+1,i)/q_sum9;
    end
    for i=1:N
        q(k+1,i)=q1(k+1,i)*q2(k+1,i)*q3(k+1,i)*q4(k+1,i)*q5(k+1,i)*q6(k+1,i)*q7(k+1,i)*q9(k+1,i);
    end
    
%   q(k+1,:)=q9(k+1,:);
    
        q_sum=sum(q(k+1,:));
    for i=1:N%总的归一化
        q(k+1,i)=q(k+1,i)/q_sum;
    end
        
    
%         if k==6
%         ch=1;
%     end
    xp_temp=xp;
    axp_temp=axp;
    diffr=axp(:,k+1,:)-xp(:,k+1,:);
    
    EE_temp=E(k+1,:);
    for i=1:N
        uu=rand;
        qq_sum=0;
        for j=1:N
            qq_sum=qq_sum+q(k+1,j);%%%%-hua
            if qq_sum>=uu
               xp(:,k+1,i)=xp_temp(:,k+1,j);
                axp(:,k+1,i)=axp_temp(:,k+1,j);
                E(k+1,i)=EE_temp(j);
                break;
            end
        end
%         q(k+1,i)=1;
    end

    
    for i=1:N
        Pb(1,k+1)=Pb(1,k+1)+E(k+1,i);
%          xe(:,k+1)=xe(:,k+1)+E(k+1,i)*xp(:,k+1,i);
         axe(:,k+1)=axe(:,k+1)+E(k+1,i)*axp(:,k+1,i);
    end
    if Pb(1,k+1)~=0;
%        xe(:,k+1)=xe(:,k+1)/Pb(1,k+1);
       axe(:,k+1)=axe(:,k+1)/Pb(1,k+1);%估计状态
    end
    
    Pb(1,k+1)=Pb(1,k+1)/N;
   
    
    %%%%%%%%%%%%用来画出粒子k时刻的分布以及真实状态的位置-hua%%%%%%%%%%%%%%%%%
   % particleFigure(k,axp,ax,N,q); 
    
end   %pf过程结束

zm=zeros(3,steps);
for k=2:steps_temp
    zm(1,k)=sqrt((xe(1,k))^2+(xe(3,k))^2);
    zm(2,k)=(xe(1,k)*xe(2,k)+xe(3,k)*xe(4,k))/zm(1,k);
    zm(3,k)=atan(xe(1,k)/xe(3,k));
end
PT=zeros(1,steps);
tt=0;

for k=1:steps
    if Pb(k)>0.6;
        PT(k)=1;
        tt=tt+1;
        Fxe(:,tt)=axe(:,k);
    end
   
    
end

APb(round,:)=Pb;

lengthINter=size(intervalN1);%【1，15】
for i=1:lengthINter(2)
%     RM(i)=sqrt((xe(1,intervalN(i))-x(1,intervalN(i)))^2+(xe(3,intervalN(i))-x(3,intervalN(i)))^2);
    aRM(round,i)=sqrt((axe(1,intervalN1(i))-ax(1,intervalN1(i)))^2+(axe(3,intervalN1(i))-ax(3,intervalN1(i)))^2);
end
toc;%%%新增-hua%%%
end

for i=1:steps
    avePb(i)=mean(APb(:,i));
end

for i=1:lengthINter(2)
    aveRM(i)=mean(aRM(:,i));
end
  

% figure;
% plot(interval,zm(1,interval),'ro-',interval,z(1,interval),'b*-');title('range');
% figure;
% plot(interval,zm(2,interval),'ro-',interval,z(2,interval),'b*-');title('doppler');
% figure;
% plot(interval,zm(3,interval),'ro-',interval,z(3,interval),'b*-');title('bearing');
figure(1);
plot(1:steps_temp,avePb(1:steps_temp),'ro-');
title('300,0.5,-80,0.5-r9');
% figure;
% % plot(x(1,interval),x(3,interval),'ro-');
% % hold on
% 
% plot(xe(1,intervalN),xe(3,intervalN),'b*-',x(1,intervalN),x(3,intervalN),'ro-');
% Oriinter=interval+1;
% plot(xe(1,intervalN),xe(3,intervalN),'b*-',x(1,Oriinter),x(3,Oriinter),'ro-');

% for i=1:lengthINter(2)
%     if PT(i)
% %     RM(i)=sqrt((xe(1,intervalN(i))-x(1,intervalN(i)))^2+(xe(3,intervalN(i))-x(3,intervalN(i)))^2);
%     FRM(i)=sqrt((axe(1,intervalN1(i))-ax(1,intervalN1(i)))^2+(axe(3,intervalN1(i))-ax(3,intervalN1(i)))^2);
% end
figure(2);

    plot(aveRM);
    title('300,0.5,-80,0.5-r9');

% figure;
% plot(axe(1,intervalN1),axe(3,intervalN1),'b*-',ax(1,intervalN1),ax(3,intervalN1),'ro-');
% figure;
%  plot(Fxe(1,:),Fxe(3,:),'b*-',ax(1,intervalN1),ax(3,intervalN1),'ro-');


















