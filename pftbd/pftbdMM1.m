%两个目标
clear all;
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
% time1=time1_end-time1_begin;% time interval of the first turn
% w1=angle1*pi/180/(time1);% angle velocity of the first turn
% F_ct1=[1 sin(w1*T)/w1 0 -(1-cos(w1*T))/w1 % transition matrix of ct model
%     0 cos(w1*T) 0 -sin(w1*T)
%     0 (1-cos(w1*T))/w1 1 sin(w1*T)/w1
%     0 sin(w1*T) 0 cos(w1*T) ];
% time2=time2_end-time2_begin;% time interval of the second turn
% w2=angle2*pi/180/(time2);% angle velocity of the second turn
% F_ct2=[1 sin(w2*T)/w2 0 -(1-cos(w2*T))/w2 % transition matrix of ct model
%     0 cos(w2*T) 0 -sin(w2*T)
%     0 (1-cos(w2*T))/w2 1 sin(w2*T)/w2
%     0 sin(w2*T) 0 cos(w2*T) ];
G=[T^2/2 0; % process noise transition matrix, assuming the noise distribution is the same for x and y axis
    T 0
    0 T^2/2
    0 T];
delta_v=0.01; %process noise standard deviation here, equal to 0;
Qv=G*delta_v^2*G';% process noise covariance

t_appear=[2 5];% target appearing time
t_disappear=[20 20]; % target disappearing time
%interval=t_appear:t_disappear;% exist time of target
Oinix1=[300,0.3,-30,0.2]; % initial state %一个坐标系下的目标初始位置
Oinix2=[300,0.3,-80,0.2]; % initial state %一个坐标系下的目标初始位置
%有问题
[RT11 ,intervalN11, lengthRR11, detT11 ,x11, ST]=Predata(steps,Oinix1,t_appear(1),t_disappear(1),T,DRbegin11,LR1,F_cv,G,delta_v);
[RT12 ,intervalN12, lengthRR12, detT11 ,x12, ST]=Predata(steps,Oinix2,t_appear(2),t_disappear(2),T,DRbegin12,LR1,F_cv,G,delta_v);



ax1=x11;
ax1(1,:)=x11(1,:)+LR1(1);%恢复坐标
ax1(3,:)=x11(3,:)+LR1(2);

ax2=x12;
ax2(1,:)=x12(1,:)+LR1(1);
ax2(3,:)=x12(3,:)+LR1(2);


z11=z;
compuz11=computT(x11,intervalN11);%在雷达自己的坐标系下求三维的量测
z11(:,1:intervalN11(lengthRR11(2)))=compuz11;

z12=z;
compuz12=computT(x12,intervalN12);
z12(:,1:intervalN12(lengthRR12(2)))=compuz12;

z1=[z11;z12];

TE=zeros(2,steps);
TE(1,intervalN11(1):intervalN11(lengthRR11(2)))=1;
TE(2,intervalN12(1):intervalN12(lengthRR12(2)))=1;

% [RT2 ,intervalN2, lengthRR2, detT2 ,x2, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin2,LR2,F_cv,G,delta_v);
% z2=z;
% compuz2=computT(x2,intervalN2);
% z2(:,1:intervalN2(lengthRR2(2)))=compuz2;
% 
% 
% [RT3 ,intervalN3, lengthRR3, detT3 ,x3, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin3,LR3,F_cv,G,delta_v);
% z3=z;
% compuz3=computT(x3,intervalN3);
% z3(:,1:intervalN3(lengthRR3(2)))=compuz3;
% 
% [RT4 ,intervalN4, lengthRR4, detT4 ,x4, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin4,LR4,F_cv,G,delta_v);
% z4=z;
% compuz4=computT(x4,intervalN4);
% z4(:,1:intervalN4(lengthRR4(2)))=compuz4;
% 
% [RT5 ,intervalN5, lengthRR5, detT5 ,x5, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin5,LR5,F_cv,G,delta_v);
% z5=z;
% compuz5=computT(x5,intervalN5);
% z5(:,1:intervalN5(lengthRR5(2)))=compuz5;
% 
% [RT6 ,intervalN6, lengthRR6, detT6 ,x6, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin6,LR6,F_cv,G,delta_v);
% z6=z;
% compuz6=computT(x6,intervalN6);
% z6(:,1:intervalN6(lengthRR6(2)))=compuz6;
% 
% [RT7 ,intervalN7, lengthRR7, detT7 ,x7, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin7,LR7,F_cv,G,delta_v);
% z7=z;
% compuz7=computT(x7,intervalN7);
% z7(:,1:intervalN7(lengthRR7(2)))=compuz7;
% 
% [RT9 ,intervalN9, lengthRR9, detT9 ,x9, ST]=Predata(steps,Oinix,t_appear,t_disappear,T,DRbegin9,LR9,F_cv,G,delta_v);
% z9=z;
% compuz9=computT(x9,intervalN9);
% z9(:,1:intervalN9(lengthRR9(2)))=compuz9;

% z(1,intervalN1)=sqrt(x1(1,intervalN1).^2+x1(3,intervalN1).^2); % z1=range between target and sensor, sensor is located in the origin
% z(2,intervalN1)=(x1(1,intervalN1).*x1(2,intervalN1)+x1(3,intervalN1).*x1(4,intervalN1))./z(1,intervalN1); % z2=doppler,  note!!!: there is a jump in doppler measure
% z(3,intervalN1)=atan(x1(1,intervalN1)./x1(3,intervalN1)); % z3=bearing, from North (positive Y axis) to the trajectory 


% z(1,intervalN)=sqrt(x(1,intervalN).^2+x(3,intervalN).^2); % z1=range between target and sensor, sensor is located in the origin
% z(2,intervalN)=(x(1,intervalN).*x(2,intervalN)+x(3,intervalN).*x(4,intervalN))./z(1,intervalN); % z2=doppler,  note!!!: there is a jump in doppler measure
% z(3,intervalN)=atan(x(1,intervalN)./x(3,intervalN)); % z3=bearing, from North (positive Y axis) to the trajectory 



interval11=intervalN11;
interval12=intervalN12;
% ----------------trajectory end and test --------------------
% subplot(2,2,1);
% plot(x(1,interval),x(3,interval),'b*-');
% subplot(2,2,2);
% plot(interval,z(1,interval),'r*-');
% subplot(2,2,3);
% plot(interval,z(2,interval),'k*-');
% subplot(2,2,4);
% plot(interval,z(3,interval),'g*-');
% figure;
% subplot(2,1,1);
% plot(interval,x(2,interval),'r*-');
% subplot(2,1,2);
% plot(interval,x(4,interval),'b*-');
% -------------  simulate amplitude data ------------
%计算真的量测

Amp=sqrt(32); % signal amplitude of target; 4(3DB);8(6DB);16(9DB);32(12DB),assuming the amplitude is a constant
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
Nr=length(r); % the number of each ceil bin ; in here: 51*41*21
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
[z_rdb11]=computRealTMul(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1,d,b,r,Dd,Lb,fai,TE);
%[z_rdb11]=computRealT(interval1,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z,d,b,r,Dd,Lb,fai,steps)
%目标二的单独量测
[z_rdb12]=computRealT(interval12,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z12,d,b,r,Dd,Lb,fai,steps)
% [z_rdb11]=computRealTMul(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1(1:3,:),d,b,r,Dd,Lb,fai,TE(1,:));
% [z_rdb12]=computRealTMul(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1(4:6,:),d,b,r,Dd,Lb,fai,TE(2,:));

% [h_rdb2 ,z_rdb_I2,z_rdb_Q2, z_rdb2]=computRealT(intervalN2,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z2,d,b,r,Dd,Lb,fai,steps);
% [h_rdb3, z_rdb_I3, z_rdb_Q3, z_rdb3]=computRealT(intervalN3,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z3,d,b3,r,Dd,Lb,fai,steps);
% [h_rdb4, z_rdb_I4 ,z_rdb_Q4 ,z_rdb4]=computRealT(intervalN4,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z4,d,b4,r,Dd,Lb,fai,steps);
% [h_rdb5, z_rdb_I5 ,z_rdb_Q5 ,z_rdb5]=computRealT(intervalN5,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z5,d,b5,r,Dd,Lb,fai,steps);
% [h_rdb6, z_rdb_I6 ,z_rdb_Q6 ,z_rdb6]=computRealT(intervalN6,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z6,d,b6,r,Dd,Lb,fai,steps);
% [h_rdb7, z_rdb_I7 ,z_rdb_Q7 ,z_rdb7]=computRealT(intervalN7,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z7,d,b7,r,Dd,Lb,fai,steps);
% [h_rdb9, z_rdb_I9 ,z_rdb_Q9 ,z_rdb9]=computRealT(intervalN9,Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z9,d,b9,r,Dd,Lb,fai,steps);
% for k=1:interval(1)-1
%     for i=1:Nr % i corresponding to r
%         for j=1:Nd % j corresponding to d
%             for m=1:Nb % m corresponding to b
%                 n_I=delta_n*randn;
%                 n_Q=delta_n*randn;
%                 z_rdb_I(i,j,m,k)=n_I;
%                 z_rdb_Q(i,j,m,k)=n_Q;
% %                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);%amplitude
%                 z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;% power
%             end
%         end
%     end
% end
% for k=interval(1):interval(end)
%     for i=1:Nr
%         for j=1:Nd
%             for m=1:Nb
%                 n_I=delta_n*randn;
%                 n_Q=delta_n*randn;
%                 h_rdb(i,j,m,k)=exp(-0.5*Lr*(r(i)-z(1,k))^2/Dr-0.5*Ld*(d(j)-z(2,k))^2/Dd-0.5*Lb*(b(m)-z(3,k))^2/Db);
%                 z_rdb_I(i,j,m,k)=Amp*h_rdb(i,j,m,k)*cos(fai)+n_I;
%                 z_rdb_Q(i,j,m,k)=Amp*h_rdb(i,j,m,k)*sin(fai)+n_Q;
% %                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);%amplitude
%                 z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;
%             end
%         end
%     end
% end
% for k=interval(end)+1:steps
%     for i=1:Nr % i corresponding to r
%         for j=1:Nd % j corresponding to d
%             for m=1:Nb % m corresponding to b
%                 n_I=delta_n*randn;
%                 n_Q=delta_n*randn;
%                 z_rdb_I(i,j,m,k)=n_I;
%                 z_rdb_Q(i,j,m,k)=n_Q;
% %                 z_rdb(i,j,m,k)=sqrt(z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2);
%                 z_rdb(i,j,m,k)=z_rdb_I(i,j,m,k)^2+z_rdb_Q(i,j,m,k)^2;
%             end
%         end
%     end
% end



% -------------- simulate data end and test-----------------
% ---------test output range-doppler-azimuth-amplitude image ----------
%
% ------------------------ SIR based TBD ----------------

N=1000; % paricle number
xmin=250;    xmax=340; % local state space 
ymin=-100;   ymax=-10;
vxmin=-2;   vxmax=2.5;
vymin=-2;   vymax=2;

xp1=zeros(4,steps,N); % state particle
zp1=zeros(3,steps,N); % range,doppler,bearing particle
xp2=zeros(4,steps,N); % state particle
zp2=zeros(3,steps,N); % range,doppler,bearing particle

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
end
if Pb1(1,1)~=0;
   xe1(:,1)=xe1(:,1)/Pb1(1,1);
end
Pb1(1,1)=Pb1(1,1)/N;
axp1=xp1;
% ------- test ----------------
steps_temp1=steps;%测试开始，采样次数
% ------- test ----------------


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
   %%%%%%%%%%%%%%%%%%% 
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
%   q(k+1,i)=1;
    end
    %xp为重采样完粒子群
    wwww=1;
%     figure(1)
    
%     plot(xp1(1,k+1,:),xp1(3,k+1,:));
    
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
    
    
%%%%%%%%%%%%%%%%%%开始检测目标二
    if Pb1(1,k+1)>=0.8
      for i=1:N % initialization
    xp2(1,k,i)=xmin+(xmax-xmin)*rand;
    xp2(2,k,i)=vxmin+(vxmax-vxmin)*rand;
    xp2(3,k,i)=ymin+(ymax-ymin)*rand;
    xp2(4,k,i)=vymin+(vymax-vymin)*rand;
    % -------examine if it is located in the valid space of r,d,b ,if not, resample---------
    %qq(1,i)=1/N;
    if rand<0.1
        E2(k,i)=1;
    else
        E2(k,i)=0;
    end
    xe2(:,k+1)=xe2(:,1)+E2(k,i)*xp2(:,k,i); % initial estimation
    axe2(:,k+1)=axe2(:,1)+E2(k,i)*xp2(:,k,i);
    Pb2(1,k)=Pb2(1,k)+E2(k,i); % initial probability
    end 
if Pb2(1,k)~=0;
   xe2(:,k)=xe2(:,k)/Pb2(1,k);
end
Pb2(1,k)=Pb2(1,k)/N;
axp2=xp2;
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
 %重点修改
 qq2(k+1,i)=computweight(xp2(:,k+1,i),LR1,ST(k+1),RT11(k+1),E2(k+1,i),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,z_rdb12(:,:,:,k+1),r,d,b,Nb,Ld);
    end
        q_sum2=sum(qq2(k+1,:));
    for i=1:N
        qq2(k+1,i)=qq2(k+1,i)/q_sum2;
    end
        

    xp_temp2=xp2;
    axp_temp2=axp2;
    diffr=axp2(:,k+1,:)-xp2(:,k+1,:);
    EE_temp2=E2(k+1,:);
    for i=1:N
        uuu=rand;
        qq_sum2=0;
        for j=1:N
            qq_sum2=qq_sum2+qq2(k+1,j);
            if qq_sum2>=uuu
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
%          xe(:,k+1)=xe(:,k+1)+E(k+1,i)*xp(:,k+1,i);
         axe2(:,k+1)=axe2(:,k+1)+E2(k+1,i)*axp2(:,k+1,i);
    end

    if Pb2(1,k+1)~=0;
%        xe(:,k+1)=xe(:,k+1)/Pb(1,k+1);
       axe2(:,k+1)=axe2(:,k+1)/Pb2(1,k+1);
    end
    
    Pb2(1,k+1)=Pb2(1,k+1)/N;  
      
      
    end
    
    
    
    
    
    
    
    
    
    
    
end



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
    aRM2(round,i)=sqrt((axe2(1,intervalN12(i))-ax2(1,intervalN12(i)))^2+(axe1(3,intervalN12(i))-ax2(3,intervalN12(i)))^2);
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
% figure;
% plot(interval,zm(1,interval),'ro-',interval,z(1,interval),'b*-');title('range');
% figure;
% plot(interval,zm(2,interval),'ro-',interval,z(2,interval),'b*-');title('doppler');
% figure;
% plot(interval,zm(3,interval),'ro-',interval,z(3,interval),'b*-');title('bearing');
figure(1);
plot(1:steps_temp1,avePb1(1:steps_temp1),'ro-');
% hold on 
% plot(1:steps_temp2,avePb2(1:steps_temp2),'bo-');


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

 plot(aveRM1,'k-');
mean(aveRM1)
figure(3);
plot(axe1(1,intervalN11),axe1(3,intervalN11),'b*-',ax1(1,intervalN11),ax1(3,intervalN11),'ro-',ax2(1,intervalN12),ax2(3,intervalN12),'ko-');
hold on 
plot(axe2(1,intervalN11),axe2(3,intervalN11),'r*-',ax1(1,intervalN11),ax1(3,intervalN11),'ro-',ax2(1,intervalN12),ax2(3,intervalN12),'ko-');




% figure;
%  plot(Fxe(1,:),Fxe(3,:),'b*-',ax(1,intervalN1),ax(3,intervalN1),'ro-');


















