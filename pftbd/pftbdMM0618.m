% Objective:
%       implement the basic pftbd, maneuvering target
%       the measurement data is amplitude of each pixel bin(range, doppler,bearing) using point spread function 
%       a more practical function,i.e., window DFT processing measurement is planed to implement later
%       target appears and disappears in the scenario
%       ---------
%
%       conclusion: for point spread function, 
%                   1) to deal with maneuvering target, the process noise std should big enough with single model;
%                   2) if a parameter will jump, such as doppler, its ceil
%                   should big in simulation 
%                   3) a better tracking performance in maintain period
%                   often means a worse performance in the initial period with the same parameters
%                   4) the bluring parameter may be chosen as the inverse
%                   of its corresponding ceil scale,
%                   5) in simulation, we can choose a  proper ceil scale ,
%                   in practice, we need to change the bluring 
clear all;
close all;
%   -------------simulate the trajectory of target----------------
steps=20; % measurement length
T=1;% scan period T=1s
MC=1;
yuzhi=0.6;
for round=1:MC
    round
x=zeros(4,steps); % true state variable x=[x,vx,y,vy];
z=zeros(3,steps); % measurement variable z=[range,doppler,bearing];
time1_begin=5; % first turn start
time1_end=20; % first turn end

LR1=[0 -100];
DRbegin11=0;%时间差
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

t_appear=[2 2];% target appearing time
t_disappear=[10 10]; % target disappearing time
%interval=t_appear:t_disappear;% exist time of target
Oinix1=[300,0.3,-80,0.2]; % initial state %一个坐标系下的目标初始位置
Oinix2=[300,0.3,-30,0.2]; % initial state %一个坐标系下的目标初始位置
%有问题
[RT11 ,intervalN11, lengthRR11, detT11 ,x11, ST]=Predata(steps,Oinix1,t_appear(1),t_disappear(1),T,DRbegin11,LR1,F_cv,G,delta_v);
[RT12 ,intervalN12, lengthRR12, detT11 ,x12, ST]=Predata(steps,Oinix2,t_appear(2),t_disappear(2),T,DRbegin12,LR1,F_cv,G,delta_v);



ax1=x11;
ax1(1,:)=x11(1,:)+LR1(1);
ax1(3,:)=x11(3,:)+LR1(2);

ax2=x12;
ax2(1,:)=x12(1,:)+LR1(1);
ax2(3,:)=x12(3,:)+LR1(2);


z11=z;
compuz11=computT(x11,intervalN11);
z11(:,1:intervalN11(lengthRR11(2)))=compuz11;

z12=z;
compuz12=computT(x12,intervalN12);
z12(:,1:intervalN12(lengthRR12(2)))=compuz12;

z1=[z11;z12];

TE=zeros(2,steps);
TE(1,intervalN11(1):intervalN11(lengthRR11(2)))=1;
TE(2,intervalN12(1):intervalN12(lengthRR12(2)))=1;





interval=intervalN11;

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
[z_rdb11]=computRealTMul(Nr,Nd,Nb,delta_n,Lr,Dr,Ld,Db,Amp,z1,d,b,r,Dd,Lb,fai,TE);




% -------------- simulate data end and test-----------------
% ---------test output range-doppler-azimuth-amplitude image ----------
%
% ------------------------ SIR based TBD ----------------

N=1000; % paricle number
xmin=250;    xmax=340; % local state space 
ymin=-100;    ymax=-10;
vxmin=-1;   vxmax=1;
vymin=-1;   vymax=1;

xp1=zeros(4,steps,N); % state particle
zp1=zeros(3,steps,N); % range,doppler,bearing particle
h_rdb_p1=zeros(steps,N); % blurring prediction
z_rdb_p1=zeros(steps,N); % amplitude prediction
xe1=zeros(4,steps); % state estimate
axe1=zeros(4,steps);

Pb1=zeros(1,steps); % existing probability
% E1=zeros(steps,N); % existing variable
% qq1=ones(steps,N); % weight 

EE1=[0.9 0.1 % E11,E10,E01,E00
    0.1 0.9];

delta_v_p1=1.5; % process noise std, CV model, using the same F_cv and G % it is a very important parameter for high maneuvering target
[xp1, E1,qq]=creat(xmin,xmax,vxmax,vxmin,vymin,vymax,ymin,ymax,N);%初始化种群
%之后处理的内容
for i=1:N % initializatio 
    xe1(:,1)=xe1(:,1)+E1(i)*xp1(:,i); % initial estimation
    axe1(:,1)=axe1(:,1)+E1(i)*xp1(:,i);
    Pb1(1,1)=Pb1(1,1)+E1(i); % initial probability
end


if Pb1(1,1)~=0;
   xe1(:,1)=xe1(:,1)/Pb1(1,1);
end
Pb1(1,1)=Pb1(1,1)/N;

%之后处理的内容

axp1=xp1;
% ------- test ----------------
steps_temp1=steps;%测试开始，采样次数
% ------- test ----------------

TarNum=0;
xp=[];
menxian=10;
for k=1:steps_temp1-1
ma=size(xp);
axe1=[];
Allxp=[];
AllE=[];

%[pb,newaxp1,fE1,axe]=PDFM(xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp1,E1,z_rdb11(:,:,:,k+1),LR1,ST(k+1),RT11(k+1),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
%以下为原有目标的跟踪维持程序
% newParNum=[];
% newxp1=[];
% newE=[];
% xe1=[];
% pb=[];
newTarNum=0;%本时刻目标数
if TarNum~=0
ParNum=N*ones(1,TarNum);
% if k==12
%     cc=1;
% end
[newTar,newParNum,newxp1,newE,xe1,Pb,newTarNum]=PDFM4(yuzhi,ParNum,TarNum,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp,E1,z_rdb11(:,:,:,k),LR1,ST(k),RT11(k),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
%[pb,newxp1,fE1,xe1]=PDFM2SingT(yuzhi,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp1,E1,z_rdb11(:,:,:,k+1),LR1,ST(k+1),RT11(k+1),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
% Pb1(1,k+1)=pb;  %k+1时刻的目标出现概率     
axe1=xe1;%k+1时刻的目标位置
Allxp=newxp1;%跟踪粒子群
AllE1=newE;%跟踪粒子群信息
% if pb<0.6
%     TarNum=0;
%     break;
% end
%newTarNum
end
%以上为原有目标跟踪程序


    %以下为新目标探测程序
    if ma(2)==0||TarNum~=0%短路或
        [Dxp, DE1,qq]=creat(xmin,xmax,vxmax,vxmin,vymin,vymax,ymin,ymax,N);
    else
        Dxp=xp;
        DE1=E1;
    end
%     if TarNum==0 %上一时刻没有发现目标
%         %使用上一周的粒子群
%         
%     else%上一时刻有目标
%         %新生粒子群
%         [xp, E1,qq]=creat(xmin,xmax,vxmax,vxmin,vymin,vymax,ymin,ymax,N);%
%     end 
    %输入的粒子群不同
    [Pb,Detectxp,DetectE,Detectxe,DetectTarNum]=PDFDetect3(Dxp,DE1,yuzhi,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,z_rdb11(:,:,:,k),LR1,ST(k),RT11(k),Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
    %以下为与已有目标判断是否为同一目标
    if newTarNum~=0%如果已有目标
        for p=1:DetectTarNum
            flag=0;%不是新目标
            for q=1:TarNum
                dis(p,q)=sqrt((Detectxe(1,p)-xe1(1,q))^2+(Detectxe(3,p)-xe1(3,q))^2);%计算新目标与既有目标之间的距离
                if dis(p,q)<menxian%不是新目标
                    flag=1;
                    break;
                end
            end
          
            if flag==0%是新目标
                %保存新目标到跟踪目标相关信息
                   axe1=[axe1 Detectxe(:,p)];%估计值的记录
                   AllE1=[AllE1; DetectE((p-1)*N+1:p*N)];%目标存在概率的记录
                   Allxp=[Allxp;Detectxp((p-1)*4+1:p*4,:)]; 
                   newTarNum=newTarNum+1;
                   %Pb的保存未编
            end
        end
    else %如果之前没有目标，把所有新探测出的目标信息都储存起来
       
            Allxp=Detectxp;
            AllE1=DetectE;
            axe1=Detectxe;
            newTarNum=DetectTarNum;
        
    end
    
    %以上为与已有目标判断是否为同一目标
    TarNum=newTarNum;%这段程序的处理需要再注意
    vvv=1;
   xp=[];
   E1=[];
   xp=Allxp;
   E1=AllE1;
   
   
   
   
%Pb1(1,k)=pb; 
    %     ParNum=N;
%     if TarNum~=0
%         Pb1(1,k+1)=pb;  %k+1时刻的目标出现概率     
%         axe1(:,k+1)=xe1;%k+1时刻的目标位置
%         xp1=newxp1;
%         E1=fE1;
%      end
%     %以上为新目标探测程序
% end
%新的点迹与航迹的关联
% if TarNum~=0 %需要做关联
%    for p=1:newTarNum
% %新的点迹与航迹的关联，最近邻法
%       for q=1:TarNum
%         TarDis(p,q)=sqrt((axe1(p,1)-lastxe1(q,1))^2+(axe1(p,3)-lastxe1(q,3))^2);
%       end
%    end
%    
%    for p=1:newTarNum
%        [mp(p,:) np(p,:)]=sort(TarDis(p,:));
%        
%        
%    end
% end
   
TarNum=newTarNum;
end




















% zm1=zeros(3,steps);
% for k=2:steps_temp1
%     zm1(1,k)=sqrt((xe1(1,k))^2+(xe1(3,k))^2);
%     zm1(2,k)=(xe1(1,k)*xe1(2,k)+xe1(3,k)*xe1(4,k))/zm1(1,k);
%     zm1(3,k)=atan(xe1(1,k)/xe1(3,k));
% end
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

lengthINter1=size(intervalN11);
for i=1:lengthINter1(2)
%     RM(i)=sqrt((xe(1,intervalN(i))-x(1,intervalN(i)))^2+(xe(3,intervalN(i))-x(3,intervalN(i)))^2);
    aRM1(round,i)=sqrt((axe1(1,intervalN11(i))-ax1(1,intervalN11(i)))^2+(axe1(3,intervalN11(i))-ax1(3,intervalN11(i)))^2);
end
end

for i=1:steps
    avePb1(i)=mean(APb1(:,i));
end

for i=1:lengthINter1(2)
    aveRM1(i)=mean(aRM1(:,i));
end
  

% figure;
% plot(interval,zm(1,interval),'ro-',interval,z(1,interval),'b*-');title('range');
% figure;
% plot(interval,zm(2,interval),'ro-',interval,z(2,interval),'b*-');title('doppler');
% figure;
% plot(interval,zm(3,interval),'ro-',interval,z(3,interval),'b*-');title('bearing');
figure;
plot(1:steps_temp1,avePb1(1:steps_temp1),'ro-');
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
figure;

    plot(aveRM1);
mean(aveRM1)
figure;
plot(axe1(1,intervalN11),axe1(3,intervalN11),'b*-',ax1(1,intervalN11),ax1(3,intervalN11),'ro-',ax2(1,intervalN12),ax2(3,intervalN12),'ko-');
% figure;
%  plot(Fxe(1,:),Fxe(3,:),'b*-',ax(1,intervalN1),ax(3,intervalN1),'ro-');


















