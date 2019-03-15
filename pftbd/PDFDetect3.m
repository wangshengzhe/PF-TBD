function [Pb,Detectxp,DetectE,Detectaxe,TarNum]=PDFDetect3(xp,E1,yuzhi,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,z_rdb11,LR1,ST,RT,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1)


%qqq=test(xp);
delL=2;
delV=0.1;
%以下程序为新目标探测程序
TarNum=0;

    Pb=[];%目标存在概率
    Detectxp=[];%跟踪粒子群
    DetectE=[];%粒子状态
    Detectaxe=[];%目标位置
[pb,newxp1,fE1,xe1]=PDFM2SingT(yuzhi,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp,E1,z_rdb11,LR1,ST,RT,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
% if pb>yuzhi%发现目标
%     TarNum=1;
% end
if pb>yuzhi%发现目标
   [clustercenter, clustern,newxp,TarParNum,newE,newI]=meanshiftPar(newxp1,fE1);%将重采样后的粒子群聚类区分目标
    newxp=newxp';
    TarNum=clustern;

    %编写各目标粒子群求目标位置
  
    for i=1:clustern
        %每个粒子群补全粒子个数，对目标状态进行估计
        sumNum=TarParNum(i);
        if i==1
            subxp=newxp(:,1:TarParNum(i));
            subE=newE(1:TarParNum(i));
        else
            subxp=newxp(:,sum(TarParNum(1,1:i-1))+1:sum(TarParNum(1,1:i)));
            subE=newE(sum(TarParNum(1,1:i-1))+1:sum(TarParNum(1,1:i)));
        end
        while sumNum<N%粒子个数少于规定大小
            %在原粒子群中随机复制粒子
           a=randperm(TarParNum(i));
           newPar=subxp(:,a(1));
           newParE=subE(a(1));
           newPar(1)=newPar(1)+delL*rand;
           newPar(2)=newPar(2)+delV*rand;
           newPar(3)=newPar(3)+delL*rand;
           newPar(4)=newPar(4)+delV*rand;
           subxp=[subxp newPar];
           subE=[subE newParE];
           sumNum=sumNum+1;
        end
         subpb=0; 
         subaxe1=zeros(4,1);
        for j=1:N
         subpb=subpb+subE(j);
         subaxe1=subaxe1+subE(j)*subxp(:,j);
        end
        if subpb~=0;
           subaxe1=subaxe1/subpb;
        end
        subpb=subpb/N;
        %每个目标跟踪完毕
        Pb=[Pb subpb];
        Detectxp=[Detectxp;subxp];
        DetectE=[DetectE subE];
        Detectaxe=[Detectaxe subaxe1];
    end
    
   %
else
    Detectxp=newxp1;
    DetectE=fE1;
    TarParNum=[];
    Pb=pb;
    
    
end

    xxxxx=1;
  