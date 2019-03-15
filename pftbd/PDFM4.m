function [newTar,newParNum,newxp,newE,axe1,Pb,NewTarNum]=PDFM4(yuzhi,ParNum,TarNum,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp,E,z_rdb11,LR1,ST,RT,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1)
newTar=ones(1,TarNum);
newTarNum=0;
%以下程序为多个跟踪粒子群跟踪多个目标
%TarNum为上周期跟踪目标数
%ParNUm为每个跟踪粒子群的粒子个数
Cnum=0;
E1=[];
newE=[];
Nmax=0;
newxp=[];
newParNum=[];
axe1=[];
NewTarNum=0;
flag=0;
for i=1:TarNum
    xp1=[];
    temp=[];
    xp1=xp((i-1)*4+1:(i-1)*4+4,1:ParNum(i));
    E1=E(Cnum+1:Cnum+ParNum(i));
    Cnum=Cnum+ParNum(i);
      [pb,newxp1,fE1,xe1]=PDFM2SingT(yuzhi,xmin,xmax,vxmin,vxmax,ymin,ymax,vymin,vymax,N,xp1,E1,z_rdb11,LR1,ST,RT,Dr,Db,Dd,Nr,Nd,Amp,Lr,Lb,delta_n,r,d,b,Nb,Ld,F_cv,G,delta_v_p1);
   
   if pb>yuzhi
       flag=1;
         Pb(i)=pb;
         NewTarNum=NewTarNum+1;
        if N>Nmax
            if Nmax==0
                newxp=[newxp;newxp1];
            else
                pp=size(newxp);
                temp=zeros(pp(1),N-Nmax);%
           newxp=[newxp temp];
           newxp=[newxp;newxp1];
            end
           Nmax=N;
        %保存本次跟踪粒子群大小
    else
        temp=zeros(4,Nmax-N);
        newxp1=[newxp1 temp];
        newxp=[newxp;newxp1];
    end
    newParNum=[newParNum;N];
     newE=[newE;fE1];
     axe1=[axe1 xe1];
     newTarNum=newTarNum+1;
     newTar(i)=1;

%         www=1;
   end
end
if flag==0
    Pb=0;
end
  