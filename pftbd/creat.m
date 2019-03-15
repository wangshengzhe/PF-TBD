function [xp1, E1,qq]=creat(xmin,xmax,vxmax,vxmin,vymin,vymax,ymin,ymax,N)
%
%   初始化种群
for i=1:N
    xp1(1,i)=xmin+(xmax-xmin)*rand;
    xp1(2,i)=vxmin+(vxmax-vxmin)*rand;
    xp1(3,i)=ymin+(ymax-ymin)*rand;
    xp1(4,i)=vymin+(vymax-vymin)*rand;
    qq(1,i)=1/N;
    if rand<0.5%0.1
        E1(i)=1;
    else
        E1(i)=0;
    end
end

end

