function www=pardist(xp,N)%画出粒子分布
% clear all
www=[];

for i=1:N
%     x=xp(1,1,i);
%     y=xp(3,1,i);
%     l=[x;y];
    www=[www xp(:,1,i)];
end
figure()
plot(www(1,:),www(3,:),'*');

% www=pardist(xp1(:,k+1,:),N);
