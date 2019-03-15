function  particleFigure(k,axp,ax,N,q)


     
   if k>=15&&k<=17
         xp_1 = zeros(1,N);
         xp_3 = zeros(1,N);
        for i = 1:N
             xp_1(i) = axp(1,k,i);
             xp_3(i) = axp(3,k,i);
         end
         figure;
         set(gca,'FontSize',12);
         hold on
         plot(ax(1, k), ax(3, k), 'r.', 'markersize',30)   %系统状态位置
         hold on;
         hs = sort(q(k,:));%将权重按升序存储
         level1 = find(q(k,:)<=hs(200));%将权重属于第一等级的粒子的索引存储
         level2 = find(q(k,:)>hs(200)&q(k,:)<=hs(400));
         level3 = find(q(k,:)>hs(400)&q(k,:)<=hs(600));
         level4 = find(q(k,:)>hs(600)&q(k,:)<=hs(800));
         level5 = find(q(k,:)>hs(800)&q(k,:)<=hs(1000));
         plot(xp_1(level1),xp_3(level1), 'g.','markersize',10);   %各个粒子位置
         hold on;
         plot(xp_1(level2),xp_3(level2), 'b.','markersize',10);   %各个粒子位置
         hold on;
         plot(xp_1(level3),xp_3(level3), 'r.','markersize',10);   %各个粒子位置
         hold on;
         plot(xp_1(level4),xp_3(level4), 'y.','markersize',10);   %各个粒子位置
         hold on;
         plot(xp_1(level5),xp_3(level5), 'k.','markersize',10);   %各个粒子位置
%          axis([290,306,-95,-70]);
%          set(gca,'XTick',290:0.5:306);
%          set(gca,'YTick',-95:1:-70);
         grid on;
         legend('True State', 'Particles');
   end
