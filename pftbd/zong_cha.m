for k=1:steps
  pai=0;  
for i=1:Nr % i corresponding to r
        for j=1:Nd % j corresponding to d
            for m=1:Nb % m corresponding to b
             pai=pai+1;   
z_cha(:,pai)=z_rdb12(i,j,m,k);
z_sum(:,pai)=z_rdb11(i,j,m,k);
            end
        end
end
  mm=1:1:26*7*51;
  figure;
  plot(mm,z_cha(:,mm),'r-')
  hold on;
  plot(mm,z_sum(:,mm),'b--')
  grid on
end