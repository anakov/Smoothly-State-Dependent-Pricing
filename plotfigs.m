% Plots figures of steady-state objects

set(0,'DefaultFigureWindowStyle','docked')  % docks all figures

figure
% midriplot
colormap([0.73 0.83 0.96])
step = (hbound-lbound)/(length(prob)-1);
bar(lbound:step:hbound,prob,1)
title('Distribution of non-zero price changes')
xlabel('Size of price changes')
ylabel('Density of price changes')
xlim([lbound hbound])

figure
colormap([0.73 0.83 0.96])
bar(hazd)
xlabel('Months elapsed since last price adjustment')
title('Hazard rate: h(k) = f(k)/s(k-1)')
xlim([0.5 horizon+0.5])
ylim([0 0.5])

figure
colormap([0.73 0.83 0.96])
bar(MeanAbsPchangeTime)
xlabel('Months elapsed since last price adjustment')
title('Mean absolute price change as a function of time since adjustment')
xlim([0.5 horizon+0.5])


figure
plot(sgrid,optflexprice-log(wbar)-markup,':g')  
hold on
plot(sgrid,pstar-log(wbar)-markup,'r')
title('Optimal price policy')
xlabel('Log (inverse) productivity')
ylabel('Log target relative price')
axis tight

figure
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,D/MedV)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
 else
   plot(Pgrid-log(wbar)-markup,D/MedV)
   xlabel('Log relative price')
end
title('Adjustment gain')
axis tight

figure
if nums>1
%   map=colormap;
%   colormap([1 1 1; map])
   mesh(sgrid,Pgrid-log(wbar)-markup,PdistSmooth)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,Pdist)
   xlabel('Log relative price')
end
title('Stationary density of firms')
axis tight
%zl=zlim;

figure
%map=colormap;
%colormap([1 1 1; map])
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,Pdisteroded)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,Pdisteroded)
   xlabel('Log relative price')
end
title('Density of firms after mc shock and inflation')
axis tight
zl=zlim;
%zlim(zl)

figure
%map=colormap;
%colormap([1 1 1; map])
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,(1-lambda).*Pdisteroded)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,(1-lambda).*Pdisteroded)
   xlabel('Log relative price')
end
title('Density of non-adjusting firms')
axis tight
zlim(zl)

figure
%map=colormap;
%colormap([1 1 1; map])
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,lambda.*Pdisteroded)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
else
   plot(Pgrid-log(wbar)-markup,lambda.*Pdisteroded)
   xlabel('Log relative price')
end
title('Density of adjusting firms')
axis tight

figure
plot(Dgrid/MedV*100,Lvalues)
title('Lambda as a function of D')
xlabel('Loss from inaction (in % of firm''s median value)')
ylabel('Probability of adjustment')
ylim([0 1])

figure
plot(Dgrid2/MedV*100,Lvalues2)
title('Lambda as a function of L (relevant range)')
xlabel('Loss from inaction (in % of firm''s median value)')
ylabel('Probability of adjustment')
hold on
axis tight
xlim([0 1])


figure
bar(Dgrid2(1:end-1)/MedV*100,diff(Lvalues2))
title('Lambda pdf')
xlabel('Loss from inaction (in % of firm''s median value)')
ylabel('Probability of adjustment')
xlim([0 1])
   
if CalvoMenuMetric<1-10*eps^.5
   if nums>1
      figure
      mesh(sgrid,Pgrid-log(wbar)-markup,lambdaSmooth)
      xlabel('Log (inverse) productivity')
      ylabel('Log relative price')
   else
      plot(Pgrid-log(wbar)-markup,lambda)
      xlabel('Log relative price')
   end
   title('Probability of adjustment Lambda')
   axis tight
end

figure
bar(0:Dstep:1,density_D,1)
title('Realized losses and hazard function')
xlabel('Non-adjustment loss L as % of median V')
ylabel('Density of firms')
hold on
plot(Dgrid2/MedV*100,Lvalues2,'r')
xlim([-Dstep/2 1+Dstep/2])
ylim([0 0.5])

figure
bar(0:Dstep:1,density_Derod,1)
title('Potential losses and hazard function')
xlabel('Non-adjustment loss L as % of median V')
ylabel('Density of firms')
hold on
plot(Dgrid2/MedV*100,Lvalues2,'r')
% plot(0:Dstep:1,Lvalues3.*density_Derod,'r')
xlim([-Dstep/2 1+Dstep/2])
ylim([0 1])
 
figure
if nums>1
   mesh(sgrid,Pgrid-log(wbar)-markup,V)
   xlabel('Log (inverse) productivity')
   ylabel('Log relative price')
 else
   plot(Pgrid-log(wbar)-markup,V)
   xlabel('Log relative price')
end
title('Value function')
axis tight

% figure
% if nums>1
%    mesh(sgrid,Pgrid,PAYOFFMAT)
%    xlabel('Log real marginal cost')
%    ylabel('Log relative price')
%  else
%    plot(Pgrid,PAYOFFMAT)
%    xlabel('Log relative price')
% end
% title('Profit function')
% axis tight

% if CalvoMenuMetric<1-2*eps^.5
% 
% % figure
% % for cols=1:length(pstar)
% %     plot(100*(Pgrid-pstar(cols)),lambda(:,cols))
% %     hold on
% %     xlabel('Percent distance from optimal price')
% %     ylabel('Probability of adjustment')
% %     title('Probability of price increase (left) or decrease (right)')
% % end
% % hold off
% 
% figure
% %for cols=1:length(pstar)
%     plot(100*(log(PMAT')-pstar'*ones(1,length(Pgrid))),lambda')
% %    hold on
% %    pause(0.2)
%     xlabel('Percent distance from optimal price')
%     ylabel('Probability of adjustment')
%     title('Probability of price increase (left) or decrease (right)')
% %end
% %hold off
%     
% end

% figure
% bar(0:npstep:1,density_lambda,1)
% title('Density of adjustment probabilities')
% xlabel('Adjustment probability lambda')
% ylabel('Density of firms (after shock)')
% xlim([-npstep/2 1+npstep/2])
