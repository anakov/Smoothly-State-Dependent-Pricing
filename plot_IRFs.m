% Plots impulse-response functions

set(0,'DefaultFigureWindowStyle','docked')  % docks all figures

% alternating color for irf plots
linescolorGrid = {'r.-','gs-','bo-','k*'};       % linecol = 1, 2, 3 for 'r.-','gs-','bo-'
if adjtype==0 || ~exist('linecol','var') || linecol>3; 
    linecol=1;
elseif adjtype==2
    linecol=2;
elseif adjtype==3
    linecol=3;
else
    linecol=linecol+1;
end;

linescolor=linescolorGrid{linecol};
MarkerSize = 5;

figure(2)
subplot(4,3,1)
plot(scalefactor*(d_A_path+d_R_path),linescolor,'MarkerSize', MarkerSize)  % plots the non-zero impulse
title('Shock process')
hold on
%legend('Calvo','SDSP')

subplot(4,3,2)
plot(scalefactor*(PI_path-mu),linescolor,'MarkerSize',MarkerSize)
title('Inflation')
y = ylim;
hold on

subplot(4,3,3)
plot(scalefactor*(R_path-Rss),linescolor,'MarkerSize', MarkerSize) % this is the monthly interest rate
title('Nominal interest rate')
hold on

subplot(4,3,4)
plot(intensive_margin_path,linescolor,'MarkerSize', MarkerSize)
title('Intensive margin')
ylim(y)
hold on


subplot(4,3,5)
plot(extensive_margin_path,linescolor,'MarkerSize', MarkerSize)
title('Extensive margin')
hold on
ylim(y)
xlim([0 TT])

subplot(4,3,6)
plot(selection_effect_path,linescolor,'MarkerSize', MarkerSize)
title('Selection effect')
hold on
ylim(y)
xlim([0 TT])

subplot(4,3,7)
plot(scalefactor*(C_path/Cbar-1),linescolor,'MarkerSize', MarkerSize)
title('Consumption')
y = ylim;
hold on

subplot(4,3,8)
plot(scalefactor*(labor_path/Nbar-1),linescolor,'MarkerSize', MarkerSize)
title('Labor')
ylim(y)
hold on

subplot(4,3,9)
plot(scalefactor*(delta_wedge/WeightPriceDispers-1),linescolor,'MarkerSize', MarkerSize)
title('Price dispersion')
ylim(y)
hold on

subplot(4,3,10)
plot(scalefactor*ex_ante_real_interest_rate,linescolor,'MarkerSize', MarkerSize)
title('Real interest rate')
xlabel('Months')
hold on

subplot(4,3,11)
plot(scalefactor*(w_path/wbar-1),linescolor,'MarkerSize', MarkerSize)
title('Real wage')
xlabel('Months')
hold on

subplot(4,3,12)
plot(scalefactor*(m_path/mbar-1),linescolor,'MarkerSize', MarkerSize)
title('Real money holdings')
xlabel('Months')
hold on   
