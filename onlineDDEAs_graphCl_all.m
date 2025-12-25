
close all;  
clear all;
clc;  
DE=load('result/500_DEairfoil runs=30.mat'); 
SHPSO=load('result/500_SHPSOairfoil runs=30.mat');
ESAO=load('result/500_ESAOairfoil runs=30.mat');
DESO=load('result/500_DESOairfoil runs=30.mat');
TSDDEO=load('result/500_TSDDEOairfoil runs=30.mat');
ESA=load('result/500_ESAairfoil runs=30.mat');

p0 = [0.0216 0.3445 0.07912 -0.6448 0.17 -0.033797 0.6748 0 0 0 0];
Npanel=20000;
uinf=1;
AOA=5*pi/180;
genNo=10;

DE_p = DE.fittest
SHPSO_p = SHPSO.fittest
ESAO_p = ESAO.fittest
DESO_p = DESO.fittest
TSDDEO_p = TSDDEO.fittest
ESA_p = ESA.fittest

figure
% %plotting the original airfoil vs. the evolved (optimized)
% plotairfoil(DE_p,'k')
% hold on
% plotairfoil(SHPSO_p,'y')
% hold on
% plotairfoil(ESAO_p,'g')
% hold on
% plotairfoil(DESO_p,'r')
% hold on
% plotairfoil(TSDDEO_p,'b')
% hold on
% plotairfoil(ESA_p,'b.')
% hold on
% legend('DE','ESAO','SHPSO','DESO', 'TSDDEO', 'ESA')
plotairfoil(p0,'k')
hold on
plotairfoil(DE_p,'b')
hold on
plotairfoil(ESA_p,'r.')
hold on
legend('p0','DE', 'ESA')
xlabel('X/C')
ylabel('Y/C')
title('Airfoil shape')

figure
%plotting the original airfoil vs. the evolved (optimized)

hold on
plotairfoil(p0,'k')
plotairfoil(DE_p,'g.-')
plotairfoil(ESA_p,'r.-')
% axis('equal')
legend('Original','DE','ESA')
xlabel('X/C')
ylabel('Y/C')
title('Airfoil shape')

