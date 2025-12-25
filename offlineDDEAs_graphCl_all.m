
close all;  
clear all;
clc;  

DDEASE=load('result/offline/DDEASE.mat'); 
BDDEALDG=load('result/offline/BDDEALDG.mat'); 
MSDDEO=load('result/offline/MSDDEO.mat'); 
MSEA=load('result/offline/MSEA.mat'); 

p0 = [0.0216 0.3445 0.07912 -0.6448 0.17 -0.033797 0.6748 0 0 0 0];
Npanel=20000;
uinf=1;
AOA=5*pi/180;
genNo=10;

[~,i] = sort(DDEASE.clfittest);
sort_result = DDEASE.fittest(i,:);
DDEASE_p = sort_result(10,:);

[~,i] = sort(BDDEALDG.clfittest);
sort_result = BDDEALDG.fittest(i,:);
BDDEALDG_p = sort_result(15,:);

[~,i] = sort(MSDDEO.clfittest);
sort_result = MSDDEO.fittest(i,:);
MSDDEO_p = sort_result(15,:);

[~,i] = sort(MSEA.clfittest);
sort_result = MSEA.fittest(i,:);
MSEA_p = sort_result(15,:);

figure
%plotting the original airfoil vs. the evolved (optimized)

hold on
plotairfoil(p0,'k')
plotairfoil(DDEASE_p,'y-')
plotairfoil(BDDEALDG_p,'b-')
plotairfoil(MSDDEO_p,'g.-')
plotairfoil(MSEA_p,'r.-')
% axis('equal')
legend('Original','DDEASE','BDDEALDG','MSDDEO','MSEA')
xlabel('X/C')
ylabel('Y/C')
title('Airfoil shape')

