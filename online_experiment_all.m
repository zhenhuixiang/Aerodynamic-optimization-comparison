clear;
clc;
%initial individual
% p0=[0.0155 0.296632 0.060015 -0.4515 0.296632 -0.06055 0.453 0 0.001260 0 7.36]; %%NACA 0012
p0=[0.0216 0.3445 0.07912 -0.6448 0.17 -0.033797 0.6748 0 0 0 0];  %%NACA 2411
range=[0.0015 0.025 0.015 -0.01 0.02 -0.015 0.075 0 0 -0.175 0.05];
%range=[0.02 0.023 0.32 0.37 0.077 0.08 -0.63 -0.65 0.15 0.19 -0.02 -0.05 0.6 0.75 0 0 0 0 -4.55 -4.9 15 15.1];
%Solver parameters
Npanel=200;
uinf=1;
AOA=5*pi/180;
genNo=10;       %number of generations
%Genetic solution
runs = 30;
funcArr={'DEairfoil'; 'SHPSOairfoil';  'ESAOairfoil'; 'DESOairfoil'; 'TSDDEOairfoil';   'ESAairfoil';  }; % functions

for j = 1: length(funcArr)
    fname = cell2mat(funcArr(j)); 
    FUN=@(genNo,p0,range,uinf,AOA,Npanel) feval(fname,genNo,p0,range,uinf,AOA,Npanel); 
    for i = 1:runs
        [cloriginal,clfittest,fittest,gfs]=FUN(genNo,p0,range,uinf,AOA,Npanel);
        gsamp1(i,:)=[-gfs(1:500)]; 
    end

    % Evaluation index
    worst_samp   = min(gsamp1(:,end));
    best_samp = max(gsamp1(:,end));
    mean_samp   = mean(gsamp1(:,end));
    median_samp = median(gsamp1(:,end));
    std_samp    = std(gsamp1(:,end));
    out1        = [best_samp,worst_samp,mean_samp,median_samp,std_samp];
    gsamp1_ave  = mean(gsamp1,1);
    gsamp1_log  = log10(gsamp1_ave);   
    gsamplog    = log10(gsamp1);  % for ackley function
    save(strcat('result/',num2str(500), '_', fname,' runs=', num2str(runs)));

%     %ploting and graphing
%     fprintf(' Original   Cl= %f \n Optimized  Cl= %f \n',cloriginal,clfittest)
%     figure
%     graphCl(fittest,Npanel,uinf,5,'k');
%     graphCl(p0,Npanel,uinf,5,'r');
%     legend('Optimized','original')
%     xlabel('AOA')
%     ylabel('Cl')
%     line([-100 100],[0,0],'color','k','LineWidth',1)
%     line([0,0],[-10 10],'color','k','LineWidth',1)
%     title('Coeffecient of lift vs. Angle of attack')
%     grid on
end