function [cloriginal,clfittest,fittest,gfs]=SHPSOairfoil(genNo,p0,range,uinf,AOA,Npanel) 
% 参数调整
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
UB = max(p0-range,p0+range);
LB = min(p0-range,p0+range);
NP = 50;
initial_sample_size=50; 
mf=500;    % maximum number of exact evaluation -----> user defined: 1000(30D,50D,100D) / 8000(30D) /
Dim = 11;

% function [ gsamp1 ,time_cost] = RUN_SHPSO(runs, Dim, FUN, Xmin, Xmax, findex)
time_begin=tic;
% rand('seed',sum(100*clock)); 
warning('off');
addpath(genpath(pwd));

D=Dim;       % dimension -----> user defined: 30 / 50 / 100 /
% 
sn1=1;  
gfs=zeros(1,fix(mf/sn1));   %-------sampling point according to fitness evaluation for ploting the convergence curve
% 
pop_size=50;                % population size                
CE=zeros(mf,2);             % achive the exact fitness

   
runs=1;
for r=1:runs
    % Initialization procedure
    fitcount=0;
    sam=LB+(UB-LB).*lhsdesign(initial_sample_size,D);      % initial LHS samples
    fit=zeros(1,initial_sample_size);
    for i=1:initial_sample_size
        x = sam(i,:);
        fitness = FUN(x);
        fit(i)=fitness;  
        fitcount=fitcount+1;
        CE(fitcount,:)=[fitcount,fit(i)];
        if mod (fitcount,sn1)==0
            cs1=fitcount/sn1; gfs(1,cs1)=min(CE(1:fitcount,2));
        end
    end
    hisx=sam; hisf=fit;                                          % archive all the exact evaluated samples
%    
    % main loop
    [r,fitcount]         
    [~,sidx]=sort(fit);                                         % 对个体按适应值由小到大进行排序
    sam=sam(sidx,:);  fit=fit(sidx);                            % 排序后的样本
    lox=sam(1:pop_size,:);  lof=fit(1:pop_size);                % 选择前ss个最优样本
    % -------------- P S O search --------------------- 
    psam=lox; efit=lof;                                         % PSO初始种群
    ps=size(psam,1);                                            % 种群规模
    [fpso,hisx,hisf,fitcount,CE,gfs,pbest,pbestval,formula1,formula2]= SHPSO_clean_version_for_share(FUN,D,ps,LB,UB,psam,efit,hisx,hisf,fitcount,mf,CE,sn1,gfs); 
    %------------End P S O search-----------------------
    
    fittest = pbest(1,:);
    clfittest =-pbestval(1);
    if clfittest==cloriginal
        fittest=p0;
    end
%     %plotting the original airfoil vs. the evolved (optimized)
%     plotairfoil(fittest,'k')
%     axis('equal')
%     hold on
%     plotairfoil(p0,'r')
%     legend('Optimized','original')
%     xlabel('X/C')
%     ylabel('Y/C')
%     title('Airfoil shape')  
%     
%     
%     fprintf('Best fitness (PSO-final): %e\n',min(hisf));        % 注：这个是采样点，并非当前最优点
% % 
%     gsamp1(r,:)=gfs;
    
%     usage(r,:)=[formula1,formula2]; 

end    
% usage_ave=mean(usage); % for SHPSO
%%%%%%%%%%%%%%%%%%%%% Output options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
