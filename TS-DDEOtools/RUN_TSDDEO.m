%% This matlab code was written by Huixiang Zhen
%% Please refer with all questions, comments, bug reports, etc. to zhenhuixiang@cug.edu.cn

function [ gsamp1 ,time_cost] = RUN_TSDDEO(runs, D, FUN, LB, UB, fname)
time_begin=tic;
warning('off');
addpath(genpath(pwd));

for r=1:runs
    % main loop
    fprintf('\n');
    disp(['函数：', fname,' 运行次数', num2str(r)]);  
    fprintf('\n');
    [hisx,hisf,fitcount,mf,CE,gfs]= TSDDEO(FUN,D,LB,UB); 

    fprintf('Best fitness (PSO-final): %e\n',min(hisf));        % 注：这个是采样点，并非当前最优点
    gsamp1(r,:)=gfs(1:mf);
end    

%%%%%%%%%%%%%%%%%%%%% Output options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% usage_ave=mean(usage); 
best_samp   = min(gsamp1(:,end));
worst_samp  = max(gsamp1(:,end));
samp_mean   = mean(gsamp1(:,end));
samp_median = median(gsamp1(:,end));
std_samp    = std(gsamp1(:,end));
out1        = [best_samp,worst_samp,samp_mean,samp_median,std_samp];
gsamp1_ave  = mean(gsamp1,1);
gsamp1_log  = log10(gsamp1_ave);
gsamplog    = log10(gsamp1);

% Time Complexity
time_cost=toc(time_begin);
save(strcat('result/NFE',num2str(mf),'_',fname,' runs=',num2str(runs),' Dim=',num2str(D)));
end
