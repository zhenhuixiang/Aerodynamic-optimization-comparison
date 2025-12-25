function [cloriginal,clfittest,fittest]=SHDEairfoil(genNo,p0,range,uinf,AOA,Npanel) 
% 参数调整
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
UB = max(p0-range,p0+range);
LB = min(p0-range,p0+range);
NP = 30;
Dim = 11;
Max_NFEs = 300;


% function [ gsamp1 ,time_cost] = RUN_DSSE(FUN, runs, Dim, findex, LB, UB, opt_f, err) 
time_begin=tic; 
warning('off'); 
%% 1.Algorithm parameters setting
runs = 1;
NP=30;                                      % Population size   
CE=zeros(Max_NFEs,2);                       % Achive the exact fitness
sn1=1;                                      % Compression parameters of convergence process
gfs=zeros(1,fix(Max_NFEs/sn1));             % Sampling point according to fitness evaluation for ploting the convergence curve

initial_sample_size=50;   
gs = 50;

%% 2.Runs Algorithm
for r=1:runs
    NFEs = 0;
    show =1;
    hisx = [];
    hisf = [];
   
    % 1. Initialization procedure
    sam=LB+(UB-LB).*lhsdesign(initial_sample_size,Dim);         % Initial LHS samples
    for i = 1:initial_sample_size   
        fit(i)=FUN(sam(i,:));
        NFEs=NFEs+1;
        if NFEs <= Max_NFEs 
            CE(NFEs,:)=[NFEs,fit(i)];
            if mod (NFEs,sn1)==0
                cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
            end
        end     
    end
    hisx=[hisx; sam]; hisf=[hisf fit]';                        % archive all the fitness evaluated samples
    [lhf, id]=sort(hisf);  
    POP = hisx(id(1:NP),:); POP_fitness = lhf(id(1:NP));
    
    gen = 0;
    Strategy = 1;
    while NFEs < Max_NFEs   
        gen = gen + 1;
        if Strategy == 1
            % Sample training samples, R B F network search and save
            [ghf,id]=sort(hisf);               % 使用历史数据库中最优的若干个体构造全局代理模型
            ghf=ghf(1:gs);     ghx=hisx(id(1:gs),:);
            ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
            spr=max(max(ghxd))/(Dim*gs)^(1/Dim);
            net=newrbe(ghx',ghf',spr);        % newrbe
            modelFUN=@(x) sim(net,x');
            maxgen=50*Dim; minerror=1e-6;
            [NewSample,~] = DE(Dim,maxgen,modelFUN,minerror,ghx);    
            NewSampleFitness=FUN(NewSample); 
            NFEs = NFEs + 1;
            if NFEs <= Max_NFEs 
                CE(NFEs,:)=[NFEs,NewSampleFitness];
                if mod (NFEs,sn1)==0
                    cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
                end
            end
            hisx = [hisx; NewSample];
            hisf = [hisf; NewSampleFitness];

            % update model
            [ghf,id]=sort(hisf);               % 使用历史数据库中最优的若干个体构造全局代理模型
            ghf=ghf(1:gs);     ghx=hisx(id(1:gs),:);
            ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
            spr=max(max(ghxd))/(Dim*gs)^(1/Dim);
            net=newrbe(ghx',ghf',spr);        % newrbe
            modelFUN=@(x) sim(net,x');

            % select
            [~, pointIndex] = max(POP_fitness); 
            pointFitness = POP_fitness(pointIndex);
            if pointFitness > NewSampleFitness
                POP(pointIndex,:) = NewSample;
                POP_fitness(pointIndex) = NewSampleFitness;
            end
            fprintf(1,'Iteration: %d,  No.evaluation: %d,  NewSampleFitness: %e,  \n',gen,NFEs,NewSampleFitness);
        end
        
        if Strategy == 2
            % DE 
            [~,id]=sort(POP_fitness);  
            pbest = POP(id(1),:);
            pbestFitness = POP_fitness(id(1));
            P = POP;
            
            CR = 0.9;
            for i=1:NP  
                F1 = rand(1)*0.6+0.3;
                F2 = rand(1)*0.6+0.3;
               %% mutation
                k0=randi([1,NP]);
                while(k0==i)
                    k0=randi([1,NP]);   
                end
                P1=P(k0,:);
                k1=randi([1,NP]);
                while(k1==i||k1==k0)
                    k1=randi([1,NP]);
                end
                P2=P(k1,:);
                k2=randi([1,NP]);
                while(k2==i||k2==k1||k2==k0)
                    k2=randi([1,NP]);
                end
                P3=P(k2,:);
                V(i,:)=P1+F1.*(P2-P3)+F2.*(pbest-P1);   

                %% bound
                for j=1:Dim
                  if (V(i,j)>UB(1,j)||V(i,j)<LB(1,j))
                     V(i,j)=LB(1,j)+rand*(UB(1,j)-LB(1,j));         
                  end
                end

                %% crossover
                jrand=randi([1,Dim]); 
                for j=1:Dim
                    k3=rand;
                    if(k3<=CR||j==jrand)
                        U(i,j)=V(i,j);
                    else
                        U(i,j)=P(i,j);      % 子代种群
                    end
                end
            end
            UFitness = modelFUN(U)';
            ssk = 0;
            for i = 1:NP
%                 pointNum = floor(5*NFEs/Max_NFEs);
                pointNum = floor(16*NFEs/Max_NFEs);
                pointIndex = [unidrnd(NP,1,pointNum) i];
                Ui_Fitness = repmat(UFitness(i),pointNum+1,1);           
                Pointsfitness = POP_fitness(pointIndex);
                if Ui_Fitness < Pointsfitness
                    UFitness(i) = FUN(U(i,:));
                    NFEs = NFEs + 1;
                    if NFEs <= Max_NFEs 
                        CE(NFEs,:)=[NFEs,UFitness(i)];
                        if mod (NFEs,sn1)==0
                            cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
                        end
                    end
                    ssk = ssk + 1;
                    hisx = [hisx; U(i,:)];
                    hisf = [hisf; UFitness(i)];
                    if UFitness(i) < POP_fitness(i)
%                         [~, pointIndex] = max(POP_fitness); 
                        POP(i,:) = U(i,:);
                        POP_fitness(i) = UFitness(i);
                    end
                end
            end
            if ssk == 0
                [~, pointIndex] = min(UFitness); 
                POP(pointIndex,:) = U(pointIndex,:);
                POP_fitness(pointIndex) = UFitness(pointIndex);
                ssk = ssk + 1;
            end
            fprintf(1,'Iteration: %d,  No.evaluation: %d,  Best: %e,  No.prescreen data: %d\n',gen,NFEs,-bestFitness,ssk);
        end
        Strategy = 2;
        if mod(gen,2) == 1
            Strategy = 1;
        end
        [~,id]=sort(hisf);               % 使用历史数据库中最优的若干个体构造全局代理模型
        bestP = hisx(id(1),:);
        bestFitness = hisf(id(1),:);
        
    end
    %---------------- End search --------------------
    endNFEs = NFEs;
    fittest = bestP;
    clfittest =-bestFitness;
    if clfittest==cloriginal
        fittest=p0;
    end
    %plotting the original airfoil vs. the evolved (optimized)
    plotairfoil(fittest,'k')
    axis('equal')
    hold on
    plotairfoil(p0,'r')
    legend('Optimized','original')
    xlabel('X/C')
    ylabel('Y/C')
    title('Airfoil shape')
    fprintf('Best fitness: %e\n',min(hisf));    
    gsamp1(r,:)=gfs;                    % record number of runs and optimization sample every time   
end    

%% 3.output
%%%%%%%%%%%%%%%%%%%%%%% Output options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluation index
best_samp   = min(gsamp1(:,end));
worst_samp  = max(gsamp1(:,end));
mean_samp   = mean(gsamp1(:,end));
median_samp = median(gsamp1(:,end));
std_samp    = std(gsamp1(:,end));
out1        = [best_samp,worst_samp,mean_samp,median_samp,std_samp];

% Convergence process (log10 and log are different!)
gsamp1_ave  = mean(gsamp1,1);
gsamp1_log  = log10(gsamp1_ave);
gsamplog    = log10(gsamp1);  

% index Compressed convergence process
for j=1:Max_NFEs
    if mod(j,sn1)==0
        j1=j/sn1; gener_samp1(j1)=j;
    end
end
