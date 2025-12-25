%% ================ DE ======================================%%%%
%%  
function [bestP,bestFitness] = DE(Dim, Max_NFEs,FUN,minerror,ghx)
    time_begin=tic;
%     disp('DE optimzer search');
    % e.g., d=20; maxfe=1000;
    % d: dimensionality
    % maxfe: maximal number of fitness evaluations
    n = Dim; NP = 50; flag_er=0;
    %parameter setting
    %parameter initiliaztion
    lu = [min(ghx); max(ghx)];
    LowerBound = lu(1, :);
    UpperBound = lu(2, :);
    
    %% 初始化
    G=1;%设置迭代器（当前迭代代数）
    F=0.5; 
    CR=0.5;
    UB=UpperBound;
    LB=LowerBound;
    
    UB=repmat((UB),NP,1);
    LB=repmat((LB),NP,1);
    
    P=(UB-LB).*rand(NP,Dim)+LB;%随机产生初始种群个体
    
    fitnessP=FUN(P);%计算种群个体适应值
        
    NFEs=NP;%记录适应度函数调用次数
    [fitnessBestP,indexBestP]=min(fitnessP);
    bestP=P(indexBestP,:);
    recRMSE(1:NP)=fitnessP;
    
    %% 总体大循环
    while NFEs<Max_NFEs
        %%变异操作+交叉操作
        fitnessBestP_old = fitnessBestP;
        for i=1:NP  
            
            %从当前种群P中随机选出P1，P2，P3
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
            
            V(i,:)=P1+F.*(P2-P3);   
            
            %%边界处理
%             V(i,:) = max(V(i,:), LowerBound);
%             V(i,:) = min(V(i,:), UpperBound);
            
            for j=1:Dim
              if (V(i,j)>UB(i,j)||V(i,j)<LB(i,j))
                 V(i,j)=LB(i,j)+rand*(UB(i,j)-LB(i,j));         
              end
            end
            
            %%交叉操作
            jrand=randi([1,Dim]); 
            for j=1:Dim
                k3=rand;
                if(k3<=CR||j==jrand)
                    U(i,j)=V(i,j);
                else
                    U(i,j)=P(i,j);
                end
            end
        end
        
        time_begin=tic;
        fitnessU=FUN(U);%计算种群个体适应值   向量计算比For循环快多了，调用一次FUN算一个个体的评价，和算种群所有的个体评价，速度是一样的
        time_cost=toc(time_begin);
        % disp(['local search time_cost=' num2str(time_cost)]);
        NFEs=NFEs+NP;
        
        for i=1:NP 
            %%选择操作
            if(fitnessU(i)<fitnessP(i))
                P(i,:)=U(i,:);
                fitnessP(i)=fitnessU(i);
                if(fitnessU(i)<fitnessBestP)
                   fitnessBestP=fitnessU(i);
                   bestP=U(i,:);
                end
            end
            recRMSE(NFEs)=fitnessP(i);
        end
        error=abs(fitnessBestP_old-fitnessBestP); 
        if error <= minerror   % minerror终止条件
            flag_er=flag_er+1;
        else
            flag_er=0;
        end
        if flag_er >=10
            break;
        end
        % disp(['Iteration ' num2str(G) '   Best Cost = ' num2str(fitnessBestP) '   NFE=' num2str(NFEs)]);
        G=G+1;
    end
    bestFitness=fitnessBestP;
    endNFEs = NFEs;
    
    time_cost=toc(time_begin);
    % disp(['local search time_cost=' num2str(time_cost)]);
end
    