function [cloriginal,clfittest,fittest,gfs]=DEairfoil(genNo,p0,range,uinf,AOA,Npanel) 
% function [recRMSE, bestP, bestFitness, endNFEs]=DEairfoil(Dim, NP, Max_NFEs, UpperBound, LowerBound, FUN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function optimizes an airfoil shape based on the coeffecient of lift
%value
%genNo      number of generations to mate
%p0         Original airfoil to oprimize
%range      Randomizer range to vary the PARSEC parameters
%uinf       flow free stream velocity
%AOA        airfoil angle of attack
%Npanel     number of panels for the fitness function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
UpperBound = max(p0-range,p0+range);
LowerBound = min(p0-range,p0+range);
NP = 50;
Dim = 11;
Max_NFEs = 500;
%%genetic parameters
    %% 初始化
    G=1;%设置迭代器（当前迭代代数）
    F=0.5; 
    CR=0.5;
    NFEs =0;
    CE = [];
    sn1=1;  
    CE=zeros(Max_NFEs,2);
    UB=repmat((UpperBound),NP,1);
    LB=repmat((LowerBound),NP,1);
    
    P=(UB-LB).*rand(NP,Dim)+LB;%随机产生初始种群个体
    for i=1:NP
        fitnessP(i)=FUN(P(i,:));%计算种群个体适应值
        NFEs=NFEs+1;
        if NFEs <= Max_NFEs 
            CE(NFEs,:)=[NFEs,fitnessP(i)];
            if mod (NFEs,sn1)==0
                cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
            end
        end 
    end
    
    %记录适应度函数调用次数
    [fitnessBestP,indexBestP]=min(fitnessP);
    bestP=P(indexBestP,:);
    recRMSE(1:NP)=fitnessP;
    
    %% 总体大循环
    while NFEs<Max_NFEs
        %%变异操作+交叉操作
        %constraining the coeffecient of lift
        
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
            V(i,:) = max(V(i,:), LowerBound);
            V(i,:) = min(V(i,:), UpperBound);
            
%             for j=1:Dim
%               if (V(i,j)>UB(i,j)||V(i,j)<LB(i,j))
%                  V(i,j)=LB(i,j)+rand*(UB(i,j)-LB(i,j));         
%               end
%             end
            
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
            fitnessU(i)=FUN(U(i,:));%计算种群个体适应值
            NFEs=NFEs+1;
            
            %%选择操作
            if(fitnessU(i)<fitnessP(i))
                P(i,:)=U(i,:);
                fitnessP(i)=fitnessU(i);
 
                if(fitnessU(i)<fitnessBestP)
                   fitnessBestP=fitnessU(i);
                   bestP=U(i,:);
                end
            end
            
            % update CE for plotting
            CE(NFEs,:)=[NFEs,fitnessP(i)];
            if mod (NFEs,sn1)==0
                cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
            end
            
            recRMSE(NFEs)=fitnessP(i);
        end
        for i=1:length(fitnessU)
            if -fitnessU(i)<=cloriginal
                    fitnessU(i)=-cloriginal;
            end
        end
        disp(['Iteration ' num2str(G) '   Best Cl = ' num2str(-fitnessBestP) '   NFE=' num2str(NFEs)]);
        G=G+1;
        
    end
    gfs = gfs(1:Max_NFEs);
    bestFitness=fitnessBestP;
    endNFEs = NFEs;
    fittest = bestP;
    clfittest =-bestFitness;
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
end


    