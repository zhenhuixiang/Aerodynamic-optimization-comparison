function [cloriginal,clfittest,fittest,gfs]=ESAOairfoil(genNo,p0,range,uinf,AOA,Npanel) 
% 参数调整
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
U_bound = max(p0-range,p0+range);
L_bound = min(p0-range,p0+range);
NP = 50;
initial_sample_size=50;   
Dim = 11;
Max_NFEs = 500;

G = 0;                                    %---- Number of Generations
CE=zeros(Max_NFEs,2);                       % Achive the exact fitness


gs = initial_sample_size;

% DE parameter
F=0.8;
CR=0.8;

% number of local search samples
ls = 2*Dim;

% number of calling Strategys 
Strategy1 = 0;
Strategy2 = 0;
runs = 1;
CE=zeros(Max_NFEs,2);                       % Achive the exact fitness
sn1=1;                                      % Compression parameters of convergence process
gfs=zeros(1,fix(Max_NFEs/sn1));             % Sampling point according to fitness evaluation for ploting the convergence curve

for r=1:runs
    NFEs = 0;
    show =0;
    hisx = [];
    hisf = [];
   
    % 1. Initialization procedure
    sam=L_bound+(U_bound-L_bound).*lhsdesign(initial_sample_size,Dim);         % Initial LHS samples
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
%% Main loop
% update P and E
P = hisx(1:NP,:); E = hisf(1:NP); % 每一代都把最好的个体保留下来了
while NFEs < Max_NFEs
    G = G + 1; disp(['iterations: ' num2str(G)]); 
    if Strategy == 1                           
        % Strategy 1: Surrogate screening
        Strategy1 = Strategy1 + 1; 
        
        % 1.1.Evolutionary operation
        LB = repmat((L_bound),length(E),1); 
        UB = repmat((U_bound),length(E),1);
        NP_last = NP; 
        NP_next = NP;
        
        U = DEoperating(P,NP_last,NP_next,Dim,hisx,F,CR,UB,LB);
        
        % 1.2.Surrogate screening
        % build global surrogate model
        [ghf,id]=sort(hisf);                            % sort history data 
        gs=length(ghf(1:end));                          % global sample number  
        
        ghx=hisx(id(1:gs),:);  ghf=ghf(1:gs);           % model_1 training samples
        ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
        spr=max(max(ghxd))/(Dim*gs)^(1/Dim);
        if spr <= 0
            spr = 1
        end
        net=newrbe(ghx',ghf',spr);        % newrbe - build RBF Neural Networks
        modelFUN=@(x) sim(net,x');
        
        % surrogate model evaluate
        fitnessModel=modelFUN(U); 
        
        % sort and screening
        [fit,sidx]=sort(fitnessModel);         % sort point based on fitness, get point indexs
        sam=U(sidx,:);                         % sorted sample
        candidate_position = sam(1,:); 
        fit0 = E(sidx(1));
    elseif  Strategy == 2
        % Strategy 2: Surrogate sampling
        Strategy2 = Strategy2 + 1;
        
        % 2.1 Build local surrogate model 
        [lhf, id]=sort(hisf);               
        lhx=hisx(id(1:ls),:); lhf = hisf(id(1:ls));
        lhxd=real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        
        spr=max(max(lhxd))/(Dim*ls)^(1/Dim);
        if spr <= 0
            spr = 1
        end
        net=newrbe(lhx',lhf',spr);        % newrbe
        LocalModelFUN=@(x) sim(net,x');
        
        % 2.2 Find a optimum of surrogate model by optimizer
        maxgen=10000; minerror=1e-20;
        [candidate_position,~] = DE(Dim, maxgen, LocalModelFUN, minerror, lhx); 
        fit0 = lhf(1);
    end
    
    % Judge Repeat Sample
    [~,ih,~]=intersect(hisx,candidate_position,'rows'); % 
    if isempty(ih)~=1
        disp(['Sample Repeat and Delete it']);
        continue;
    end
    
    % evaluate candidate
    candidate_fit=FUN(candidate_position);
    NFEs = NFEs + 1; 
    if show 
        disp([' Strategy=' num2str(Strategy) ]);
        disp([' candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
    end
    
    if Strategy == 1
        % save candidate into dataset, and sort dataset
        if candidate_fit <= E(sidx(1))
            P(sidx(1),:) = candidate_position;
            E(sidx(1)) = candidate_fit;
        end
        if candidate_fit < hisf(1)
            P = [P; candidate_position];
            E = [E; candidate_fit];
            NP = NP + 1;
        end
    end
    % save candidate into dataset, and sort dataset
    hisx=[hisx; candidate_position];  hisf=[hisf; candidate_fit];   % update history database 
    [hisf,sidx]=sort(hisf);                                         % sort point based on fitness, get point indexs
    hisx=hisx(sidx,:); 
    
    StrategyChange = 1;
    % update BestEvaluation and BestPoint
    if  candidate_fit <= fit0
        BP = candidate_position;
        BE = candidate_fit;
        disp(['Best Cost(global search) = ' num2str(BE) ' NFE=' num2str(NFEs)]);
        StrategyChange = 0;
    end
    % 策略转换
    if StrategyChange
        Strategy = 3-Strategy;
    end
    
    % update CE for plotting
    CE(NFEs,:)=[NFEs,candidate_fit];
    if mod (NFEs,sn1)==0
        cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
    end
    
end
% output
bestFitness=BE;
bestP=BP;

fittest = bestP;
clfittest =-bestFitness;
if clfittest==cloriginal
    fittest=p0;
end
% %plotting the original airfoil vs. the evolved (optimized)
% plotairfoil(fittest,'k')
% axis('equal')
% hold on
% plotairfoil(p0,'r')
% legend('Optimized','original')
% xlabel('X/C')
% ylabel('Y/C')
% title('Airfoil shape')
% fprintf('Best fitness: %e\n',min(hisf));    
end