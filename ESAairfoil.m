%% Core code of ESA
function [cloriginal,clfittest,fittest,gfs]=ESAairfoil(genNo,p0,range,uinf,AOA,Npanel) 
% 参数调整
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
U_bound = max(p0-range,p0+range);
L_bound = min(p0-range,p0+range);
NP = 50;
initial_sample_size=50; 
Dim = 11;

warning off
format short;
format compact; 

%% Parameters setting
NFE=0;
MaxNFE = 500;
CE=zeros(MaxNFE,2);             % achive the exact fitness
sn1=1;                          % The degree of data compression when saving data
gen = 1;
if Dim < 100
    initial_sample_size=100;    % < 100 dimension
elseif Dim >= 100
    initial_sample_size=150;    % >= 100 dimension                       
end
gfs=zeros(1,fix(MaxNFE/sn1));   % sampling point according to fitness evaluation for ploting the convergence curve

% -------------- initial LHS samples --------------------
sam=repmat(L_bound,initial_sample_size,1)+(repmat(U_bound,initial_sample_size,1)-repmat(L_bound,initial_sample_size,1)).*lhsdesign(initial_sample_size,Dim);
fit=zeros(1,initial_sample_size);
for i=1:initial_sample_size
    x = sam(i,:);
    fitness = FUN(x);
    fit(i)=fitness;  
    NFE=NFE+1;
    CE(NFE,:)=[NFE,fit(i)];
    if mod (NFE,sn1)==0
        cs1=NFE/sn1; gfs(1,cs1)=min(CE(1:NFE,2));
    end
end
% build database
hx=sam; hf=fit;                                            % archive all the exact evaluated samples
[~,sidx]=sort(hf);                                         % The individuals were sorted according to their fitness values
hx=hx(sidx,:);  hf=fit(sidx);                              % Sorted DB
% ------------------------------------------------------

% DE parameter
NP = 50;
F=0.5;
CR=0.9;

% number of local search samples
ls = 25 + Dim;
if ls > 60
    ls = 60;
end
ls2 = 100;

% number of calling Strategys 
Strategy1 = 0;
Strategy2 = 0;
Strategy3 = 0;
Strategy4 = 0;

% RL setting
State =  randperm(8,1); 
Action =  randperm(4,1); 
alp=0.1;                    
gamma = 0.9;                
Q_Agent=zeros(State,Action);
Q_Agent = [ 0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;
            0.25 0.25 0.25 0.25;];

% show evaluated candidate each generation
show = 0; 

%% Main loop
while NFE <= MaxNFE
    disp(['NFE: ' num2str(NFE)]); 
    State_old = State;
    state_offset = 0;
    R = 0; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update state and action 
    Qvalue1=Q_Agent(State,:); 
    temp=exp(Qvalue1);
    ratio=cumsum(temp)/sum(temp); 
    jtemp=find(rand(1)<ratio);
    Action=jtemp(1);    
    Strategy = Action; % action
    
    % log data 
    log_Q_Agent{NFE} = Q_Agent;
    log_State{NFE} = State;
    log_ratio{NFE} = ratio;
    log_Action{NFE} = Action;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update state and action 
    
    % Actions
    if Strategy == 1           % Strategy 1: Surrogate screening     
        Strategy1 = Strategy1 + 1; 
        LB = repmat((L_bound),NP,1); 
        UB = repmat((U_bound),NP,1);
        P = hx(1:NP,:); E = hf(1:NP); % update P and E
        U = DEoperating(P,NP,NP,Dim,hx,F,CR,UB,LB);
        [ghf,id]=sort(hf);                            % sort history data 
        gs=length(ghf(1:end));                        % global sample number  
        if gs > 300
            gs = 300;
        end
        ghx=hx(id(1:gs),:);  ghf=ghf(1:gs);           % model_1 training samples
        ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
        spr=max(max(ghxd))/(Dim*gs)^(1/Dim);
        net=newrbe(ghx',ghf,spr);                       % newrbe - build RBF Neural Networks
        modelFUN=@(x) sim(net,x');                      % build global surrogate model
        fitnessModel=modelFUN(U);                       % surrogate model evaluate
        [fit,sidx]=sort(fitnessModel);                  % sort point based on fitness, get point indexs
        sam=U(sidx,:);                                  % sorted sample
        candidate_position = sam(1,:);                  % screening a candidate
    elseif  Strategy == 2       % Strategy 2: Surrogate sampling
        Strategy2 = Strategy2 + 1;
        [~, id]=sort(hf);
        lhx=hx(id(1:ls),:); lhf = hf(id(1:ls));
        lhxd=real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        spr=max(max(lhxd))/(Dim*ls)^(1/Dim);
        net=newrbe(lhx',lhf,spr);        % newrbe
        LocalModelFUN=@(x) sim(net,x');  % build local surrogate model 
        a1 = 100;
        Max_NFE=a1*Dim+30000; minerror=1e-20;
        [candidate_position,~] = DE(Dim, Max_NFE, LocalModelFUN, minerror, lhx); % find a optimum of surrogate model by optimizer
    elseif  Strategy == 3      % Strategy 3: FC 
        Strategy3 = Strategy3 + 1;
        [lhf, id]=sort(hf);               
        lhx=hx(id(1:ls2),:); lhf = hf(id(1:ls2));  % lhx 全局代理的样本自变量，lhf 全局代理的样本因变量
        lhxd=real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        spr=max(max(lhxd))/(Dim*ls2)^(1/Dim); 
        net=newrbe(lhx',lhf,spr);        % newrbe
        LocalModelFUN=@(x) sim(net,x');
        [candidate_position] = full_crossover(LocalModelFUN,lhx); 
    elseif  Strategy == 4      % Strategy 4: STR 
        Strategy4 = Strategy4 + 1;
        [lhf, id]=sort(hf);  
        idx = id(1:ls2);
        lhx=hx(idx,:); lhf = hf(idx);  % lhx 全局代理的样本自变量，lhf 全局代理的样本因变量
        [newdata_x, newdata_f] = STR(FUN,lhx,lhf); 
        % record
        num_c = size(newdata_f,1);
        best_candidate_fit = inf;
        best_candidate_position = [];
        for a = 1:num_c
            NFE = NFE + 1; 
            candidate_position = newdata_x(a,:);
            candidate_fit = newdata_f(a);
            if show 
                disp(['Strategy=' num2str(Strategy) ]);
                disp(['candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
            end
            % save candidate into dataset, and sort dataset
            hx=[hx; candidate_position];  hf=[hf, candidate_fit];   % update htory database 
            [hf,idx]=sort(hf);                                         % sort point based on fitness, get point indexs
            hx=hx(idx,:); 
            % update BestEvaluation and BestPoint
            if  candidate_fit <= hf(1) 
                state_offset = 1; R=1;
                BP = candidate_position;
                BE = candidate_fit;
                disp(['Best Cost(Strategy ' num2str(Strategy)  ') = ' num2str(BE) ' NFE=' num2str(NFE)]);
            end
            % update CE for plotting
            CE(NFE,:)=[NFE,candidate_fit];
            if mod (NFE,sn1)==0
                cs1=NFE/sn1; gfs(1,cs1)=min(CE(1:NFE,2));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
        State=2*Strategy+state_offset-1; % next state
        temp=max(Q_Agent(State,:));
        Q_Agent(State_old,Action)=(1-alp)*Q_Agent(State_old, Action)+alp*(R+gamma*temp); 
        State_old = State;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
    end
    
    % judge Repeat Sample 
    [~,ih,~]=intersect(hx,candidate_position,'rows'); 
    if isempty(ih)~=1
        disp(['Sample Repeat and Delete it']);
        continue;
    end
    % evaluate candidate
    candidate_fit=FUN(candidate_position);
    NFE = NFE + 1; 
    if show 
        disp(['Strategy=' num2str(Strategy) ]);
        disp(['candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
    end
    % save candidate into dataset, and sort dataset
    hx=[hx; candidate_position];  hf=[hf, candidate_fit];   % update htory database 
    [hf,sidx]=sort(hf);                                         % sort point based on fitness, get point indexs
    hx=hx(sidx,:);  
    % update BestEvaluation and BestPoint
    if  candidate_fit <= hf(1) 
        state_offset = 1; R=1;
        BP = candidate_position;
        BE = candidate_fit;
        disp(['Best Cost(Strategy ' num2str(Strategy)  ') = ' num2str(BE) ' NFE=' num2str(NFE)]);
    end
    % update CE for plotting
    CE(NFE,:)=[NFE,candidate_fit];
    if mod (NFE,sn1)==0
        cs1=NFE/sn1; gfs(1,cs1)=min(CE(1:NFE,2));
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
    State=2*Strategy+state_offset-1; % next state 
    temp=max(Q_Agent(State,:));
    Q_Agent(State_old,Action)=(1-alp)*Q_Agent(State_old, Action)+alp*(R+gamma*temp); 
    State_old = State;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update agent
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
% fprintf('Best fitness: %e\n',min(hf));    