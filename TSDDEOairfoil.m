function [cloriginal,clfittest,fittest,gfs]=TSDDEOairfoil(genNo,p0,range,uinf,AOA,Npanel) 
% 参数调整
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
U_bound = max(p0-range,p0+range);
L_bound = min(p0-range,p0+range);
NP = 50;
Dim = 11;
MaxFEs = 500;
initial_sample_size=50;    % < 100 dimension

D = Dim;
LB = L_bound;
UB = U_bound;


% ---------------  Initialization --------------------
NFEs=0;
ps=50;                          % population size
CE=zeros(MaxFEs,2);             % achive the exact fitness
sn1=1;                          % 保存数据时的数据压缩程度
gen = 1;
show = 0;

gfs=zeros(1,fix(MaxFEs/sn1));   % sampling point according to fitness evaluation for ploting the convergence curve

% -------------- initial LHS samples --------------------
sam=repmat(LB,initial_sample_size,1)+(repmat(UB,initial_sample_size,1)-repmat(LB,initial_sample_size,1)).*lhsdesign(initial_sample_size,D);
fit=zeros(1,initial_sample_size);
for i=1:initial_sample_size
    x = sam(i,:);
    fitness = FUN(x);
    fit(i)=fitness;  
    NFEs=NFEs+1;
    CE(NFEs,:)=[NFEs,fit(i)];
    if mod (NFEs,sn1)==0
        cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
    end
end
% build database
hx=sam; hf=fit;                                            % archive all the exact evaluated samples
[~,sidx]=sort(hf);                                         % 对个体按适应值由小到大进行排序
hx=hx(sidx,:);  hf=fit(sidx);                              % 排序后的样本

% -------------- 建模训练样本数 -----------------------
gs = initial_sample_size;                   
ls = 50;

% -------------- first stage parameters ---------------------
cc=[2.05 2.05];                      % acceleration constants  
iwt=0.7298;                          

VRmax = UB;
VRmin = LB;
mv=0.5*(VRmax-VRmin);
VRmin_ = VRmin;
VRmax_ = VRmax;
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);            %---位置上界,防溢出
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;                          %---速度上界
vel=Vmin+2.*Vmax.*rand(ps,D);        % initialize the velocity of the particles
pos=hx(1:ps,:); e=hf(1:ps);                       % initialize position and fitness
pbest=pos;  pbestval=e;              % initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);              % initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);
besty=1e200;
bestp=zeros(1,D);

% -------------- secend stage 参数 --------------------- 
F = 0.5;
CR = 0.5;

%% Main loop
while NFEs < MaxFEs
    % Sample training samples
    [ghf,id]=sort(hf);               % 使用历史数据库中最优的若干个体构造全局代理模型
    if NFEs < 200 
        % R B F network
        ghf=ghf(1:gs);     ghx=hx(id(1:gs),:);
        ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
        spr=max(max(ghxd))/(D*gs)^(1/D);
        net=newrbe(ghx',ghf,spr);        
        modelFUN=@(x) sim(net,x');

        % Record the old best position and fitness
        besty_old=besty;
        bestp_old=bestp;
        
        % Find optimum of surrogate model by SLPSO
        maxgen=50*D; minerror=1e-8;
        [bestp,~] = SLPSO(D,maxgen,modelFUN,minerror,ghx);    % 利用SLPSO寻找全局模型的最优解
        [~,ih,ip]=intersect(hx,bestp,'rows');
        if isempty(ih)==1                               % 如果得到的不是重复数据
            besty=FUN(bestp); 
            hx=[hx;bestp];  hf=[hf,besty];              % update history database 
            % Exact evaluate the model optimum
            NFEs = NFEs + 1;
            if NFEs <= MaxFEs
                CE(NFEs,:)=[NFEs,besty];
                if mod (NFEs,sn1)==0
                    cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
                end
            end
            
            % Update the RBF model                                        
            if ghf(end) > besty                             % 选取最优样本集建立全局模型
                [ghf,id]=sort(hf);                          % 使用历史数据库中最优的若干个体构造全局代理模型
                ghf=ghf(1:gs); ghx=hx(id(1:gs),:);
                % R B F network
                ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));  
                spr=max(max(ghxd))/(D*gs)^(1/D);
                net=newrbe(ghx',ghf,spr);       % newrbe
                modelFUN=@(x) sim(net,x'); 
            end
        else
            besty = hf(ih,:);
        end
        
        % Record the new best position and fitness of surrogate
        besty_new=besty;
        bestp_new=bestp;
        % Update model optimum
        if besty_new < besty_old      
            besty=besty_new;    
            bestp=bestp_new;
            bestprep=repmat(bestp_new,ps,1);
        else
            besty=besty_old;
            bestp=bestp_old;
            bestprep=repmat(bestp_old,ps,1);
        end
        
        if besty < gbestval                             % 代理模型的近似最优小于当前代PSO种群的gbest,使用besty引导粒子向模型最优位置飞行
            [~,ip,~]=intersect(pbest,gbest,'rows');
            pbest(ip,:) =bestp;
            pbestval(ip)=besty;                         % 更新gbest对应的pbest
            gbestrep=bestprep;
            aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        else
            aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        end

        vel=iwt.*(vel+aa);                              % 带收缩因子的PSO算法
        vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
        vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
        pos=pos+vel;
        pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
            +(pos<VRmin).*(VRmin+0.25.*(VRmax-VRmin).*rand(ps,D))+(pos>VRmax).*(VRmax-0.25.*(VRmax-VRmin).*rand(ps,D));

       %% Fitness estimation of new population
        e=modelFUN(pos)';
        
        candidx=find(e' < pbestval);
        pos_trmem=pos(candidx, :);                 % e-pbest strategy    
        [~,ih,ip]=intersect(hx,pos_trmem,'rows');
        if ~isempty(ip)==1        
            pos_trmem(ip,:)=[];
            e(candidx(ip))=hf(ih);
            candidx(ip)=[];
        end
        
        % Exact evaluate the prescreened candidates
        ssk=size(pos_trmem,1);
        e_trmem=zeros(1,ssk);        
        for k=1:ssk
            e_trmem(k)=FUN(pos_trmem(k,:));
            NFEs = NFEs + 1;
            if NFEs <= MaxFEs 
                CE(NFEs,:)=[NFEs,e_trmem(k)];
                if mod (NFEs,sn1)==0
                    cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
                end
            end
            hx=[hx;pos_trmem(k,:)];hf=[hf,e_trmem(k)];        % update history database 
            kp=candidx(k);
            if e_trmem(k)<pbestval(kp)                        % update pbest
                pbest(kp,:)=pos_trmem(k,:);
                pbestval(kp)=e_trmem(k); 
            end
        end  

        % Update gbest
        [gbestval,tmp]=min(pbestval); 
        gbest=pbest(tmp,:);
        gbestrep=repmat(gbest,ps,1);                          % update the gbest    
        bestfit=min([gbestval,besty]);
        fprintf(1,'Iteration: %d,  No.evaluation: %d,  Best: %e,  No.prescreen data: %d\n',gen,NFEs,bestfit,ssk);
        gen = gen + 1;
        
    else
        
        disp(['iterations: ' num2str(gen)]); 
        % 1.1.Evolutionary operation
        NP_last = 50; 
        NP_next = 50;
        
        VRmin2=repmat(VRmin_,NP_last,1);
        VRmax2=repmat(VRmax_,NP_last,1);            %---位置上界,防溢出
        
        % update P and E
        P = hx(1:NP_last,:); E = hf(1:NP_last); % 每一代都把最好的个体保留下来了
        
        U = DEoperating(P,NP_last,NP_last,D,hx,F,CR,VRmax2,VRmin2);
        
        % 1.2.Surrogate screening
        % build global surrogate model
        [ghf,id]=sort(hf);                            % sort history data 
        ghx=hx(id(1:gs),:);  ghf=ghf(1:gs);           % model_1 training samples
        ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
        spr=max(max(ghxd))/(D*gs)^(1/D);
        net=newrbe(ghx',ghf,spr);        % newrbe - build RBF Neural Networks
        modelFUN=@(x) sim(net,x');
        
        % surrogate model evaluate
        fitnessModel=modelFUN(U); 
        
        % sort and screening
        [fit,sidx]=sort(fitnessModel);         % sort point based on fitness, get point indexs
        sam=U(sidx,:);                         % sorted sample
        candidate_position = sam(1,:); 
        
         % Judge Repeat Sample
        [~,ih,~]=intersect(hx,candidate_position,'rows'); % 
        if isempty(ih)~=1
            disp(['Sample Repeat and Delete it']);
            continue;
        end

        % evaluate candidate
        candidate_fit=FUN(candidate_position);
        NFEs = NFEs + 1; 
        if show 
            disp([' candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
        end

        % save candidate into dataset, and sort dataset
        hx=[hx; candidate_position];  hf=[hf, candidate_fit];   % update history database 
        [hf,sidx]=sort(hf);                                         % sort point based on fitness, get point indexs
        hx=hx(sidx,:); 

        % update BestEvaluation and BestPoint
        if  candidate_fit <= hf(1) 
            BP = candidate_position;
            BE = candidate_fit;
            disp(['Best Cost(1) = ' num2str(BE) ' NFE=' num2str(NFEs)]);
        end

        % update CE for plotting
        CE(NFEs,:)=[NFEs,candidate_fit];
        if mod (NFEs,sn1)==0
            cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
        end
        
        [lhf, id]=sort(hf);               
        lhx=hx(id(1:ls),:); lhf = hf(id(1:ls));
        lhxd=real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
        
        spr=max(max(lhxd))/(D*ls)^(1/D);
        net=newrbe(lhx',lhf,spr);        % newrbe
        LocalModelFUN=@(x) sim(net,x');
        
        % 2.2 Find a optimum of surrogate model by optimizer
         maxgen=50*D; minerror=1e-8;
        [candidate_position,~] = DE(D, maxgen, LocalModelFUN, minerror, lhx); 
        
         % Judge Repeat Sample
        [~,ih,~]=intersect(hx,candidate_position,'rows'); % 
        if isempty(ih)~=1
            disp(['Sample Repeat and Delete it']);
            continue;
        end

        % evaluate candidate
        candidate_fit=FUN(candidate_position);
        NFEs = NFEs + 1; 
        if show 
            disp([' candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
        end

        % save candidate into dataset, and sort dataset
        hx=[hx; candidate_position];  hf=[hf, candidate_fit];   % update history database 
        [hf,sidx]=sort(hf);                                         % sort point based on fitness, get point indexs
        hx=hx(sidx,:); 

        % update BestEvaluation and BestPoint
        if  candidate_fit <= hf(1) 
            BP = candidate_position;
            BE = candidate_fit;
            disp(['Best Cost(2) = ' num2str(BE) ' NFE=' num2str(NFEs)]);
        end

        % update CE for plotting
        CE(NFEs,:)=[NFEs,candidate_fit];
        if mod (NFEs,sn1)==0
            cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
        end
        
        if rand < 0.2
            ls2 = gs;
            [lhf, id]=sort(hf);               
            lhx=hx(id(1:ls2),:); lhf = hf(id(1:ls2));  % lhx 全局代理的样本自变量，lhf 全局代理的样本因变量
            lhxd=real(sqrt(lhx.^2*ones(size(lhx'))+ones(size(lhx))*(lhx').^2-2*lhx*(lhx')));
            spr=max(max(lhxd))/(D*ls2)^(1/D); 
            if spr <= 0
                spr = 1;
            end
            net=newrbe(lhx',lhf,spr);        % newrbe
            LocalModelFUN=@(x) sim(net,x');
            [candidate_position] = full_crossover(LocalModelFUN,lhx); 

            % Judge Repeat Sample
            [~,ih,~]=intersect(hx,candidate_position,'rows'); % 
            if isempty(ih)~=1
                disp(['Sample Repeat and Delete it']);
                continue;
            end

            % evaluate candidate
            candidate_fit=FUN(candidate_position);
            NFEs = NFEs + 1; 
            if show 
                disp([' candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
            end

            % save candidate into dataset, and sort dataset
            hx=[hx; candidate_position];  hf=[hf, candidate_fit];       % update history database 
            [hf,sidx]=sort(hf);                                         % sort point based on fitness, get point indexs
            hx=hx(sidx,:); 

            % update BestEvaluation and BestPoint
            if  candidate_fit <= hf(1) 
                BP = candidate_position;
                BE = candidate_fit;
                disp(['Best Cost(3) = ' num2str(BE) ' NFE=' num2str(NFEs)]);
            end

            % update CE for plotting
            CE(NFEs,:)=[NFEs,candidate_fit];
            if mod (NFEs,sn1)==0
                cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
            end
        end
        gen = gen + 1;
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
gfs = gfs(1:500);
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
end
