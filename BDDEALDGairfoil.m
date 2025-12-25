function [cloriginal,clfittest,fittest]=BDDEALDGairfoil(genNo,p0,range,uinf,AOA,Npanel) 
    
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
UB = max(p0-range,p0+range);
LB = min(p0-range,p0+range);
Dim = 11;
c = Dim;

dims = Dim;
initial_sample_size = 11*dims;
sam=repmat(LB,initial_sample_size,1)+(repmat(UB,initial_sample_size,1)-repmat(LB,initial_sample_size,1)).*lhsdesign(initial_sample_size,dims);
for i = 1:initial_sample_size
    fitness(i)=FUN(sam(i,:));
end
hx=sam; hf=fitness';
DATA = [hx, hf];
L = DATA;
bu = UB;
bd = LB;

rand('state',sum(100*clock));
nc=c;%Number of neurons of RBF models
T=100;%Number of RBF models for DSG
D = size(L,2)-1;

%Build Model Pool
l = norm(bu - bd)/sqrt(D) * 1e-6 ;
[~,idx] = min(L(:,c+1));
Xb = L(idx,:);

[ W,B,C,S] = LDG_BS( L,c,nc,T,l);


%-------------------------------------------
gmax=100;
pc=1;%Crossover Probability
pm=1/c;%Mutation Probability
n=100;%Population Sixe
% Online Optimization-------------------------------------------
tic;
POP = initialize_pop(n,c,bu,bd);

%RBF Predictors 
Y= RBF_Ensemble_predictor( W,B,C,S,POP,c );
POP=[POP,Y];
g=1;
gbest=[];
while g<=gmax
    if g~=1
        POP=POP(:,1:c);
        Y= RBF_Ensemble_predictor(W,B,C,S,POP,c );
        POP=[POP,Y];
    end
    %Variations    
    NPOP1=SBX( POP,bu,bd,pc,n );
    [ Y ] = RBF_Ensemble_predictor( W,B,C,S,NPOP1,c );
    NPOP1=[NPOP1,Y];
    NPOP2=mutation(POP,bu,bd,pm,n);
    [ Y ] = RBF_Ensemble_predictor( W,B,C,S,NPOP2,c );
    NPOP2=[NPOP2,Y];
    POP=[POP;NPOP1;NPOP2];
    %Model Combination
    YAVE=mean(POP(:,c+1:end),2);
    [A,Is]=sort(YAVE);
    POP=[POP(Is(1:n),1:c)];
    g=g+1;
    P= POP(1,1:c);
    gbest=[gbest;P];
    
end

predict = @(x) feval("RBF_Ensemble_predictor",W,B,C,S,x,c);

toc;
time=toc;

%% output
bestFitness=FUN(P);
bestP=P;

fittest = bestP;
clfittest =-bestFitness;
if clfittest==cloriginal
    fittest=p0;
end
%% plotting the original airfoil vs. the evolved (optimized)
% plotairfoil(fittest,'k')
% axis('equal')
% hold on
% plotairfoil(p0,'r')
% legend('Optimized','original')
% xlabel('X/C')
% ylabel('Y/C')
% title('Airfoil shape')
% fprintf('Best fitness: %e\n',clfittest);    
end