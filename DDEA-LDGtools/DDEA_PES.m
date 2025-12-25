function [time,P,gbest, predict] =DDEA_PES(c,L,bu,bd)

rand('state',sum(100*clock));
nc=c;%Number of neurons of RBF models
T=200;%Number of RBF models for DSG
Q=100;%Ensemble size
D = size(L,2)-1;

%Build Model Pool
l = norm(bu - bd)/sqrt(D) * 1e-6 ;
[~,idx] = min(L(:,c+1));
Xb = L(idx,:);

[ W,B,C,S] = DSG( L,c,nc,T,l);
[ I] = SE(W,B,C,S,Xb,c,Q); % selected RBF Models


%-------------------------------------------
gmax=100;
pc=1;%Crossover Probability
pm=1/c;%Mutation Probability
n=100;%Population Sixe
% Online Optimization-------------------------------------------
tic;
POP = initialize_pop(n,c,bu,bd);

%RBF Predictors 
Y= RBF_Ensemble_predictor( W(I,:),B(I),C(:,:,I),S(:,I),POP,c );
POP=[POP,Y];
g=1;
gbest=[];
while g<=gmax
    if g~=1
        POP=POP(:,1:c);
        Y= RBF_Ensemble_predictor(W(I,:),B(I),C(:,:,I),S(:,I),POP,c );
        POP=[POP,Y];
    end
    %Variations    
    NPOP1=SBX( POP,bu,bd,pc,n );
    [ Y ] = RBF_Ensemble_predictor( W(I,:),B(I),C(:,:,I),S(:,I),NPOP1,c );
    NPOP1=[NPOP1,Y];
    NPOP2=mutation(POP,bu,bd,pm,n);
    [ Y ] = RBF_Ensemble_predictor( W(I,:),B(I),C(:,:,I),S(:,I),NPOP2,c );
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

predict = @(x) mean(feval("RBF_Ensemble_predictor",W(I,:),B(I),C(:,:,I),S(:,I),x,c),2);

toc;
time=toc;

end

