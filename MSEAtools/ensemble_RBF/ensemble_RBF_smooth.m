function [time, predict] =ensemble_RBF_smooth(c,L,bu,bd)
tic;
rand('state',sum(100*clock));
nc=size(L,1);%Number of neurons of RBF models
T=10;%Number of RBF models for DSG
Q=10;%Ensemble size
D = size(L,2)-1;

%Build Model Pool
l = norm(bu - bd)/sqrt(D) * 1e-6 ;
[~,idx] = min(L(:,c+1));
Xb = L(idx,:);

% [ W,B,C,S] = DSG_supsample( L,c,nc,T,l);
[ W,B,C,S] = DSG_subsample_smooth( L,c,nc,T,l);
[ I] = SE(W,B,C,S,Xb,c,Q); % selected RBF Models
predict = @(x) mean(feval("RBF_Ensemble_predictor",W(I,:),B(I),C(:,:,I),S(:,I),x,c),2);
% predict = @(x) feval("RBF_Ensemble_predictor",W(I,:),B(I),C(:,:,I),S(:,I),x,c);
toc;
time=toc;

end

