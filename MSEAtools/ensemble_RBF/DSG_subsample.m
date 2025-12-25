function [ W,B,C,Sp] = DSG_subsample( L,c,nc,T,l)
% Usage: [ W,B,C,S] = RBF_EnsembleUN( L,c,nc,T)
%Build RBF Model Pool
% Input:
% L             - Offline Data with c Decision Variables and Exact Objective Value
% c             - Number of Decision Variables
% nc            - Number of neurons of RBF models
% T             - %Number of RBF models
%
% Output: 
% W             - Weights of RBF Models
% B             - Bais of RBF Models
% C             - Centers of RBF Models
% S             - Widths of RBF models
%
    %%%%    Authors:    Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty
    %%%%    University of Surrey, UK and Taiyuan University of Science and Technology, China.
    %%%%    EMAIL:      wanghanding.patch@gmail.com
    %%%%    WEBSITE:    https://sites.google.com/site/handingwanghomepage
    %%%%    DATE:       May 2018
%------------------------------------------------------------------------
%This code is part of the program that produces the results in the following paper:

%Handing Wang, Yaochu Jin, Chaoli Sun, John Doherty, Offline data-driven evolutionary optimization using selective surrogate ensembles, IEEE Transactions on Evolutionary Computation, Accepted.

%You are free to use it for non-commercial purposes. However, we do not offer any forms of guanrantee or warranty associated with the code. We would appreciate your acknowledgement.
%------------------------------------------------------------------------


S_size = floor(size(L,1)/10);
nc = size(L,1) - S_size;

[ W2,B2,Centers,Spreads ] = RBF( L(:,1:c),L(:,c+1),nc);
Y = RBF_predictor( W2,B2,Centers,Spreads,L(:,1:c));
diff = abs(Y - L(:,c+1));
[diff_sort,Idx] = sort(diff,'descend');
S = L(Idx(1:S_size),1:c+1);

W=zeros(T,nc);
B=zeros(1,T);
C=zeros(c,nc,T);
Sp=zeros(nc,T);
DTS = zeros(size(L,1)-S_size,c+1,T);

% W(1,:)=W2;
% B(1)=B2;
% C(:,:,1)=Centers;
% Sp(:,1)=Spreads;

D = size(L,2)-1;
k = norm(ones(1,D)); %生成的最长的随机数的长度

    %Build Model Pool
    for j=1:T
        del = randperm(size(L,1),S_size);
        K = L;
        K(del,:) = [];
        DTS(:,:,j) = [K];
        [ W2,B2,Centers,Spreads ] = RBF( DTS(:,1:c,j),DTS(:,c+1,j),nc);
        W(j,:)=W2;
        B(j)=B2;
        C(:,:,j)=Centers;
        Sp(:,j)=Spreads;
    end

end

