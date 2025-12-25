% load offline data
load('DDSOP1.mat')
load('DDSOP2.mat')
load('DDSOP3.mat')

% the boundary of decision variable
LB = [ 0.0001, 0.0001, 0.0001, 0.0001];
UB = [ 100, 100, 300, 100];
Dim = length(LB);

DDSOP1(:, end) = -DDSOP1(:, end);
[time,P,gbest ] = DDEA_SE(Dim,DDSOP1,UB,LB);
P