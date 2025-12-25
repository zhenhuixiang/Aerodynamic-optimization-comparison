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

%%
model_num = 4;  
t1 = clock; 

% prepare good fitness data and center point
[DB_f,id]=sort(hf); 
DB_x = hx(id,:);
DB = [DB_x DB_f];
DS = length(hf);
GDS = floor(DS*0.2); % top 20% in the DB to form GD
if GDS < 10 
    GDS = 10;
end
if GDS > DS 
    GDS = DS;
end
GD_x = DB_x(1:GDS,:);
center = mean(GD_x,1);

% train and test dataset
E = zeros(1,model_num); % Model Error Criterion
R = zeros(1,model_num); % Distance Deviation Criterion
test_subset = 10; % number of subset
for j = 1:test_subset
    % dataset divide
    index = randperm(GDS,floor(GDS/2));
    test_hf = DB_f(index); test_hx = DB_x(index,:);
    train_hf = DB_f;  train_hx = DB_x;
    train_hf(index) = [];  train_hx(index,:) = [];
    DS_train = length(train_hf);

    % RBF parameters
    ghxd=real(sqrt(train_hx.^2*ones(size(train_hx'))+ones(size(train_hx))*(train_hx').^2-2*train_hx*(train_hx')));
    spr=max(max(ghxd))/(dims*DS_train)^(1/dims);

    % train models
    for i = 1:model_num
        % RBF network
        h = 4*(i-1);
        spr= spr * 2^h;
        net=newrbe(train_hx',train_hf',spr);
        RBF_FUN=@(x) sim(net,x');
        % test error
        error = abs(test_hf - RBF_FUN(test_hx)');
        E(i) = E(i) + sum(error);
    end
end
E = E/test_subset;

% build models including all samples for predict
predict = zeros(1,model_num);
predict_pos = [];
for i = 1:model_num
    % RBF network use all data
    ghxd=real(sqrt(hx.^2*ones(size(hx'))+ones(size(hx))*(hx').^2-2*hx*(hx')));
    spr=max(max(ghxd))/(dims*DS)^(1/dims);
    h = 4*(i-1);
    spr= spr * 2^h;
    net=newrbe(hx',hf',spr);
    RBF_FUN=['RBF_FUN',num2str(i)];
    eval([RBF_FUN,'=@(x) sim(net,x'');']);
    eval(['RBF_FUN=',RBF_FUN,';']);

    % model predict
    maxgen = 30000+300*dims;
    minerror = 1e-10;
    RBF_FUN_Arr{i} = RBF_FUN;
    [best_pos,bestever] = SLPSO(dims, maxgen, RBF_FUN, minerror, GD_x);

    % obtain real fitness to verify effect of this algorithm
    predict_pos(i,:) = best_pos;
    predict(i) = FUN(best_pos);
end


T = 10;
for k = 1:T
    R_temp = pdist2(DB_x(k,:), predict_pos);
    R = R_temp + R;
end

R_N = mapminmax(R);
E_N = mapminmax(E);
loss = (R_N + E_N)/2;
[~, index]=min(loss);

model_selected = index; % selected model index
pos = predict_pos(index,:); % predicted solution
pred = predict(index); % fitness of solution


totaltime = etime(clock,t1);

%% output
bestFitness=pred;
bestP=pos;

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