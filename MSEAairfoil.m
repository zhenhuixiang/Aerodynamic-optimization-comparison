function [cloriginal,clfittest,fittest]=MSEAairfoil(genNo,p0,range,uinf,AOA,Npanel) 
    
addpath(genpath(pwd));

[cloriginal,~]=solver(p0,uinf,AOA,Npanel);
FUN=@(x) feval("solver",x,uinf,AOA,Npanel)*(-1);  
UB = max(p0-range,p0+range);
LB = min(p0-range,p0+range);
Dim = 11;

%
p_t = 1;                           % time of dividing test and train data
p_top = 0.2;                       % test data selected from p_top data
p_r = 0.5;                         

% Êý¾Ý
initial_sample_size = 11*Dim;
sam=repmat(LB,initial_sample_size,1)+(repmat(UB,initial_sample_size,1)-repmat(LB,initial_sample_size,1)).*rand(initial_sample_size,Dim);

for i = 1:initial_sample_size
    fitness(i)=FUN(sam(i,:));
end

hx=sam; hf=fitness';
DATA = [hx, hf];                            % generated offline Data
t1 = clock;                                 % time start

% 1.2 Initialize elite data number
[sort_hf,id]=sort(hf); 
sort_hx = hx(id,:);
nd = length(hf);
num_elite_data = floor(nd*p_r);             % elite data

% 1.3 Initialize model pool number and selection criteria
model_num = 4;
test_error = zeros(1,model_num);
co_Evaluate = zeros(1,model_num);
sort_error = zeros(1,model_num);
R = zeros(1,model_num);

% 2.Selection criteria [test_error][sort_error]
test_subset_num = p_t;
for j2 = 1:test_subset_num
    % 2.1 Divided dataset to train and test dataset
    topnum = floor(0.1*initial_sample_size); 
    index = randperm(floor(p_top*initial_sample_size),topnum);
    test_hf = sort_hf(index); test_hx = sort_hx(index,:); % hf rather than ghf
    train_hf = sort_hf;  train_hx = sort_hx;
    train_hf(index) = [];  train_hx(index,:) = [];

    % 2.2 build models based on train data and compute test_error and sort_error
    % RBF-L networks [model 1]
    i = 1;
    ghxd=real(sqrt(train_hx.^2*ones(size(train_hx'))+ones(size(train_hx))*(train_hx').^2-2*train_hx*(train_hx')));
    spr=max(max(ghxd))/(Dim*nd)^(1/Dim);
    net=newrbe(train_hx',train_hf',spr);
    RBF_FUN=@(x) sim(net,x');
    RBF_FUN_Arr{i} = RBF_FUN;
    error = abs(test_hf - RBF_FUN(test_hx)'); % test error
    test_error(i) = test_error(i) + sum(error);
    [~,id1] = sort(test_hf);
    [~,id2] = sort(RBF_FUN(test_hx)');
    sort_error(i) = sort_error(i)+ sum(abs(id1-id2));

    % RBF-G networks [model 2]
    i = 2;
    ghxd=real(sqrt(train_hx.^2*ones(size(train_hx'))+ones(size(train_hx))*(train_hx').^2-2*train_hx*(train_hx')));
    spr=max(max(ghxd))/(Dim*nd)^(1/Dim)*2^12;
    net=newrbe(train_hx',train_hf',spr);
    RBF_FUN=@(x) sim(net,x');
    RBF_FUN_Arr{i} = RBF_FUN;
    error = abs(test_hf - RBF_FUN(test_hx)'); % test error
    test_error(i) = test_error(i) + sum(error);
    [~,id1] = sort(test_hf);
    [~,id2] = sort(RBF_FUN(test_hx)');
    sort_error(i) = sort_error(i)+ sum(abs(id1-id2));

    % Ensemble of RBF-L and RBF-G [model 3]
    i = 3;
    RBF_FUN_Arr_ensemble = RBF_FUN_Arr([1,2]);
    ensemble = @(x)ensemble_RBF_FUN(x,RBF_FUN_Arr_ensemble);
    RBF_FUN_Arr{3} = ensemble;
    error = abs(test_hf - ensemble(test_hx)');
    test_error(i) = test_error(i) + sum(error);
    [~,id1] = sort(test_hf);
    [~,id2] = sort(ensemble(test_hx)');
    sort_error(i) = sort_error(i)+ sum(abs(id1-id2));

    % Ensemble model of many RBFs [model 4]
    i = 4;
    [ensemble_RBF_time,ensemble_model] = ensemble_RBF_smooth(Dim,[train_hx, train_hf],UB,LB);
    RBF_FUN_Arr{4} = ensemble_model;
    error = abs(test_hf - ensemble_model(test_hx));
    test_error(i) = test_error(i) + sum(error); 
    [~,id1] = sort(test_hf);
    [~,id2] = sort(ensemble_model(test_hx));
    sort_error(i) = sort_error(i)+ sum(abs(id1-id2)); 
end
test_error = test_error/test_subset_num; %%%%%%%%%%%%%%%%% Selection criteria 1 [test_error]
sort_error = sort_error/test_subset_num; %%%%%%%%%%%%%%%%% Selection criteria 2 [sort_error]

% 2.3 Predict by build models with all data [search]
% RBF network use all data
elite_data = sort_hx(1:num_elite_data,:);
predict_pos = [];
predict = [];
% Search predicted optimal based on surrogate RBF-L [model 1]
i = 1;
ghxd=real(sqrt(hx.^2*ones(size(hx'))+ones(size(hx))*(hx').^2-2*hx*(hx')));
spr=max(max(ghxd))/(Dim*nd)^(1/Dim);
net=newrbe(hx',hf',spr);
RBF_FUN=['RBF_FUN',num2str(i)];
eval([RBF_FUN,'=@(x) sim(net,x'');']);
eval(['RBF_FUN=',RBF_FUN,';']);
maxgen = 200000+3000*Dim;
minerror = 1e-20;
RBF_FUN_Arr{i} = RBF_FUN;
[best_pos,bestever] = JADE(Dim, maxgen, RBF_FUN, minerror, elite_data);
predict_pos(i,:) = best_pos;
predict(i) = FUN(best_pos);

% Search predicted optimal based on surrogate RBF-G [model 2]
i = 2;
ghxd=real(sqrt(hx.^2*ones(size(hx'))+ones(size(hx))*(hx').^2-2*hx*(hx')));
spr=max(max(ghxd))/(Dim*nd)^(1/Dim)*2^12;
net=newrbe(hx',hf',spr);
RBF_FUN=['RBF_FUN',num2str(i)];
eval([RBF_FUN,'=@(x) sim(net,x'');']);
eval(['RBF_FUN=',RBF_FUN,';']);
maxgen = 200000+3000*Dim;
minerror = 1e-20;
RBF_FUN_Arr{i} = RBF_FUN;
[best_pos,bestever] = JADE(Dim, maxgen, RBF_FUN, minerror, elite_data);
predict_pos(i,:) = best_pos;
predict(i) = FUN(best_pos);

% Search predicted optimal based on ensemble model of RBF-L and RBF-G [model 3] 
i = 3;
maxgen = 200000+3000*Dim; 
minerror = 1e-20;
RBF_FUN_Arr_ensemble = RBF_FUN_Arr([1,2]);
ensemble = @(x)ensemble_RBF_FUN(x,RBF_FUN_Arr_ensemble);
RBF_FUN_Arr{i} = ensemble;
[best_pos, bestever] = JADE(Dim, maxgen, ensemble, minerror, elite_data);
predict_pos(i,:) = best_pos;
predict(i) = FUN(best_pos);

% Search predicted optimal based on ensemble model of many RBFs [model 4]
i = 4;
[ensemble_RBF_time,ensemble_model] = ensemble_RBF_smooth(Dim,[hx,hf],UB,LB);
maxgen = 200000+3000*Dim; 
minerror = 1e-20;
RBF_FUN_Arr{i} = ensemble_model;
[best_pos,bestever] = JADE(Dim, maxgen, ensemble_model, minerror, elite_data);
predict_pos(i,:) = best_pos;
predict(i) = FUN(best_pos);

% calculate distance to center point based on predicted position
for k = 1:10
    if k > num_elite_data
        break;
    end
    R_temp = pdist2(elite_data(k,:), predict_pos);
    R = R_temp + R;             %%%%%%%%%%%%%%%%%% Selection criteria 3 [R]
end

% calculate mutual evaluation based on predicted position
for m = 1:model_num
    RBF_FUN = RBF_FUN_Arr{m};
    co_E_temp(m,:) = RBF_FUN(predict_pos);
    co_E_temp(m,m) = 0;
end
co_Evaluate = sum(co_E_temp,1); %%%%%%%%%%%%%%%%%% Selection criteria 4 [co_E]

% 2.4 normalize to [0,1]
R_N = mapminmax(R, 0, 1);
co_Evaluate_N = mapminmax(co_Evaluate, 0, 1);
test_error_N = mapminmax(test_error, 0, 1); 
sort_error_N = mapminmax(sort_error, 0, 1);

% 2.5 calculate loss based different weight
predict_N = mapminmax(predict, 0, 1);

loss_arr(1,:) = [0 1 1 1]; % model 1
loss_arr(2,:) = [1 0 1 1]; % model 2
loss_arr(3,:) = [1 1 0 1]; % model 3
loss_arr(4,:) = [1 1 1 0]; % model 4
loss_arr(5,:) = 1 * test_error_N + 0.05 * co_Evaluate_N;

% 2.6 select and result
for i = 1: size(loss_arr,1)
    [~, loss(i)] = min(loss_arr(i,:));
    model_selected =  loss(i);
    model_used(i) = model_selected;
    model_rank(i) = sum(predict <= predict(model_selected));
    pred(i) = predict(model_selected);
    result(i) = predict(model_selected);
    pos(i,:) = predict_pos(model_selected,:);
end
totaltime = etime(clock,t1);
time = totaltime;

%% output
bestFitness=pred(end);
bestP=pos(end,:);

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