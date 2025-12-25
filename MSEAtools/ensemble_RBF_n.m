function [time, predict] = ensemble_RBF_n(c,data,bu,bd)
tic;

nd=size(data,1); % Number of neurons of RBF models
dim = size(data,2)-1;
Q=10; % Ensemble size

for i = 1:Q
    S_size = floor(size(data,1)/10);
    del = randperm(size(data,1),S_size);
    train_data = data;
    train_data(del,:) = [];
    train_hf = train_data(:,dim+1);  train_hx = train_data(:,1:dim);
    ghxd=real(sqrt(train_hx.^2*ones(size(train_hx'))+ones(size(train_hx))*(train_hx').^2-2*train_hx*(train_hx')));
    spr=max(max(ghxd))/(dim *nd)^(1/dim );
    net=newrbe(train_hx',train_hf',spr);
    RBF_FUN=@(x) sim(net,x');
    RBF_FUN_Arr{i} = RBF_FUN;
end

predict = @(x)ensemble_RBF_FUN(x,RBF_FUN_Arr);

toc;
time=toc;

end