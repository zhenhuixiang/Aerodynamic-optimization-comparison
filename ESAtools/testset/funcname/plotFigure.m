close all;
%% 1
figure;
fname = 'ACKLEY'
[Xmin, Xmax] = variable_domain(fname);
arr = [1:1000]/1000;
x = (arr*(Xmax-Xmin) + Xmin)'
y = ACKLEY(x)'

samples_x = [-30.7, -20.3, -5.5, 10.2, 26.8]'
samples_y = ACKLEY(samples_x)'
data = [samples_x samples_y]
%
plot(x,y)
hold on
scatter(samples_x,samples_y)

% R B F network
D =1; gs=5;
ghf=samples_y;     ghx=samples_x;
ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
spr=max(max(ghxd))/(D*gs)^(1/D);
net=newrbe(ghx',ghf',spr);        
modelFUN=@(x) sim(net,x');
y_hat = modelFUN(x)

hold on 
plot(x,y_hat)
[m,p] = min(y_hat)
new_sample_x = x(p)
new_sample_y = ACKLEY(new_sample_x)

%% 2
figure;
plot(x,y)
hold on
samples_y = [samples_y; new_sample_y] 
samples_x = [samples_x; new_sample_x] 
Xmin = min(samples_x);
Xmax = max(samples_x);
x = (arr*(Xmax-Xmin) + Xmin)'
[~,id]=sort(samples_y);                          
D =1; gs=5;
gx=samples_x(id(1:gs),:);  gy=samples_y(id(1:gs));           % model_1 training samples
scatter(gx,gy)
% R B F network
ghf=gy;     ghx=gx;
ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
spr=max(max(ghxd))/(D*gs)^(1/D);
net=newrbe(ghx',ghf',spr);        
modelFUN=@(x) sim(net,x');
y_hat = modelFUN(x)

hold on 
plot(x,y_hat)
[m,p] = min(y_hat)
new_sample_x = x(p)
new_sample_y = ACKLEY(new_sample_x)

%% 3
figure;
plot(x,y)
hold on
samples_y = [samples_y; new_sample_y] 
samples_x = [samples_x; new_sample_x] 
Xmin = min(samples_x);
Xmax = max(samples_x);
x = (arr*(Xmax-Xmin) + Xmin)'
[~,id]=sort(samples_y);
D =1; gs=5;
gx=samples_x(id(1:gs),:);  gy=samples_y(id(1:gs),:);           % model_1 training samples
scatter(gx,gy)
% R B F network
ghf=gy;     ghx=gx;
ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
spr=max(max(ghxd))/(D*gs)^(1/D);
net=newrbe(ghx',ghf',spr);        
modelFUN=@(x) sim(net,x');
y_hat = modelFUN(x)

hold on 
plot(x,y_hat)
[m,p] = min(y_hat)
new_sample_x = x(p)
new_sample_y = ACKLEY(new_sample_x)

%% 4
figure;
plot(x,y)
hold on
samples_y = [samples_y; new_sample_y] 
samples_x = [samples_x; new_sample_x] 
Xmin = min(samples_x);
Xmax = max(samples_x);
x = (arr*(Xmax-Xmin) + Xmin)'
[~,id]=sort(samples_y);                        
D =1; gs=5;
gx=samples_x(id(1:gs),:);  gy=samples_y(id(1:gs),:);           % model_1 training samples
scatter(gx,gy)
% R B F network
ghf=gy;     ghx=gx;
ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
spr=max(max(ghxd))/(D*gs)^(1/D);
net=newrbe(ghx',ghf',spr);        
modelFUN=@(x) sim(net,x');
y_hat = modelFUN(x)

hold on 
plot(x,y_hat)
[m,p] = min(y_hat)
new_sample_x = x(p)
new_sample_y = ACKLEY(new_sample_x)

%% 5
figure;
plot(x,y)
hold on
samples_y = [samples_y; new_sample_y] 
samples_x = [samples_x; new_sample_x] 
Xmin = min(samples_x);
Xmax = max(samples_x);
x = (arr*(Xmax-Xmin) + Xmin)'
[~,id]=sort(samples_y);                        
D =1; gs=5;
gx=samples_x(id(1:gs),:);  gy=samples_y(id(1:gs),:);
scatter(gx,gy)
% R B F network
ghf=gy;     ghx=gx;
ghxd=real(sqrt(ghx.^2*ones(size(ghx'))+ones(size(ghx))*(ghx').^2-2*ghx*(ghx')));
spr=max(max(ghxd))/(D*gs)^(1/D);
net=newrbe(ghx',ghf',spr);        
modelFUN=@(x) sim(net,x');
y_hat = modelFUN(x)

hold on 
plot(x,y_hat)
[m,p] = min(y_hat)
new_sample_x = x(p)
new_sample_y = ACKLEY(new_sample_x)

