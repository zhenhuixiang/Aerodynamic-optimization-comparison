%% Experiments for Paper
clc;
clear;
addpath(genpath(pwd));
% Test Functions
% 'ROSENBROCK'; 'ACKLEY'; 'GRIEWANK'; 'Ellipsoid'; 'CEC05_f10';'CEC05_f16'; 'CEC05_f19';
TestFuns = { 'GRIEWANK'};    % The objective functions to be tested
dims = [10];                        % Dimensions to be tested
Runs = 20;                              % Number of runs

d = size(dims,2);
o = length(TestFuns);

bag_gsamp1 = {};                        % pack gsamp1
bag_time_cost = {};                     % pack time cost

% runs according to dims and objs.
for i = 1:d
    for j = 1:o
        fname = cell2mat(TestFuns(j));                  
        FUN=@(x) feval(fname,x); 
        [Xmin, Xmax] = variable_domain(fname); 
        LB = repmat((Xmin),1,dims(i));             
        UB = repmat((Xmax),1,dims(i),1);
        [ gsamp1,time_cost] = RUN_TSDDEO(Runs,dims(i),FUN, LB, UB, fname);
        
        % Each line contains result of a objective function with different dimensions
        bag_result(j,i) = {gsamp1(end)};
        bag_gsamp1(j,i) = {gsamp1};
        bag_time_cost(j,i) = {time_cost};
    end
end
save Result                     