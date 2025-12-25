% 测试DDEA在测试函数上的效果
clc;
clear;
addpath(genpath(pwd));
warning off all

funcArr={'Rastrigin'; 'ROSENBROCK'; 'ACKLEY'; 'GRIEWANK'; 'Ellipsoid'; }; 
dimsArr = [10]; 
runs = 5;
o = length(funcArr);
d = size(dimsArr,2);

for i = 1:d
	dims = dimsArr(i);
    for j = 1:o   
        fname = cell2mat(funcArr(j));  
        result = [];
        time = [];
        for r = 1:runs
            % funcArr = {'Rastrigin'; 'ROSENBROCK'; 'ACKLEY'; 'GRIEWANK'; 'Ellipsoid'; 'CEC05_f10'; 'CEC05_f16';  'CEC05_f19'; };    % The objective functions to be tested
            FUN=@(x) feval(fname,x);
            [Xmin, Xmax] = variable_domain(fname);
            LB = repmat((Xmin),1,dims);
            UB = repmat((Xmax),1,dims,1);
            initial_sample_size = 11*dims;
            sam=repmat(LB,initial_sample_size,1)+(repmat(UB,initial_sample_size,1)-repmat(LB,initial_sample_size,1)).*lhsdesign(initial_sample_size,dims);
            fitness = FUN(sam);
            hx=sam; hf=fitness';
            DATA = [hx, hf];

            t1 = clock;

            % the boundary of decision variable
            Dim = length(LB);
            [time,P,gbest,predict] = DDEA_SE(Dim,DATA,UB,LB);

            totaltime = etime(clock,t1);
            result(r) = FUN(P);
            time(r) = totaltime;
        end
        best_result = min(result)
        worst_result = max(result)
        mean_result = mean(result)
        median_result = median(result);
        std_result    = std(result);
        mean_time = mean(time);
        out1 = [best_result,worst_result,mean_result,median_result,std_result];
        save(strcat('result/funcname/',fname,' runs=',num2str(runs),' Dim=',num2str(dims)));
    end
end