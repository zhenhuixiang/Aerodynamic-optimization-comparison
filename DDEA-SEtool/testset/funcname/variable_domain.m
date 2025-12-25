%%%*********************************************************************************************%%%
%% Implementation of Surrogate-assisted Hierarchical Particle Swarm Optimization (SHPSO)
%% H. Yu, Y. Tan, J. Zeng, C. Sun, Y. Jin, Surrogate-assisted hierarchical 
%% particle swarm optimization, Information Sciences, 454-455 (2018) 59-72.
%%%*********************************************************************************************%%%
%% This paper and this code should be referenced whenever they are used to 
%% generate results for the user's own research. 
%%%*********************************************************************************************%%%
%% This matlab code was written by Haibo Yu
%% Please refer with all questions, comments, bug reports, etc. to tyustyuhaibo@126.com
% 
%% Test functions' variable domain
function [Xmin, Xmax] = variable_domain(fname)
    switch fname
        case 'GRIEWANK'
            Xmin=-600; Xmax=600;
        case 'ACKLEY'
            Xmin=-32.768; Xmax=32.768;
        case 'ROSENBROCK'
            Xmin=-2.048; Xmax=2.048;
        case 'RASTRIGIN'
            Xmin=-5.12; Xmax=5.12;
        case 'ELLIPSOID'
            Xmin=-5.12; Xmax=5.12;
        case 'S_GRIEWANK'
            Xmin=-600; Xmax=600;
        case 'S_ACKLEY'
            Xmin=-32.768; Xmax=32.768;
        case 'S_ROSENBROCK'
            Xmin=-2.048; Xmax=2.048;
        case 'S_RASTRIGIN'
            Xmin=-5.12; Xmax=5.12;
        case 'S_ELLIPSOID'
            Xmin=-5.12; Xmax=5.12;
        case 'CEC05_F10'
            Xmin=-5; Xmax=5;
        case 'CEC05_F16'
            Xmin=-5; Xmax=5;
        case 'CEC05_F19'
            Xmin=-5; Xmax=5;
        otherwise
            Xmin=-5; Xmax=5;
    end
end
% 变量上下限
%% For 'FITNESS' 
% Xmin=-32.768;Xmax=32.768;%---Ackley 
% Xmin=-600;Xmax=600;%---Griewank
% Xmin=-2.048;Xmax=2.048;%---Rosenbrock
% Xmin=-5.12;Xmax=5.12;%---Ellipsoid 

%% For 'benchmark_func'
% Xmin=-5;Xmax=5; % cec05func 10
% Xmin=-5;Xmax=5; % cec05func 16
% Xmin=-5;Xmax=5; % cec05func 19
