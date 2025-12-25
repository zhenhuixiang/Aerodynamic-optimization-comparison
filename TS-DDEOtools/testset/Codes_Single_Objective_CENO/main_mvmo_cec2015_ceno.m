%==========================================================================
%    MVMO-based solution for CEC 2015 Special Session and Competition  
%    on Bound Constrained Single-Objective Computationally Expensive 
%    Numerical Optimization
%==========================================================================

%                                  Reference
%--------------------------------------------------------------------------
% [1] Q. Chen, B. Liu, Q. Zhang, J. J. Liang, P. N. Suganthan, and B. Y. Qu, 
%     "PProblem Definitions and Evaluation Criteria for CEC 2015 Special
%      Session and Competition on Bound Constrained Single-Objective
%      Computationally Expensive Numerical Optimization," Nov. 2014
% [Online] Available at http://www.ntu.edu.sg/home/epnsugan/
%--------------------------------------------------------------------------

%By: Dr.-Ing. José L. Rueda & Prof. István Erlich
%30.10.2014 


matlabpool close force local

close all
clear all
clc

global proc

optrepnum=1;
best_score=Inf;
best_run_no=0;

for icont=1:optrepnum
    


proc.fbest_stats=zeros(15,5);
proc.mvmo_time=zeros(15,1);

%For algorithm Complexity
tStart1 = cputime;
for i=1:1000000
    x= 0.55d0+i;
    x=x+x; 
    x=x/2;
    x=x*x; 
    x=sqrt(x); 
    x=log(x);
    x=exp(x); 
    y=x/(x+2);
end
T0 = cputime-tStart1;


ranges = [-100,100];


%                Preliminary definitions & Parallelization
%--------------------------------------------------------------------------
refresh=100; %Printing step
algorithm_name='mvmo';
algorithm_hd=str2func(algorithm_name);
test_bed_OPF_hd=str2func('test_func_calc');
args{1}=refresh;
args{2}=algorithm_name;
args{3}=10; %30; %Problem dimension
args{4}=50*args{3}; %Maximum function evaluations
% args{5}=repmat(-100,1,args{3}); %Min limits of control variables 
% args{6}=repmat(100,1,args{3});%Max limits of control variables  
args{7}=[100:100:1500]'; %Theoretical optimum value
args{8}=20;%20; %Optimization runs
args{9}=0;%Printing 1=yes, 0=no

run_in_parallel=0; %1 = yes,  0 = no

v=ver;
% Whether or not MATLAB parallel computing toolbox is installed.
toolbox_installed=any(strcmp('Parallel Computing Toolbox',{v.Name}));
% If yes, the cluster 
% is deactivated at first.
if toolbox_installed
    isOpen=matlabpool('size')>0;
    if isOpen
        matlabpool close
    end
% If not, trails will be 
% processed in sequence.
else
    run_in_parallel=0;
end

% Activation of cluster 
% consisting of NumWorkers 
% computational cores/threads.
if run_in_parallel
    NumWorkers=7;
    local_sched=findResource('scheduler','type','local');
    local_sched.ClusterSize=NumWorkers;
    isOpen=matlabpool('size')>0;
    if ~isOpen
        matlabpool(NumWorkers);
    end
end
%--------------------------------------------------------------------------


%                           Running optimization
%--------------------------------------------------------------------------
fnum_count=zeros(15,1);


for func_num=1:15 %From 1 to 15 - Problem number
    
    lu_lim = ranges;
    args{5}=lu_lim(1)*ones(1,args{3}); %Min limits of control variables
    args{6}=lu_lim(2)*ones(1,args{3}); %Max limits of control variables 
        
    test_bed_OPF_hd(func_num,1,icont,args); 
    parfor krun=1:args{8}
        test_bed_OPF_hd(func_num,0,icont,args);
        % Call to your implementation.
        feval(algorithm_hd,test_bed_OPF_hd,func_num,krun,icont,args);
%         fprintf('Run %d finished.\n',krun);
    end
    test_bed_OPF_hd(func_num,987,icont,args);   
       
  
    fnum_count(func_num,:)=func_num;  
    
end


%--------------------------------------------------------------------------


% Parallelization
% Deactivation of 
% shared-memory session.
if run_in_parallel
    isOpen=matlabpool('size')>0;
	if isOpen
        matlabpool close
	end
end


tot_score=sum(proc.fbest_stats(:,3))+sum(proc.fbest_stats(:,4));

if icont==1
    best_score=tot_score;
    best_run_no=icont;
    fbest_stats=[fnum_count, proc.fbest_stats];
    mvmo_time=[fnum_count, proc.mvmo_time/T0];
    save(proc.problem_filename_stats, 'fbest_stats','-ASCII')
    save(proc.problem_filename_time, 'mvmo_time','-ASCII')
    fprintf('Total score  far for %d -D Problems is %17.5E\n', args{3}, tot_score); 
else
    if tot_score<best_score
        rmdir([num2str(best_run_no) '_output_data_' args{2}],'s') 
        delete([num2str(best_run_no) '_output_data_' args{2} '.zip'])
        delete([num2str(best_run_no) '_' args{2} '_statistics_' num2str(args{3})  '.txt'])
        delete([num2str(best_run_no) '_' args{2} '_time_' num2str(args{3})  '.txt'])
        best_score=tot_score;
        best_run_no=icont;  
        fbest_stats=[fnum_count, proc.fbest_stats];
        mvmo_time=[fnum_count, proc.mvmo_time/T0];
        save(proc.problem_filename_stats, 'fbest_stats','-ASCII')
        save(proc.problem_filename_time, 'mvmo_time','-ASCII')
        fprintf('Total score  far for %d -D Problems is %17.5E\n', args{3}, tot_score);  
    else
        rmdir([num2str(icont) '_output_data_' args{2}],'s') 
        delete([num2str(icont) '_output_data_' args{2} '.zip'])
    end
end

if icont==optrepnum
        fprintf('%d -D Problems: The best repetition is %d , which entails the best overall score  is %17.5E\n', args{3},best_run_no, best_score); 
end




end
