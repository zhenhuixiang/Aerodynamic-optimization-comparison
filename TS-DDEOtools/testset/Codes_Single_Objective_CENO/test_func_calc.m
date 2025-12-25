function [f,oo_xx,g_sum,x_out]=test_func_calc(iii,jjj,kkk,args,varargin)
    global proc 
    global ps
    global res
   
    f=NaN;
    oo_xx=0;
    g_sum=0;
    x_out=NaN;

    
    proc.test_case=iii;
    proc.i_run=jjj;
    proc.refresh=args{1};
    proc.algorithm_name=args{2}; 
    proc.printing=args{9}; 


     
    if proc.i_run==987
        ps.D=args{3}; %Problem dimension
        proc.n_eval=args{4}; %Maximum function evaluations
        proc.n_run=args{8};
        proc.problem_filename_stats=[num2str(kkk) '_' proc.algorithm_name '_statistics'...
                                     '_' num2str(ps.D)  '.txt']; 
        proc.problem_filename_time=[num2str(kkk) '_' proc.algorithm_name '_time'...
                                     '_' num2str(ps.D)  '.txt'];
        
        common_str=[proc.algorithm_name '_' num2str(proc.test_case) '_'...
                    num2str(ps.D)];
        fitness=NaN*zeros(proc.n_eval,proc.n_run);
        variables=NaN*zeros(proc.n_run,ps.D);
        
        for i_run=1:proc.n_run
            common_str1=[common_str '_run_' num2str(i_run)];
            fitness(:,i_run)=dlmread([common_str1 '_fitness.txt']);
            delete([common_str1 '_fitness.txt'])
            complexity(:,i_run)=dlmread([common_str1 '_complexity.txt']);
            delete([common_str1 '_complexity.txt'])
            variables(i_run,:)=dlmread([common_str1 '_variables.txt']);
            delete([common_str1 '_variables.txt'])

        end
        
        fitness_save=fitness(round([0.01:0.01:0.09 0.1:0.1:1].*proc.n_eval),:);
        fitness_save=fitness_save';
        
        proc.fbest_stats(proc.test_case,:)=[min(fitness(end,:)),max(fitness(end,:)),...
                     median(fitness(end,:)),mean(fitness(end,:)),std(fitness(end,:))];
        proc.mvmo_time(proc.test_case,:)=mean(complexity);
        
%         fitness_save=fitness([1 100:100:end],:);
%         objective_save=objective([1 100:100:end],:);
        
        mkdir([num2str(kkk) '_output_data_' proc.algorithm_name]);
        save([cd [filesep num2str(kkk) '_output_data_' proc.algorithm_name filesep] common_str '_fitness_full.txt'],...
            'fitness','-ASCII')
        
        save([cd [filesep num2str(kkk) '_output_data_' proc.algorithm_name filesep] common_str '.txt'],...
            'fitness_save','-ASCII')

        save([cd [filesep num2str(kkk) '_output_data_' proc.algorithm_name filesep] common_str '_variables.txt'],...
            'variables','-ASCII')
        
        zip([cd filesep num2str(kkk) '_output_data_' proc.algorithm_name '.zip'],[cd filesep num2str(kkk) '_output_data_' proc.algorithm_name])
        
    end
 
    
    if (~numel(varargin)>0)&(proc.i_run==1)
        addpath(genpath(cd));
        rand('state',sum(proc.i_run*100*clock));
        if proc.printing
            fprintf('Random number generation seed initialized ');
            fprintf('according to "rand(''state'',(trial=%d)*sum(100*clock))".\n',proc.i_run);
            fprintf('Test function: %d\n',proc.test_case);
            fprintf('Dimension: %d\n',args{3});
        end
        proc.t1=cputime;
    end

    
    if (~numel(varargin)>0)&(proc.i_run==0)
        ps.D=args{3}; %Problem dimension
        proc.n_eval=args{4}; %Maximum function evaluations
        %Min and max limits of control variables  
        ps.x_min=args{5};
        ps.x_max=args{6};
%         proc.ofref=args{7}; %Theoretical optimum value
    end
    

    if (numel(varargin)>0) &(proc.i_eval==0)
        res{proc.test_case}{proc.i_run}.fitness=inf*ones(proc.n_eval,1);
        res{proc.test_case}{proc.i_run}.variables=NaN*zeros(1,ps.D);
    end
    
    
    if numel(varargin)>0    
        proc.ofref=args{7}; %Theoretical optimum value
        [f,x_out] = FE(varargin{1});
    end
    
end

function[f,x_out]=FE(x_in)
    global proc 
    global ps 
    global res

    
        f = cec15problems('eval',x_in',proc.test_case)-proc.ofref(proc.test_case);
        x_out = x_in;
        proc.i_eval=proc.i_eval+1;
        

        if proc.i_eval<=proc.n_eval&&proc.i_eval>=1
            tmp=res{proc.test_case}{proc.i_run}.fitness...
                (proc.last_improvement);
                  
               
            if f<tmp
                proc.last_improvement=proc.i_eval;
             	res{proc.test_case}{proc.i_run}.fitness...
                    (proc.last_improvement)=f;

                res{proc.test_case}{proc.i_run}.variables=x_in;

            else
                res{proc.test_case}{proc.i_run}.fitness(proc.i_eval)=...
                    res{proc.test_case}{proc.i_run}.fitness...
                        (proc.last_improvement);

            end
            
 
            
            if ((proc.i_eval==1)||(mod(proc.i_eval,proc.refresh)==0))& proc.printing
                    fprintf('Func: %5d,   Dim: %5d,   trial: %5d,   i_eval: %7d,   f_best: %12.7E \n',...
                    proc.test_case,ps.D,proc.i_run,proc.i_eval,res{proc.test_case}{proc.i_run}.fitness(proc.i_eval));
            end
        end
        
        if proc.i_eval>=proc.n_eval
            proc.finish=1;
            
            proc.finish=1;
            t2=cputime;
            if ~isfield(proc,'t1')
                proc.t1=0;
            end
            res{proc.test_case}{proc.i_run}.complexity=t2-proc.t1;

            common_str=[proc.algorithm_name '_' num2str(proc.test_case)...
                        '_' num2str(ps.D) '_run_' num2str(proc.i_run)];
            val=[];
            val=res{proc.test_case}{proc.i_run}.fitness;
        	save([common_str '_fitness.txt'],'val','-ASCII')

            val=[];
        	val=res{proc.test_case}{proc.i_run}.variables;
            save([common_str '_variables.txt'],'val','-ASCII')
            
            val=[];
            val=res{proc.test_case}{proc.i_run}.complexity;
            save([common_str '_complexity.txt'],'val','-ASCII')

            return;
        end

end