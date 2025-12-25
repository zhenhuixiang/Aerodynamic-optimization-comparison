function mvmo(fhd,iii,jjj,kkk,args)
    global proc
    global ps %disc
    global parameter table
    global changed_best
    persistent amax
         
    
    proc.finish=0;
    proc.i_eval=0;
	proc.last_improvement=1;

    % only the following parameters are required from the user    
    parameter.n_par=80       ;
    parameter.n_tosave=15 ;
    parameter.fs_factor_start=1;
    parameter.fs_factor_end=35             ;
    parameter.delta_Shape_dyn=0.1  ;
    parameter.local_prob= 00  ;
    min_eval_LS =round(0.95 *proc.n_eval);
    max_eval_LS =round(1.00*proc.n_eval+10);
    parameter.ratio_gute_max=0.9        ; 
    parameter.ratio_gute_min=0.1 ;
    parameter.n_random_ini =9  ; %round(ps.D/6);
    parameter.n_random_last=3 ;
    local_max=1   ;  % max number of local search runs allowed
    % end of parameter input

    parameter.scaling=(ps.x_max-ps.x_min);      
    indpendent_runs=parameter.n_par*2 ;
 

    %%----------------- Create initial random population ------------------   
    xx=zeros(parameter.n_par,ps.D);
    x_norm=xx;
    for iijj=1:parameter.n_par %Start initialization: x_normalized
          for jjkk=1:ps.D
              xx(iijj,jjkk)=ps.x_min (jjkk) + rand*(ps.x_max(jjkk)-ps.x_min(jjkk));
          end
%           xx(iijj,1:ps.D) = -100+200*rand(ps.D,1);
        x_norm(iijj,:)=(xx(iijj,:)-ps.x_min)./parameter.scaling;
    end % End initialization
    
    x_normalized=x_norm;
    
    %% ------------------ Initialize control parameters --------------------
    n_eval = proc.n_eval;
    n_par = parameter.n_par; 
    D=ps.D ; 
    n_to_save = parameter.n_tosave; 
    fs_factor_start = parameter.fs_factor_start;
    fs_factor_end = parameter.fs_factor_end;
    Shape_dyn = ones(n_par,D);%
    IX=ones(n_par); 
    shape = zeros(n_par,D);
    for i=1:n_par
        IX(i)=i;
        for k=1:ps.D
            Shape_dyn(i,k)=350  ;  
            shape(i,k)= 350;
        end
    end
    delta_Shape_dyn = parameter.delta_Shape_dyn; 
    local_search0_percentage=parameter.local_prob;
    ratio_gute_min=parameter.ratio_gute_min;
    ratio_gute_max=parameter.ratio_gute_max;
    n_randomly_ini = parameter.n_random_ini;
    n_randomly_last = parameter.n_random_last;
    local_i=0;    % local seurch runs counter
    
    %% --------------------- Data structure for the table ---------------------
    table.bests = NaN*zeros(parameter.n_tosave,ps.D,n_par);
    table.fitness = Inf*ones(parameter.n_tosave,1,n_par);
    table.objective = Inf*ones(parameter.n_tosave,1,n_par);
    table.feasibility = zeros(parameter.n_tosave,1,n_par);
    
    %% ----------------------------- Mapping ----------------------------------
    local_search=zeros(n_par);
    meann = x_normalized;
        meann(:)=0.5;
    meann_app = meann;
   %% ------------------------ Variable selection ----------------------------
    izm = zeros(1,n_par);
    izz = zeros(1,n_par);
    considered = true(n_par,ps.D);
    probab=ones(n_par,ps.D);
    values_noneq = zeros(n_to_save,1);
    
   %% ---------------------- Checking correctness ------------------------   
   
    if (n_randomly_last<1)
        n_randomly_last=1;
    end
    
    if (n_randomly_last>ps.D)
        n_randomly_last=ps.D;
    end
    
    if (n_randomly_ini>ps.D)
        n_randomly_ini=ps.D;
    end  
    if (n_eval<=0)
        n_eval=100000d0;
    end

    if (fs_factor_start<=0)
        fs_factor_start=1d0;
    end
    
    if (fs_factor_end<=0)
        fs_factor_end=1d0;
    end

    fs_factor0=fs_factor_start;
    
    yes_n_randomly=true;
    if (n_randomly_ini==n_randomly_last)
        n_randomly=n_randomly_last;
        yes_n_randomly=false;
    end
    
    if (n_to_save<=1)
        n_to_save=2d0;
    end
    
    if (delta_Shape_dyn<=0)
        delta_Shape_dyn=0.2d0;
    end
    
    delta_Shape_dyn0=delta_Shape_dyn;
    delta_Shape_dyn1=1.d0+delta_Shape_dyn0;
    
    yes_fs_factor=true;
    if (fs_factor_start==fs_factor_end)
        yes_fs_factor=false;
    end

    %% ----------------------------- Counters ---------------------------------

    no_in = zeros(1,n_par);
    no_inin = zeros(1,n_par);
    amax=0.d0;
    
    local_search0=local_search0_percentage/100; % Probability of local search (percentage / number of optimization variables)
    goodbad=zeros(n_par);
    firsttime=true;
    

    delta_nrandomly=n_randomly_ini-n_randomly_last;
    A=zeros(n_par,1); 
    local_search_called=0;
   
  while 1       
        %Evaluating the particles.....
            ff=proc.i_eval/n_eval ;
            ff2=ff*ff;
            vvqq=10.d0^(-(3.3d0+ff*15.d0));
                        
            %Determining the relative number of particles belonging to the group of
            %good particles
            border_gute0=ratio_gute_max-ff*(ratio_gute_max-ratio_gute_min);
            border_gute=round(n_par*border_gute0) ;
            if border_gute < 3 || n_par-border_gute < 1
                border_gute=n_par ;
             end                     
            %Selecting the subset of variables to be mutated
            if yes_n_randomly
                n_randomly_X=round(n_randomly_ini-ff*delta_nrandomly);
                n_randomly=round(n_randomly_last+rand*(n_randomly_X-n_randomly_last));
            end
            %Quadratic variation of fs_factor0
            if yes_fs_factor
                fs_factor0=fs_factor_start+ff2*(fs_factor_end-fs_factor_start); 
            end
            ipx=0;

       while ipx  < n_par 
             ipx=ipx+1 ;
              ipp=IX(ipx);
            considered(ipp,1:ps.D) = false;
           %Denormalize to the original [min, max] range 
              x_normalized(ipp,:) = ps.x_min+parameter.scaling.* x_normalized(ipp,:);
             if  local_search(ipp) && local_i <  local_max && local_search0 > 0.d0 
                [msgstr, msgid] = lastwarn ;
                TFrcond = strcmp('MATLAB:nearlySingularMatrix',msgid); % Only informative from 'fmincon' function 
                if TFrcond~=0
                    rcond_value0=str2num(msgstr(81:end-1));
                end

                [ffx,oox,~,x_normalized(ipp,:),FEVALS] = LocalSearchMVMOSH(fhd,iii,jjj,kkk,args,x_normalized(ipp,:),proc.n_eval-proc.i_eval); %Local search
                local_search(ipp)=0 ;
                local_i=local_i+1  ; 
                local_search_called=local_search_called+1;
                xt.fitness=ffx ; % Constraint handling outsourced so far
                xt.objective=oox;
                
                if xt.fitness==xt.objective
                    xt.feasibility= true; 
                else
                    xt.feasibility= false;
                end

                [~, msgid1] = lastwarn;
                TFrcond1 = strcmp('MATLAB:nearlySingularMatrix',msgid1);
             else
                [ffx,oox,~,x_normalized(ipp,:)]=feval(fhd,iii,jjj,kkk,args,x_normalized(ipp,:)); %Problem evaluation
                xt.fitness=ffx; % Constraint handling outsourced so far   
                xt.objective=oox;
                if xt.fitness==xt.objective
                    xt.feasibility= true; %
                else
                    xt.feasibility= false;
                end
            end
            
            if proc.finish
                return;
            end
            
             x_normalized(ipp,:) = (x_normalized(ipp,:)-ps.x_min)./parameter.scaling;
                
            % Store the n-best solution corresponding to the corresponding particle's archive
            Fill_solution_archive();
            meann_app(ipp,:)= meann(ipp,:);
            if  proc.i_eval  > indpendent_runs  ; %&&  border_gute < n_par;    
            %Determining the proportion of good particles
                    if changed_best || firsttime
                     A(1:n_par)=table.fitness(1,:,1:n_par);
                    if  firsttime
                      amax=0.  ;%max(A);
                      firsttime=false  ;
                    end
                    for ia=1:n_par
                      if ~table.feasibility(1,:,ia)
                        A(ia)=A(ia)+amax;
                      end
                    end 
                    [~,IX] = sort(A);
                   end
                      verybest=IX(1);
                      goodbad(IX(border_gute+1:n_par))=0;
                      goodbad(IX(1:border_gute))=1;
                    %Multi-parent strategy for bad particles  
                 if ~goodbad(ipp)    
                        [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)      ;  
                        bbb=1.1d0+(rand-0.5d0)*2.0d0 ; 
                        beta1 = 3.d0*bbb*((1.d0+2.5*ff2)*rand - (1.d0-ff2)*0.8d0 )  ;  
                        beta10=0.d0     ;
                        while beta1 ~= beta10
                          beta10=beta1 ;
                          for jx=1:ps.D
                                    ccc=table.bests(1,jx,bestp)-table.bests(1,jx,worstp);
                                   if abs(ccc) >  1.d-15  
                                    x_normalized(ipp,jx) =table.bests(1,jx,onep1)+ beta1*ccc;   
                                    else
                                     x_normalized(ipp,jx)=table.bests(1,jx,onep1)    ;
                                    end
                                 if table.bests(1,jx,bestp) < 0.85 && table.bests(1,jx,bestp) > 0.15 && rand > 0.15 
                                   while x_normalized(ipp,jx) > 1.0d0 ||  x_normalized(ipp,jx)<0.0d0    
                                   [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)      ;       
                                    ccc=table.bests(1,jx,bestp)-table.bests(1,jx,worstp);
                                    if abs(ccc) >  1.d-15    
                                     bbb=1.1d0+(rand-0.5d0)*2.2d0 ;%*1.5d0;
                                     beta1 = 3.d0*bbb*((1.d0+2.5*ff2)*rand - (1.d0-ff2)*0.8d0 )  ;  
                                     x_normalized(ipp,jx) =table.bests(1,jx,onep1)+ beta1*ccc;  
                                    else
                                     x_normalized(ipp,jx)=table.bests(1,jx,onep1)  ; 
                                    end
                                   end
                                 elseif table.bests(1,jx,bestp) > 0.85 || table.bests(1,jx,bestp) < 0.15 && rand < 0.15 
                                   while x_normalized(ipp,jx) > 1.0d0 ||  x_normalized(ipp,jx) < 0.0d0    
                                    [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)      ;   
                                    ccc=table.bests(1,jx,bestp)-table.bests(1,jx,worstp);
                                    if abs(ccc) >  1.d-15     
                                     bbb=1.1d0+(rand-0.5d0)*2.2d0 ;%*1.5d0;
                                     beta1 = 3.d0*bbb*((1.d0+2.5*ff2)*rand - (1.d0-ff2)*0.8d0 )  ;  
                                     x_normalized(ipp,jx) =table.bests(1,jx,onep1)+ beta1*ccc;    
                                    else
                                     x_normalized(ipp,jx)=table.bests(1,jx,onep1)   ;
                                    end
                                   end
                                 else
                                   if x_normalized(ipp,jx) > 1.0d0    
                                     x_normalized(ipp,jx)=1.0d0;
                                   elseif x_normalized(ipp,jx) < 0.0d0  
                                      x_normalized(ipp,jx)=0.0d0 ;
                                   end
                                 end
                          end 
                        end
                                  meann_app(ipp,1:ps.D) =x_normalized(ipp,1:ps.D); 
                     else
                         x_normalized(ipp,1:ps.D)= 0.999*table.bests(1,1:ps.D,ipp)+0.001d0*table.bests(1,1:ps.D,verybest);   % x_normalized_best(ipp,:); %Local best-based parent assignment for good particles
                   end
            else
                     x_normalized(ipp,1:ps.D)=table.bests(1,1:ps.D,ipp); %Local best-based parent assignment during independent evaluations
            end
            
 %           considered(ipp,1:ps.D) = false;
            if rand < local_search0 &&  proc.i_eval > min_eval_LS &&  proc.i_eval < max_eval_LS && ipp==bestp ; %goodbad(ipp) 
              local_search(ipp)=1;
            else
              local_search(ipp)=0;
            end    
             if local_search(ipp)< 1
             if no_inin(ipp) <=  n_to_save   
               n_randomly=ps.D ;
              end  
             VariableSelect1(); % Call random variable selection strategy
             end
            %Generation of random input for the mapping function
            for ivar=1:D
                if considered(ipp,ivar) 
                    x_normalized(ipp,ivar) = rand;
                    if shape(ipp,ivar) > 0.d0 
                      sss1 = shape(ipp,ivar);
          %            sss2=sss1;
                      delta_ddd_x=delta_Shape_dyn0*(rand-0.5d0)*2.0d0+delta_Shape_dyn1; 
                      if (sss1 > Shape_dyn(ipp,ivar))
                        Shape_dyn(ipp,ivar) = Shape_dyn(ipp,ivar)*delta_ddd_x;
                      else
                        Shape_dyn(ipp,ivar) = Shape_dyn(ipp,ivar)/delta_ddd_x;
                      end
                        if  Shape_dyn(ipp,ivar) > sss1
                            grosser=Shape_dyn(ipp,ivar);
                            kleiner=sss1;
                        else
                            kleiner=Shape_dyn(ipp,ivar);
                            grosser=sss1;
                        end   
                        if   meann_app(ipp,ivar) > 0.5d0    
                                 sss1=   grosser;
                                 sss2=   kleiner;
                        else
                                 sss1=   kleiner ;
                                 sss2=   grosser ;
                        end  
                      fs_factor=fs_factor0*(1.0d0 + ff2*rand) ; 
                      sss1=sss1*fs_factor;
                      sss2=sss2*fs_factor;        
                      x_normalized(ipp,ivar)=...
                            h_function(meann_app(ipp,ivar),sss1,sss2,x_normalized(ipp,ivar)); %Mapping function
                                           
                    end
                end
            end
        end %End n_par loop
       
    end %End while loop        

    
     %% ----------------------- Complementary functions ------------------------
     function [ffx,oox,ggx,xn_out,FEVALS] = LocalSearchMVMOSH(testcase,iii,jjj,kkk,args,xx_yy,FEsAllowed)
        global PPL GGL
        if FEsAllowed <= 0, return, end
        options=optimset('Display','off','algorithm','interior-point','UseParallel','never','MaxFunEvals',FEsAllowed,'FinDiffType','central' ) ;
        [Xsqp, FUN , ~ , output]=...
            fmincon(@(xx_yy)LSearch(xx_yy,testcase,iii,jjj,kkk,args),xx_yy,[],[],[],[],ps.x_min,ps.x_max,[],options);
        
        
        FEVALS=output.funcCount  ;
        for nvar=1:size(xx_yy,2)
            if isnan(Xsqp(1,nvar))
                Xsqp=xx_yy;
                break;
            end
        end
        
        xn_out = Xsqp;
        ffx=FUN;
        oox=PPL; 
        ggx=GGL;
    end

    function J=LSearch(xx_yy2,testcase,iii,jjj,kkk,args)
        global PPL GGL 
        [J,PPL,GGL,~] = feval(testcase,iii,jjj,kkk,args,xx_yy2);
        
    end
        
     
        function Fill_solution_archive()
        
  
        no_in(ipp) = no_in(ipp)+1;
        changed = false;
        changed_best=false;

        if no_in(ipp) ==1 % the solution coming to the table for the first time large
            table.fitness(1:n_to_save,:,ipp) = 1.e200;
            table.feasibility(1:n_to_save,:,ipp) = 0;
             table.bests(1,:,ipp)   = x_normalized(ipp,:) ;
            table.fitness(1,:,ipp) = xt.fitness ;
            table.objective(1,:,ipp) = xt.objective;
            table.feasibility(1,:,ipp) = xt.feasibility    ;  %repmat(,n_to_save,1);

            no_inin(ipp)=no_inin(ipp)+1;
            changed_best=true;
            
        else % not for the first time and check for the update
           i_position =0;
               for ij=1:n_to_save 
                   if (xt.fitness < table.fitness(ij,:,ipp) && xt.feasibility == table.feasibility(ij)) || (table.feasibility(ij,:,ipp) <  xt.feasibility)                                 
                       i_position = ij;
                       changed =true;
                       if (ij<n_to_save)
                           no_inin(ipp) = no_inin(ipp)+1; % how many times good solutions were found   
                       end
                       break;
                   end
               end
        end

        if changed   % if the new individual is better than any archived individual.
                     % Move the individuals and corresponding fitness values
                     % downward so that the individuals are sorted based on the
                     % fitness value in a descending order             
            nnnnnn = n_to_save;
            if i_position==1
              changed_best=true;
            end     
            if (no_inin(ipp) < n_to_save); nnnnnn = no_inin(ipp); end
            isdx=nnnnnn:-1:i_position+1;
            table.bests(isdx,:,ipp) = table.bests(isdx-1,:,ipp);
            table.fitness(isdx,:,ipp)= table.fitness(isdx-1,:,ipp);
            table.objective(isdx,:,ipp)= table.objective(isdx-1,:,ipp);
            table.feasibility(isdx,:,ipp)= table.feasibility(isdx-1,:,ipp);
 
            % save the new best
            table.bests(i_position,:,ipp) = x_normalized(ipp,:);
            table.fitness(i_position,:,ipp) = xt.fitness;
            table.objective(i_position,:,ipp) = xt.objective;
            table.feasibility(i_position,:,ipp) = xt.feasibility;

            % calculation of mean and variance
            if ((no_inin(ipp)>=n_to_save))
                for ivvar = 1:D
                    [meann(ipp,ivvar),shape(ipp,ivvar)] = mv_noneq(nnnnnn,table.bests(1:nnnnnn,ivvar,ipp),meann(ipp,ivvar),shape(ipp,ivvar),vvqq);
                end

            end
        end
        end
    
        function VariableSelect1()
          mode=4 ;
          n_var=ps.D;
          switch mode
            case 1
                for ii=1:n_randomly
                    isrepeat = false;
                    while ~isrepeat
                        inn=round(rand*(n_var-1))+1;
                        if (~considered(ipp,inn))
                            isrepeat = true;
                        end
                    end
                    considered(ipp,inn)=true;
                end
            case 2
                in_randomly=0;
                isrepeat = false;
                izz(ipp)=round(rand*(n_var-1))+1; %NEWWWWW
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                        izz(ipp)=n_var;
                    end
                    considered(ipp,izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly)) 
                        isrepeat = true;
                    end
                end
            case 3
                in_randomly=0;
                izm(ipp)=izm(ipp)-1;
                
                             izm(ipp)=round(rand*(n_var-1))+1;
               
                
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                izz(ipp)=izm(ipp);
                isrepeat = false;
                while ~isrepeat
                    in_randomly=in_randomly+1;
                    if (izz(ipp)< 1) 
                         izz(ipp)=n_var;
                    end
                    considered(ipp,izz(ipp))=true;
                    izz(ipp)=izz(ipp)-1;
                    if (~(in_randomly<n_randomly)) 
                        isrepeat = true;
                    end
                end   
            case 4
                izm(ipp)=izm(ipp)-1;
                if (izm(ipp)< 1) 
                    izm(ipp)=n_var;
                end
                 considered(ipp,izm(ipp))=true;
                 if (n_randomly > 1)  
                    for ii=1:n_randomly-1
                        isrepeat = false;
                        while ~isrepeat
                            inn=round(rand*(n_var-1))+1;
                            if (~considered(ipp,inn))
                                isrepeat = true;
                            end
                        end
                        considered(ipp,inn)=true;
                    end
                 end
            case 5
                  summep=sum(probab(ipp,:));
                  wahr=probab(ipp,:)/summep;
                  SS0=0.d0;
                  SS=zeros(1,n_var);
                  for imm=1:(n_var-1)
                      SS0=SS0+wahr(imm);
                      SS(imm+1)=SS0;
                  end
                  for ijlr=2:n_var
                      wahr(ijlr)=wahr(ijlr)+SS(ijlr);
                  end 
                  for ltt=1:n_randomly
                       isrepeat = false;
                       while  ~isrepeat
                        rnn=rand;
                        for irw=1:n_var
                          if considered(ipp,irw)
                           continue
                          end
                          if (irw==1)
                              unten=0.d0;
                          else
                              unten=wahr(irw-1);
                          end
                         if (rnn>=unten)&&(rnn<wahr(irw))
                              isrepeat = true;
                             considered(ipp,irw)=true;
                             break;
                          end
                        end
                       end
                  end
          end
            
        end
%
    function [vmean,vshape] = mv_noneq(nnnnnn,values,vmean,vshape,vvqq)
            iz =1; 
            values_noneq(iz)=values(1);
            for ii_jj=2:nnnnnn
                izz = iz;
                gleich = false;
                for kk_ii=1:izz
                    if  abs(values_noneq(kk_ii) - values(ii_jj)) <  vvqq ; 
                        gleich = true;
                        break;
                    end
                end
                if ~gleich;
                    iz = iz+1;
                    values_noneq(iz)=values(ii_jj);
                end
            end
            if (iz>1)
                   vmean = 0.2*vmean+0.8*sum(values_noneq(1:iz))/iz;
               values_noneq(1:iz)=values_noneq(1:iz)-vmean;
            end
           if (iz>1)
                  vv =sum(values_noneq(1:iz).^2)/iz;
                  if vv > 1.d-50 
                      vshape=-log(vv);
                  end
            end
        end 
    end

    %% Evacuated h-function
    function x = h_function(x_bar,s1,s2,x_p)
        
        H = x_bar .* (1.d0 - exp(-x_p.*s1)) + ...
            (1.0d0 - x_bar) .* exp(-(1.d0-x_p).*s2);              
        H0 = (1.d0 - x_bar) .* exp(-s2);
        H1 = x_bar .* exp(-s1);
        x = H + H1 .* x_p + H0 .* (x_p - 1.d0);
       
    end
    
    
    function [bestp,onep1,worstp] = must(IX,border_gute,n_par,ff2)
                                              
%
    iup=round(5.d0*(1.d0-ff2))+1;
    ilow=1;
    bestp=round(rand*(iup-ilow))+ ilow;
    worstp=-1 ;
    iup=15;
    ilow=0;
    while (worstp <= bestp) ||  (worstp > n_par)
     worstp=round(rand*(iup-ilow))+ ilow;  
     worstp=round(border_gute+(worstp-3));   
    end  
    iup=worstp-1;
    ilow=bestp+1;
    onep1=round(rand*(iup-ilow))+  ilow  ;    
    onep1= IX(onep1) ;
    bestp= IX(bestp) ;
    worstp= IX(worstp) ;
   end
    