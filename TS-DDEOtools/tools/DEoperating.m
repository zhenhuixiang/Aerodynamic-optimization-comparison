% DEoperating use NP_last of population generate NP_next of offspring
function [U] = DEoperating(P,NP_last,NP_next,Dim,hisx,F,CR,UB,LB)
    NP = NP_last;
    a = randperm(NP_last);
    for m=1:NP_next  
        if m > NP_last
            i = randi(NP_last);
        else
            i = a(m);
        end
        %% mutation
        k0=randi([1,NP]);
        while(k0==i)
            k0=randi([1,NP]);   
        end
        P1=P(k0,:);
        k1=randi([1,NP]);
        while(k1==i||k1==k0)
            k1=randi([1,NP]);
        end
        P2=P(k1,:);
        k2=randi([1,NP]);
        while(k2==i||k2==k1||k2==k0)
            k2=randi([1,NP]);
        end
        P3=P(k2,:);

        Xpbest = hisx(1,:);
        V(m,:)=P1+F.*(P2-P3); 
%         V(m,:)=P1+F.*(P2-P3)+F.*(Xpbest-P1);   

        %% bound
        for j=1:Dim
          if (V(m,j)>UB(i,j)||V(m,j)<LB(i,j))
             V(m,j)=LB(i,j)+rand*(UB(i,j)-LB(i,j));         
          end
        end

        %% crossover
        jrand=randi([1,Dim]); 
        for j=1:Dim
            k3=rand;
            if(k3<=CR||j==jrand)
                U(m,j)=V(m,j);
            else
                U(m,j)=P(i,j);      % 子代种群
            end
        end
    end
end