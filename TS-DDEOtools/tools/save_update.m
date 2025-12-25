    % Judge Repeat Sample
    [~,ih,~]=intersect(hx,candidate_position,'rows'); % 
    if isempty(ih)~=1
        disp(['Sample Repeat and Delete it']);
    end

    % evaluate candidate
    candidate_fit=FUN(candidate_position);
    NFEs = NFEs + 1; 
    if show 
        disp([' candidate_fit=' num2str(candidate_fit) ' point = ' num2str(candidate_position) ]);
    end

    % save candidate into dataset, and sort dataset
    hx=[hx; candidate_position];  hf=[hf, candidate_fit];   % update history database 
    [hf,sidx]=sort(hf);                                         % sort point based on fitness, get point indexs
    hx=hx(sidx,:); 

    % update BestEvaluation and BestPoint
    if  candidate_fit <= hf(1) 
        BP = candidate_position;
        BE = candidate_fit;
        disp(['Best Cost(3) = ' num2str(BE) ' NFE=' num2str(NFEs)]);
    end

    % update CE for plotting
    CE(NFEs,:)=[NFEs,candidate_fit];
    if mod (NFEs,sn1)==0
        cs1=NFEs/sn1; gfs(1,cs1)=min(CE(1:NFEs,2));
    end