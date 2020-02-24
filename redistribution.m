function lambdanew=redistribution(lambdacurrent,lambdainit,OS,OSnew,redistr)
% updates the firing propensity when an origin fires or is passively replicated
if redistr==1

    additional_lambda=sum(lambdainit); %total firing propensity 
    lambdanew=lambdacurrent; %initialize
    [forks, ori2alter]=find(OSnew==0); %which origins are in PreR

    active_lambda=0; %initialize

    for j=1:length(forks)
        active_lambda=active_lambda+lambdainit(ori2alter(j));
    end

    if active_lambda~=0
        coef=additional_lambda/active_lambda;
    elseif active_lambda==0
        coef=0;
    end

    TFP=0;
    for j=1:length(forks)
        lambdanew(forks(j),ori2alter(j)) = lambdainit(ori2alter(j))*coef;
        TFP=TFP+lambdanew(forks(j),ori2alter(j));
    end

elseif redistr==0
    
    lambdanew=lambdacurrent; %initialize
    [forks, ori2alter]=find(OSnew==0); %which origins are in PreR
    
    for j=1:length(forks)
        lambdanew(forks(j),ori2alter(j)) = lambdainit(ori2alter(j));
    end
    
end