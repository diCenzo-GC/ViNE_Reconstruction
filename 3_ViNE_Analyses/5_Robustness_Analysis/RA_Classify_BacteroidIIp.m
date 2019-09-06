%% Start Cobra and set the 'ibm_cplex' Cobra solver
initCobraToolbox;
changeCobraSolver('ibm_cplex', 'all');

%% Bacteroid IIp reactions Classification (after adding bacteroid biomass sink) based on R.A.

% load the model and add bacteroid biomass sink
clearvars
modelFileName = load ('finalNodulatedPlant.mat');
model = modelFileName.finalNodulatedPlant; % load the model
model = addReaction(model,'BacteroidIId_SINK_Biomass',{'BacteroidIId_BIOMASS[c]'},[-1],0,0,1000000,0);
model = addReaction(model,'BacteroidIIp_SINK_Biomass',{'BacteroidIIp_BIOMASS[c]'},[-1],0,0,1000000,0);
model = addReaction(model,'BacteroidIZ_SINK_Biomass',{'BacteroidIZ_BIOMASS[c]'},[-1],0,0,1000000,0);

BacIIpRxnsPos = strfind(model.rxns, 'BacteroidIIp'); % find reactions whose name contains the string 'BacteroidIIp'
BacIIpRxnsPos(cellfun('isempty',BacIIpRxnsPos))={0}; % replace every empty cell with a zero
BacIIpRxnsPos = cell2mat(BacIIpRxnsPos); % transform cell array into double array
BacIIpRxns = model.rxns(BacIIpRxnsPos==1); % list Bacteroid IIp reactions
nRxns = length(BacIIpRxns); % number of BacteroidIIp reactions
sol = optimizeCbModel(model); % optimize model maximizing growth
nPoints = 41; % number of points for robustness analysis

% make a vector for every class
SteadyRxns = cell(nRxns,1);
SynergicRxns = cell(nRxns,1);
CompetingRxns = cell(nRxns,1);

for i=1:nRxns
    RxnStr = BacIIpRxns{i}; % reaction specified by the i-th position of BacIIpRxns
    disp(RxnStr);
    optFlux = sol.x(findRxnIDs(model, RxnStr)); % store optimal flux
    disp(strcat('optimal flux =',num2str(optFlux)))
    
    % store the original bounds
    minFlux = model.lb(findRxnIDs(model,RxnStr));
    maxFlux = model.ub(findRxnIDs(model,RxnStr));
    disp(strcat('lower bound =',num2str(minFlux),'...','upper bound =',num2str(maxFlux)))
    
    % walking away from the optimal flux, find those flux values which 
    % nullify growth
    
    % start with finding the zeroing-value standing at the right of the
    % optimal flux
    increment = 10^-5;
    z=1;
    while z==1
        forcedFlux = optFlux + increment;
        modelTest = changeRxnBounds(model,RxnStr,forcedFlux,'b');
        solTest = optimizeCbModel(modelTest);
        if forcedFlux >= maxFlux
            z=0;
            rightLimit=maxFlux;
        else
            if solTest.f>0
                z=1;
                increment=increment*10;
            else
                z=0;
                rightLimit=forcedFlux;
            end          
        end
    end
    
    % then find the zeroing-value standing at the left of the optimal flux
    modelTest = model;
    solTest = optimizeCbModel(model);
    increment = 10^-5;
    z=1;
    while z==1
        forcedFlux = optFlux - increment;
        modelTest = changeRxnBounds(model,RxnStr,forcedFlux,'b');
        solTest = optimizeCbModel(modelTest);
        if forcedFlux <= minFlux
            z=0;
            leftLimit=minFlux;
        else
            if solTest.f>0
                z=1;
                increment=increment*10;
            else
                z=0;
                leftLimit=forcedFlux;
            end
        end
    end    
    clearvars modelTest solTest
    disp(strcat('left limit =',num2str(leftLimit),'...','right limit =',num2str(rightLimit)))
    
    % perform Robustness Analysis between leftLimit and RightLimit
    RangeWidth = rightLimit - leftLimit;
    Step = RangeWidth / (nPoints-1);
    RobuArray = zeros(nPoints,2);
    for j=1:nPoints
        testedFlux = leftLimit + (j-1)*Step;
        model = changeRxnBounds(model,RxnStr,testedFlux,'b');
        solNew = optimizeCbModel(model);
        RobuArray(j,1) = testedFlux;
        RobuArray(j,2) = solNew.f;
        clearvars testedFlux solNew
    end
    ctrlFlux = RobuArray(:,1); % vector of control flux values
    objFlux = RobuArray(:,2); % vector or objective flux values
    
    % find first relevant point of the Robustness Analysis array
    firstPos = find(objFlux,1,'first');
    if firstPos>1
        firstPos=firstPos-1;
    end    
    
    %find last relevant point of the Robustness Analysis array
    lastPos = find(objFlux,1,'last');
    if lastPos<length(objFlux)
        lastPos = (lastPos+1);
    end 
    
    % resize RobuArray
    RobuArray((lastPos+1):(length(RobuArray)),:)=[]; % delete last useless rows from RobuArray
    RobuArray(1:(firstPos-1),:)=[]; % delete first useless rows from RobuArray
    ctrlFlux = RobuArray(:,1);
    objFlux = RobuArray(:,2);

    % add the optimal flux point to the array (if it isn't already there)
    if isempty(find(ctrlFlux==optFlux))
        model = changeRxnBounds(model,RxnStr,optFlux,'b');
        solNew = optimizeCbModel(model);
        RobuArray(length(RobuArray)+1,1) = optFlux;
        RobuArray(length(RobuArray),2) = solNew.f;
        clearvars testedFlux solNew
    end
    
    % add the zero-point to the array, but only if it is within the bounds
    if minFlux<0 && maxFlux>0 && isempty(find(ctrlFlux==0))
        model = changeRxnBounds(model,RxnStr,0,'b');
        solNew = optimizeCbModel(model);
        RobuArray(length(RobuArray)+1,1) = 0;
        RobuArray(length(RobuArray),2) = solNew.f;
        clearvars testedFlux solNew
    end
    
    % sort RobuArray's rows based on the asceding order of the elements of 
    % the ctrlFlux vector
    ctrlFlux = RobuArray(:,1);
    [ctrlFluxSorted, Indexes]=sort(ctrlFlux);
    RobuArray=RobuArray(Indexes,:);
    clearvars ctrlFluxSorted Indexes
    disp(RobuArray)
    ctrlFlux = RobuArray(:,1);
    objFlux = RobuArray(:,2);
    
    % find the objective value corresponding to the smallest ctrlFlux value
    alphaPos = find(abs(ctrlFlux) == min(abs(ctrlFlux)));
    alphaCtrl = ctrlFlux(alphaPos);
    alphaObj = objFlux(alphaPos);
    
    % find maximum and minimum of the objFlux vector
    maxObj = max(objFlux);
    minObj = min(objFlux);
    
    % get the index of the first empty cell, for each class vector
    SteadyRxnsIPos = find(cellfun('isempty',SteadyRxns),1);
    SynergicRxnsIPos = find(cellfun('isempty',SynergicRxns),1);
    CompetingRxnsIPos = find(cellfun('isempty',CompetingRxns),1);
    
    % LET'S CLASSIFY!
    if maxObj-minObj <= 0.2*maxObj
        SteadyRxns{SteadyRxnsIPos}=RxnStr;
        disp('STEADY')
    else 
        if alphaObj >= 0.95*maxObj
            CompetingRxns{CompetingRxnsIPos}=RxnStr;
            disp('COMPETING')
        else
            SynergicRxns{SynergicRxnsIPos}=RxnStr;
            disp('SYNERGIC')
        end
    end
    
    % restore original bounds and clear variables
    model.lb(findRxnIDs(model,RxnStr)) = minFlux;
    model.ub(findRxnIDs(model,RxnStr)) = maxFlux;
    clearvars -except model* Bac* sol nPoints *Rxns            
end
clearvars -except model* *Rxns nPoints

% concatenate the class vectors and delete empty rows
ListOfRxns = [SteadyRxns SynergicRxns CompetingRxns];
SteadyRxnsIPos = find(cellfun('isempty',SteadyRxns),1);
SynergicRxnsIPos = find(cellfun('isempty',SynergicRxns),1);
CompetingRxnsIPos = find(cellfun('isempty',CompetingRxns),1);
firstEmptyRow = max([SteadyRxnsIPos SynergicRxnsIPos CompetingRxnsIPos]);
ListOfRxns(firstEmptyRow:nRxns,:) = [];

% delete empty rows in each class vector, too
SteadyRxns(SteadyRxnsIPos:nRxns)=[];
SynergicRxns(SynergicRxnsIPos:nRxns)=[];
CompetingRxns(CompetingRxnsIPos:nRxns)=[];

% count reactions for every class vector
nSteadyRxns = length(SteadyRxns);
nSynergicRxns = length(SynergicRxns);
nCompetingRxns = length(CompetingRxns);

% visualize the list of sorted reactions
headers = {'Steady Reactions', 'Synergic Reactions', 'Competing Reactions'};
ListOfRxns = vertcat(headers, ListOfRxns);
disp(ListOfRxns);
ListOfRxns(1,:)=[];

% visualize the count of reactions for each class vector
CountOfRxns = [nSteadyRxns nSynergicRxns nCompetingRxns];
CountOfRxns = num2cell(CountOfRxns);
headers_bis = {'N° Steady Reactions','N° Synergic Reactions','N° Competing Reactions'};
CountOfRxns = vertcat(headers_bis, CountOfRxns);
disp(CountOfRxns);

% convert absolute frequencies to relative frequencies
CountOfRxns = [nSteadyRxns/nRxns nSynergicRxns/nRxns nCompetingRxns/nRxns];

% make a bar plot, and, for each class, also print the ratio between number  
% of reactions for the class and total number of reactions 
barLabels =  {strcat(num2str(nSteadyRxns),'/',num2str(nRxns),'=',num2str(nSteadyRxns/nRxns)),...
    strcat(num2str(nSynergicRxns),'/',num2str(nRxns),'=',num2str(nSynergicRxns/nRxns)),...
    strcat(num2str(nCompetingRxns),'/',num2str(nRxns),'=',num2str(nCompetingRxns/nRxns))};
figure
bar(CountOfRxns)
title('Frequencies of Reaction Classes in Bacteroids of Nodule Zone IIp')
axis = gca;
axis.XTickLabel = headers;
xPoints = axis.XTick;
text(xPoints, CountOfRxns, barLabels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
clearvars -except model* *Rxns nPoints
