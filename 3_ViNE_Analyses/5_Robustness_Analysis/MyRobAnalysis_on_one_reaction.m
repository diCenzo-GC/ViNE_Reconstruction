%% Start Cobra and set the 'ibm_cplex' Cobra solver
initCobraToolbox;
changeCobraSolver('ibm_cplex', 'all');    

%%  My Robustness Analysis on one reaction 

clearvars
modelFileName = load ('finalNodulatedPlant.mat');
model = modelFileName.finalNodulatedPlant; % load the model
figure
sol = optimizeCbModel(model); % optimize model maximizing growth
nPoints=41;
    
RxnStr = 'BacteroidIId_MNXR97897';
disp(RxnStr);
RxnAbbr = strrep(RxnStr,'BacteroidIId_','BactIId-');
 
% store the original bounds and the optimal flux
optFlux = sol.x(findRxnIDs(model, RxnStr));
disp(strcat('optimal flux =',num2str(optFlux)))
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
ctrlFlux = RobuArray(:,1);
objFlux = RobuArray(:,2); % vector of objective flux values
    
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

% add the zero-point to the array, but only if it is within the bounds
if minFlux<0 && maxFlux>0 && isempty(find(ctrlFlux==0))
    model = changeRxnBounds(model,RxnStr,0,'b');
    solNew = optimizeCbModel(model);
    RobuArray(length(RobuArray)+1,1) = 0;
    RobuArray(length(RobuArray),2) = solNew.f;
    clearvars testedFlux solNew
end

% add the optimal flux point to the array
if isempty(find(ctrlFlux==optFlux))
    model = changeRxnBounds(model,RxnStr,optFlux,'b');
    solNew = optimizeCbModel(model);
    RobuArray(length(RobuArray)+1,1) = optFlux;
    RobuArray(length(RobuArray),2) = solNew.f;
    clearvars testedFlux solNew
end

% sort RobuArray's rows based on the asceding order of the elements of the
% first row
ctrlFlux = RobuArray(:,1);
[ctrlFluxSorted, Indexes]=sort(ctrlFlux);
RobuArray=RobuArray(Indexes,:);
clearvars ctrlFluxSorted Indexes
disp(RobuArray)

% make a vector out of each column of the array
ctrlFlux = RobuArray(:,1);
objFlux = RobuArray(:,2);

% find the objective value corresponding to the smallest ctrlFlux value
alphaPos = find(abs(ctrlFlux) == min(abs(ctrlFlux)));
alphaCtrl = ctrlFlux(alphaPos);
alphaObj = objFlux(alphaPos);

% find the objective value corresponding to the biggest ctrlFlux value
omegaPos = find(abs(ctrlFlux)==max(abs(ctrlFlux)));
omegaCtrl = ctrlFlux(omegaPos);
omegaObj = objFlux(omegaPos);

% find the maximum and the minimum of objFlux vector
maxObj = max(objFlux);
minObj = min(objFlux);
    
% plot the Robustness Analysis results
plot(ctrlFlux, objFlux);
hold on
xlabel('(\mumol/h/g_{DW})','FontSize',8);
ylabel('Objective');
title(RxnAbbr,'FontSize',8);
ylim([0 2])
line([0 0],[0 2],'Color','red')
hold off

    % LET'S CLASSIFY!
    if maxObj-minObj <= 0.1*maxObj
        disp('STEADY')
    else 
        if abs(mean(ctrlFlux)) <= 10^-4 && alphaObj >= 0.95*maxObj
            disp('TOXIC')
        else
            if alphaObj >= 0.95*maxObj && sign(optFlux)+sign(omegaCtrl)~=0
                disp('COMPETING')
            else
                disp('SYNERGIC')
            end
        end
    end

% restore original bounds
model.lb(findRxnIDs(model,RxnStr)) = minFlux;
model.ub(findRxnIDs(model,RxnStr)) = maxFlux;  
    
    