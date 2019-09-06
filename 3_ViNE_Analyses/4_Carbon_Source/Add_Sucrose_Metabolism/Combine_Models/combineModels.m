%% Load the models

% Load sucrose model
load('finalNodulatedPlant_sucrose.mat');
model_sucrose = finalNodulatedPlant;

% Load ViNE
load('../../finalNodulatedPlant.mat');
model_vine = finalNodulatedPlant;

%% Combine the models

% Get reaction formulas 
model_sucrose_formulas = printRxnFormula(model_sucrose, model_sucrose.rxns, false);
model_vine_formulas = printRxnFormula(model_vine, model_vine.rxns, false);

% Get the relevant information from the models
rxnNameList = vertcat(model_sucrose.rxnNames, model_vine.rxnNames);                        
rxnAbrList = vertcat(model_sucrose.rxns, model_vine.rxns);                        
rxnList = vertcat(model_sucrose, model_vine);
revFlagList = vertcat(model_sucrose.rev, model_vine.rev);
lowerBoundList = vertcat(model_sucrose.lb, model_vine.lb);
upperBoundList= vertcat(model_sucrose.ub, model_vine.ub);
subSystemList= vertcat(model_sucrose.subSystems, model_vine.subSystems);
grRuleList = vertcat(model_sucrose.grRules, model_vine.grRules);
geneNameList = vertcat(model_sucrose.genes, model_vine.genes);

% Combine the models
combinedModel = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, ...
    lowerBoundList, upperBoundList, subSystemList, grRuleList);
combinedModel.lb = transpose(combinedModel.lb);
combinedModel.ub = transpose(combinedModel.ub);
combinedModel.c = transpose(combinedModel.c);
combinedModel = changeObjective(combinedModel, 'Biomass');

% Save
save('combinedModel.mat', 'combinedModel');
save('allWorkspace.mat');
clear


