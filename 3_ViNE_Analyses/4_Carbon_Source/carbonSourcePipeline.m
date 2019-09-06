%% Make the combined model

% Change directory
cd Add_Sucrose_Metabolism

% Run pipeline
pipeline;
clear

% Change directory
cd ..

%% Precursor model growth with sucrose

% Prepare output variable
output_carbonUsage = cell(5,2);
output_carbonUsage{1,1} = 'Carbon_Source';
output_carbonUsage{1,2} = 'Biomass';
output_carbonUsage{2,1} = 'Precursor_+Sucrose';
output_carbonUsage{3,1} = 'Precursor_-Sucrose';
output_carbonUsage{4,1} = 'Reduced_+Sucrose';
output_carbonUsage{5,1} = 'Reduced_-Sucrose';

% Load the model
load('Add_Sucrose_Metabolism/Constrain_Nodule/nodulatedPlant.mat');
model = nodulatedPlant;

% Increase flux
model.ub = model.ub * 1000;
model.lb = model.lb * 1000;

% Set O2 consumption limit of nodule zone III
model.ub(findRxnIDs(model, 'TRN_OXYGEN-MOLECULE')) = 1000 * 649 * 0.02;
sol = optimizeCbModel(model);
model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) = ...
    sol.x(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE'));
model.ub(findRxnIDs(model, 'TRN_OXYGEN-MOLECULE')) = 1000000;

% Remove unwanted NoduleIII_EXCT reactions
toRemove = {'NoduleIII_EXCT_for_MNXM621_e0';...
    'NoduleIII_EXCT_for_MNXM1503_e0';'NoduleIII_EXCT_for_MNXM198_e0';...
    'NoduleIII_EXCT_for_MNXM165_e0';'NoduleIII_EXCT_for_MNXM468_e0';...
    'NoduleIII_EXCT_for_MNXM615_e0'};
model = tncore_remove_reactions(model, toRemove);
sol = optimizeCbModel(model);
reactions = model.rxns(strmatch('NoduleIII_EXCT_for', model.rxns));
toRemove = {};
for n = 1:length(reactions)
    pos = findRxnIDs(model, reactions{n});
    if abs(sol.x(pos)) < 0.00001
        toRemove = vertcat(toRemove, reactions{n});
    end
end
toRemove = setdiff(toRemove, {'NoduleIII_EXCT_for_MNXM167_e0'});
model = tncore_remove_reactions(model, toRemove);

% Record growth rate 
sol = optimizeCbModel(model);
output_carbonUsage{2,2} = sol.f;

%% Precursor model growth with C4-dicarboxylates

% Load the model
load('Add_Sucrose_Metabolism/Constrain_Nodule/nodulatedPlant.mat')
model = nodulatedPlant;

% Increase flux
model.ub = model.ub * 1000;
model.lb = model.lb * 1000;

% Set O2 consumption limit of nodule zone III
model.ub(findRxnIDs(model, 'TRN_OXYGEN-MOLECULE')) = 1000 * 649 * 0.02;
sol = optimizeCbModel(model);
model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) = ...
    sol.x(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE'));
model.ub(findRxnIDs(model, 'TRN_OXYGEN-MOLECULE')) = 1000000;

% Remove unwanted NoduleIII_EXCT reactions
toRemove = {'NoduleIII_EXCT_for_MNXM167_e0';'NoduleIII_EXCT_for_MNXM621_e0';...
    'NoduleIII_EXCT_for_MNXM1503_e0';'NoduleIII_EXCT_for_MNXM198_e0';...
    'NoduleIII_EXCT_for_MNXM165_e0';'NoduleIII_EXCT_for_MNXM468_e0';...
    'NoduleIII_EXCT_for_MNXM615_e0'};
model = tncore_remove_reactions(model, toRemove);
sol = optimizeCbModel(model);
reactions = model.rxns(strmatch('NoduleIII_EXCT_for', model.rxns));
toRemove = {};
for n = 1:length(reactions)
    pos = findRxnIDs(model, reactions{n});
    if abs(sol.x(pos)) < 0.00001
        toRemove = vertcat(toRemove, reactions{n});
    end
end
toRemove = setdiff(toRemove, {'NoduleIII_EXCT_for_MNXM25_e0';'NoduleIII_EXCT_for_MNXM98_e0';'NoduleIII_EXCT_for_MNXM93_e0'});
model = tncore_remove_reactions(model, toRemove);

% Record growth rate (oxygen limit, proton transfer)
sol = optimizeCbModel(model);
output_carbonUsage{3,2} = sol.f;

%% Compare sucrose versus C4-dicarboxylates

% Load the model
load('Add_Sucrose_Metabolism/combinedModel.mat');
model = combinedModel;

% Make version without sucrose usage
model2 = changeRxnBounds(model, 'NoduleIII_EXCT_for_MNXM167_e0', 0, 'b');

% Get FBA solutions
sol = optimizeCbModel(model);
sol2 = optimizeCbModel(model2);

% Store output
output_carbonUsage{4,2} = sol.f;
output_carbonUsage{5,2} = sol2.f;

% Clean workspace
clearvars -except output_*

%% Effect of mitochondrial terminal oxidase limitation

% Load the model
load('Add_Sucrose_Metabolism/combinedModel.mat');
model = combinedModel;

% Remove oxygen limit
model = changeRxnBounds(model, 'TNIII_OXYGEN-MOLECULE', 1000000, 'u');

% Get maximal rate of terminal oxidase
sol = optimizeCbModel(model);
O2max = sol.x(findRxnIDs(model, 'NoduleIII_MNXR140037_M'));

% Make an output variable
output_terminalOxidase = cell(101,7);

% Test effect of terminal oxidase
for n = 0:100
    testModel = changeRxnBounds(model, 'NoduleIII_MNXR140037_M', O2max * (n / 100), 'u');
    sol = optimizeCbModel(testModel);
    output_terminalOxidase{n+1,1} = O2max * (n / 100);
    output_terminalOxidase{n+1,2} = sol.f;
    output_terminalOxidase{n+1,3} = sol.x(findRxnIDs(model, 'NoduleIII_EXCT_for_MNXM167_e0'));
    output_terminalOxidase{n+1,4} = sol.x(findRxnIDs(model, 'NoduleIII_EXCT_for_MNXM25_e0')) + ...
        sol.x(findRxnIDs(model, 'NoduleIII_EXCT_for_MNXM98_e0')) + ...
        sol.x(findRxnIDs(model, 'NoduleIII_EXCT_for_MNXM93_e0'));
    output_terminalOxidase{n+1,5} = sol.x(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE'));
    output_terminalOxidase{n+1,6} = sol.x(findRxnIDs(model, 'NoduleIII_MNXR140037_M'));
    output_terminalOxidase{n+1,7} = sol.x(findRxnIDs(model, 'BacteroidIII_MNXR138955'));
end

% Add headers
headers = {'Mito_oxidase_limit', 'Biomass', 'Sucrose', ...
    'Dicarboxylates', 'Nodule_Oxygen', 'Mito_oxidase', 'Bacteroid_oxidase'};
output_terminalOxidase = vertcat(headers, output_terminalOxidase);

% Save and clear
clearvars -except output_*
save('carbonSourceOutput.mat');
clear
