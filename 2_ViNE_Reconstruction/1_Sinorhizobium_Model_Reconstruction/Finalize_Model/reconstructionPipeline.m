%% Initialize the system

% MATLAB R2016b
% Gurobi version 7.0.1
% libSBML version 5.13.0
% SBMLToolbox version 4.1.0
% COBRA Toolbox downloaded May 12, 2017
% Tn-Core Toolbox version 2.2

% Set up the system
clear all
addpath(genpath('../../Software/cobratoolbox/'));
rmpath(genpath('../../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../Software/SBMLToolbox-4.1.0/'));
addpath(genpath('../../Software/Tn-Core-v2.2/'));
initCobraToolbox;
changeCobraSolver('gurobi', 'all');

%% Finalize the model

% Import model and check growth
model3 = readCbModel('expanded_model.xml');
model3 = changeObjective(model3, 'rxnBIOMASS');
optimizeCbModel(model3)

% Remove dead-ends and export
model3 = tncore_deadends(model3);
tncore_export(model3, model3);

% Do stuff manually, including adding a NGAM reaction

%% Remove unwanted genes or GPRs

% Load the model
model = readCbModel('finalModel.xml');

% Remove the NonEssential gene
model2 = deleteModelGenes(model, 'NonEssential');
model2 = tncore_delete(model2);

% Find gene index number for Unknown
model3 = model2;
unknownID = findGeneIDs(model3,'Unknown');

% Replace the unknowns that are not needed
A = [' | x(' num2str(unknownID) ')'];
B = ['x(' num2str(unknownID) ') | '];
model3.grRules = strrep(model3.grRules,' or Unknown','');
model3.grRules = strrep(model3.grRules,'Unknown or ','');
model3.rxnNotes = strrep(model3.rxnNotes,' or Unknown','');
model3.rxnNotes = strrep(model3.rxnNotes,'Unknown or ','');
model3.rules = strrep(model3.rules,A,'');
model3.rules = strrep(model3.rules,B,'');

% Fix the rxnGeneMat
rxnGeneMatNewFull = cell(length(model3.rxns),length(model3.genes));
for n = 1:length(model3.rxns)
    if ~isempty(model3.rules{n})
        rulesTemp = strrep(model3.rules{n}, '&', '|');
        rules = strsplit(rulesTemp, '|');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            rxnGeneMatNewFull{n,str2num(rules{m})} = 1;
        end
    end
end
temp = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(temp) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
model3.rxnGeneMat = rxnGeneMatNew;

% Save the model
finalModel = model3;
save('finalModel_intermediate.mat', 'finalModel');
save('allWorkspace.mat');
clear;

%% Fix energy

% No effect on growth was seen with or without flux through the NGAM
% reaction. So, there was an issue with the model somewhere. Worked on
% fixing this manually in the directory called Test_Energy. I was able to
% fix this error.

%% Remove duplicate reactions/genes
% Identify reactions added during the automated step with at least one
% associated gene also associated with a reaction that was manually added
% to the model. Remove these reactions, where appropriate.

% Load the models
model = readCbModel('finalModel_energyFixed.xml');
tempModel = readCbModel('model_noAutomated.xml');
model = tncore_fix(model);

% Find automated reactions with genes from manual stages
[~, reactionsAll] = findRxnsFromGenes(model, tempModel.genes, 0, 1);
reactionsAll = reactionsManual(:,1);
reactionsMan = setdiff(reactionsAll, tempModel.rxns);

% Save
save('reactionsToDelete.mat');
clear;
% Now manually check the list of genes/reactions and remove as appropriate

%% Final modifications

% After removing the reactions from above, duplicate reactions were
% identified and either removed (if second reaction added by Tn-Core) or
% the GPRs merged (if both rxns were manually added). Any new dead-ends
% were removed, unnecessary exchange reactions were removed, and unused
% metabolites were deleted

% Import final model
iGD1348 = readCbModel('finalModel.xml');

% Save full final model
save('iGD1348.mat', 'iGD1348');

% Remove the NonEssential gene
model2 = deleteModelGenes(iGD1348, 'NonEssential');
model2 = tncore_delete(model2);

% Find gene index number for Unknown
model3 = model2;
unknownID = findGeneIDs(model3,'Unknown');

% Replace the unknowns that are not needed
A = [' | x(' num2str(unknownID) ')'];
B = ['x(' num2str(unknownID) ') | '];
model3.grRules = strrep(model3.grRules,' or Unknown','');
model3.grRules = strrep(model3.grRules,'Unknown or ','');
model3.rxnNotes = strrep(model3.rxnNotes,' or Unknown','');
model3.rxnNotes = strrep(model3.rxnNotes,'Unknown or ','');
model3.rules = strrep(model3.rules,A,'');
model3.rules = strrep(model3.rules,B,'');

% Fix the rxnGeneMat
rxnGeneMatNewFull = cell(length(model3.rxns),length(model3.genes));
for n = 1:length(model3.rxns)
    if ~isempty(model3.rules{n})
        rulesTemp = strrep(model3.rules{n}, '&', '|');
        rules = strsplit(rulesTemp, '|');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            rxnGeneMatNewFull{n,str2num(rules{m})} = 1;
        end
    end
end
temp = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(temp) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
model3.rxnGeneMat = rxnGeneMatNew;
iGD1348 = model3;

% Save model in MAT format
save('iGD1348_working.mat', 'iGD1348');
clear;
