%% Load and prepare the model

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Make version without oxygen limit
model2 = changeRxnBounds(model, 'TNIII_OXYGEN-MOLECULE', 1000000, 'u');

%% Prepar output variable

output_initialCosts = cell(3,7);
output_initialCosts{1,1} = 'Condition';
output_initialCosts{2,1} = 'Oxygen_limit';
output_initialCosts{3,1} = 'No_oxygen_limit';

%% Growth rate with nodule and nitrogen fixation

% Get solutions
sol = optimizeCbModel(model);
sol2 = optimizeCbModel(model2);

% Record results
output_initialCosts{1,2} = 'N2-fixation';
output_initialCosts{2,2} = sol.x(findRxnIDs(model, 'PlantBiomass'));
output_initialCosts{3,2} = sol2.x(findRxnIDs(model2, 'PlantBiomass'));

%% Growth rate with nodule and nitrogen fixation without bacteroid ATPM

% List of ATPM reactions
atpm = {'BacteroidIId_rxnNGAM'; 'BacteroidIIp_rxnNGAM'; 'BacteroidIZ_rxnNGAM'; 'BacteroidIII_rxnNGAM'};

% Remove ATPM
model3 = changeRxnBounds(model, atpm, 0, 'l');
model4 = changeRxnBounds(model2, atpm, 0, 'l');

% Get solutions
sol = optimizeCbModel(model3);
sol2 = optimizeCbModel(model4);

% Record results
output_initialCosts{1,3} = 'N2-fixation_noATPM';
output_initialCosts{2,3} = sol.x(findRxnIDs(model3, 'PlantBiomass'));
output_initialCosts{3,3} = sol2.x(findRxnIDs(model4, 'PlantBiomass'));

%% Growth rate without nodule, with nitrogen fixation, with ZIII bacteroid ATPM

% Remove nodule from first model
model3 = changeRxnMets(model, {'BiomassNodule[c]'}, {'BiomassNodule[c]'}, 'Biomass', 0);
noduleReactions = model3.rxns(strmatch('Nodule', model3.rxns));
noduleReactions2 = model3.rxns(strmatch('NoduleIII_', model3.rxns));
noduleReactions = setdiff(noduleReactions, noduleReactions2);
bacteroidReactions = model3.rxns(strmatch('Bacteroid', model3.rxns));
bacteroidReactions2 = model3.rxns(strmatch('BacteroidIII_', model3.rxns));
bacteroidReactions = setdiff(bacteroidReactions, bacteroidReactions2);
model3 = changeRxnBounds(model3, noduleReactions, 0, 'b');
model3 = changeRxnBounds(model3, bacteroidReactions, 0, 'b');

% Remove nodule from second model
model4 = changeRxnMets(model2, {'BiomassNodule[c]'}, {'BiomassNodule[c]'}, 'Biomass', 0);
noduleReactions = model4.rxns(strmatch('Nodule', model4.rxns));
noduleReactions2 = model4.rxns(strmatch('NoduleIII_', model4.rxns));
noduleReactions = setdiff(noduleReactions, noduleReactions2);
bacteroidReactions = model4.rxns(strmatch('Bacteroid', model4.rxns));
bacteroidReactions2 = model4.rxns(strmatch('BacteroidIII_', model4.rxns));
bacteroidReactions = setdiff(bacteroidReactions, bacteroidReactions2);
model4 = changeRxnBounds(model4, noduleReactions, 0, 'b');
model4 = changeRxnBounds(model4, bacteroidReactions, 0, 'b');

% Get solutions
sol = optimizeCbModel(model3);
sol2 = optimizeCbModel(model4);

% Record results
output_initialCosts{1,4} = 'N2-fixation_noNodule_wATPM';
output_initialCosts{2,4} = sol.x(findRxnIDs(model3, 'PlantBiomass'));
output_initialCosts{3,4} = sol2.x(findRxnIDs(model4, 'PlantBiomass'));

%% Growth rate without nodule, with nitrogen fixation, without ZIII bacteroid ATPM

% Remove nodule from first model
model3 = changeRxnMets(model, {'BiomassNodule[c]'}, {'BiomassNodule[c]'}, 'Biomass', 0);
noduleReactions = model3.rxns(strmatch('Nodule', model3.rxns));
noduleReactions2 = model3.rxns(strmatch('NoduleIII_', model3.rxns));
noduleReactions = setdiff(noduleReactions, noduleReactions2);
bacteroidReactions = model3.rxns(strmatch('Bacteroid', model3.rxns));
bacteroidReactions2 = model3.rxns(strmatch('BacteroidIII_', model3.rxns));
bacteroidReactions = setdiff(bacteroidReactions, bacteroidReactions2);
bacteroidReactions = vertcat(bacteroidReactions, {'BacteroidIII_rxnNGAM'});
model3 = changeRxnBounds(model3, noduleReactions, 0, 'b');
model3 = changeRxnBounds(model3, bacteroidReactions, 0, 'b');

% Remove nodule from second model
model4 = changeRxnMets(model2, {'BiomassNodule[c]'}, {'BiomassNodule[c]'}, 'Biomass', 0);
noduleReactions = model4.rxns(strmatch('Nodule', model4.rxns));
noduleReactions2 = model4.rxns(strmatch('NoduleIII_', model4.rxns));
noduleReactions = setdiff(noduleReactions, noduleReactions2);
bacteroidReactions = model4.rxns(strmatch('Bacteroid', model4.rxns));
bacteroidReactions2 = model4.rxns(strmatch('BacteroidIII_', model4.rxns));
bacteroidReactions = setdiff(bacteroidReactions, bacteroidReactions2);
bacteroidReactions = vertcat(bacteroidReactions, {'BacteroidIII_rxnNGAM'});
model4 = changeRxnBounds(model4, noduleReactions, 0, 'b');
model4 = changeRxnBounds(model4, bacteroidReactions, 0, 'b');

% Get solutions
sol = optimizeCbModel(model3);
sol2 = optimizeCbModel(model4);

% Record results
output_initialCosts{1,5} = 'N2-fixation_noNodule_woATPM';
output_initialCosts{2,5} = sol.x(findRxnIDs(model3, 'PlantBiomass'));
output_initialCosts{3,5} = sol2.x(findRxnIDs(model4, 'PlantBiomass'));

%% Growth with ammonium with nodule

% Add ammonium uptake
model3 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model4 = changeRxnBounds(model2, 'Root_TEC_AMMONIUM', 1000000, 'u');

% Get solutions
sol = optimizeCbModel(model3);
sol2 = optimizeCbModel(model4);

% Record results
output_initialCosts{1,6} = 'Ammonium';
output_initialCosts{2,6} = sol.x(findRxnIDs(model3, 'PlantBiomass'));
output_initialCosts{3,6} = sol2.x(findRxnIDs(model4, 'PlantBiomass'));

%% Growth with ammonium without nodule

% Add ammonium uptake
model3 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model4 = changeRxnBounds(model2, 'Root_TEC_AMMONIUM', 1000000, 'u');

% Remove nodule from first model
model3 = changeRxnMets(model3, {'BiomassNodule[c]'}, {'BiomassNodule[c]'}, 'Biomass', 0);
noduleReactions = model3.rxns(strmatch('Nodule', model3.rxns));
bacteroidReactions = model3.rxns(strmatch('Bacteroid', model3.rxns));
model3 = changeRxnBounds(model3, noduleReactions, 0, 'b');
model3 = changeRxnBounds(model3, bacteroidReactions, 0, 'b');

% Remove nodule from second model
model4 = changeRxnMets(model4, {'BiomassNodule[c]'}, {'BiomassNodule[c]'}, 'Biomass', 0);
noduleReactions = model4.rxns(strmatch('Nodule', model4.rxns));
bacteroidReactions = model4.rxns(strmatch('Bacteroid', model4.rxns));
model4 = changeRxnBounds(model4, noduleReactions, 0, 'b');
model4 = changeRxnBounds(model4, bacteroidReactions, 0, 'b');

% Get solutions
sol = optimizeCbModel(model3);
sol2 = optimizeCbModel(model4);

% Record results
output_initialCosts{1,7} = 'Ammonium_noNodule';
output_initialCosts{2,7} = sol.x(findRxnIDs(model3, 'PlantBiomass'));
output_initialCosts{3,7} = sol2.x(findRxnIDs(model4, 'PlantBiomass'));

%% Test effect of NGAM

% Create and add a new bacteroid NGAM reaction
newNGAM = '5670 BacteroidIId_MNXM2[c] + 5670 BacteroidIIp_MNXM2[c] + 630 BacteroidIZ_MNXM2[c] + 6300 BacteroidIII_MNXM2[c] + 5670 BacteroidIId_MNXM3[c] + 5670 BacteroidIIp_MNXM3[c] + 630 BacteroidIZ_MNXM3[c] + 6300 BacteroidIII_MNXM3[c] -> 5670 BacteroidIId_MNXM1[c] + 5670 BacteroidIIp_MNXM1[c] + 630 BacteroidIZ_MNXM1[c] + 6300 BacteroidIII_MNXM1[c] + 5670 BacteroidIId_MNXM9[c] + 5670 BacteroidIIp_MNXM9[c] + 630 BacteroidIZ_MNXM9[c] + 6300 BacteroidIII_MNXM9[c] + 5670 BacteroidIId_MNXM7[c] + 5670 BacteroidIIp_MNXM7[c] + 630 BacteroidIZ_MNXM7[c] + 6300 BacteroidIII_MNXM7[c]';
model3 = addReaction(model, 'combinedNGAM', newNGAM, [], 0, 0, 1000000, 0);
model4 = addReaction(model2, 'combinedNGAM', newNGAM, [], 0, 0, 1000000, 0);

% Turn off old NGAM reactions
model3.lb(findRxnIDs(model3, 'BacteroidIId_rxnNGAM')) = 0;
model3.lb(findRxnIDs(model3, 'BacteroidIIp_rxnNGAM')) = 0;
model3.lb(findRxnIDs(model3, 'BacteroidIZ_rxnNGAM')) = 0;
model3.lb(findRxnIDs(model3, 'BacteroidIII_rxnNGAM')) = 0;
model4.lb(findRxnIDs(model4, 'BacteroidIId_rxnNGAM')) = 0;
model4.lb(findRxnIDs(model4, 'BacteroidIIp_rxnNGAM')) = 0;
model4.lb(findRxnIDs(model4, 'BacteroidIZ_rxnNGAM')) = 0;
model4.lb(findRxnIDs(model4, 'BacteroidIII_rxnNGAM')) = 0;

% Run the pareto analyses
output_effectNGAM_1 = paretoFrontier(model3, 'combinedNGAM', 'Biomass', 100);
output_effectNGAM_2 = paretoFrontier(model4, 'combinedNGAM', 'Biomass', 100);

%%  Clear workspace and save

clearvars -except output_*
