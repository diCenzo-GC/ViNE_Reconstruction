%% Load the models

load('noduleModel.mat');
load('plantModel.mat');
load('medicagoModel');
load('melilotiModel');

%% Prepare the input models

prepareModels;

%% Make four nodule models, one per zone, and a zone I from medicago model

zoneI = medicagoModel;
zoneIId = noduleModel;
zoneIIp = noduleModel;
zoneIZ = noduleModel;
zoneIII = noduleModel;

%% Rename reactions, metabolites, and genes in the nodule zones

renameFeatures;

%% Delete all exchange reactions from the nodule zones

deleteExchange;

%% Add the nodule to the plant

% Get reaction formulas for all models
zoneIFormulas = printRxnFormula(zoneI_noEx, zoneI_noEx.rxns, false);
zoneIIdFormulas = printRxnFormula(zoneIId_noEx, zoneIId_noEx.rxns, false);
zoneIIpFormulas = printRxnFormula(zoneIIp_noEx, zoneIIp_noEx.rxns, false);
zoneIZFormulas = printRxnFormula(zoneIZ_noEx, zoneIZ_noEx.rxns, false);
zoneIIIFormulas = printRxnFormula(zoneIII_noEx, zoneIII_noEx.rxns, false);
plantFormulas = printRxnFormula(plantModel, plantModel.rxns, false);

% Get the relevant information from all models
rxnNameList = vertcat(plantModel.rxnNames, zoneI_noEx.rxnNames, zoneIId_noEx.rxnNames, ...
    zoneIIp_noEx.rxnNames, zoneIZ_noEx.rxnNames, zoneIII_noEx.rxnNames);                        
rxnAbrList = vertcat(plantModel.rxns, zoneI_noEx.rxns, zoneIId_noEx.rxns, ...
    zoneIIp_noEx.rxns, zoneIZ_noEx.rxns, zoneIII_noEx.rxns);                        
rxnList = vertcat(plantFormulas, zoneIFormulas, zoneIIdFormulas, zoneIIpFormulas, ...
    zoneIZFormulas, zoneIIIFormulas);
revFlagList = vertcat(plantModel.rev, zoneI_noEx.rev, zoneIId_noEx.rev, ...
    zoneIIp_noEx.rev, zoneIZ_noEx.rev, zoneIII_noEx.rev);
lowerBoundList = vertcat(plantModel.lb, zoneI_noEx.lb, zoneIId_noEx.lb, ...
    zoneIIp_noEx.lb, zoneIZ_noEx.lb, zoneIII_noEx.lb);
upperBoundList= vertcat(plantModel.ub, zoneI_noEx.ub, zoneIId_noEx.ub, ...
    zoneIIp_noEx.ub, zoneIZ_noEx.ub, zoneIII_noEx.ub);
subSystemList= vertcat(plantModel.subSystems, zoneI_noEx.subSystems, zoneIId_noEx.subSystems, ...
    zoneIIp_noEx.subSystems, zoneIZ_noEx.subSystems, zoneIII_noEx.subSystems);
grRuleList = vertcat(plantModel.grRules, zoneI_noEx.grRules, zoneIId_noEx.grRules, ...
    zoneIIp_noEx.grRules, zoneIZ_noEx.grRules, zoneIII_noEx.grRules);
geneNameList = vertcat(plantModel.genes, zoneI_noEx.genes, zoneIId_noEx.genes, ...
    zoneIIp_noEx.genes, zoneIZ_noEx.genes, zoneIII_noEx.genes);

% Combine the models
nodulatedPlantDisconnected = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, ...
    lowerBoundList, upperBoundList, subSystemList, grRuleList);
nodulatedPlantDisconnected.lb = transpose(nodulatedPlantDisconnected.lb);
nodulatedPlantDisconnected.ub = transpose(nodulatedPlantDisconnected.ub);
nodulatedPlantDisconnected.c = transpose(nodulatedPlantDisconnected.c);
nodulatedPlantDisconnected = changeObjective(nodulatedPlantDisconnected, 'Biomass');

%% Connect nodule metabolism to plant metabolism

% Transfer of protons from plant to bacteroid in zone III
nodulatedPlantDisconnected = addReaction(nodulatedPlantDisconnected, 'BacteroidIII_protConv', ...
    'BacteroidIII_cpd00067[e] -> BacteroidIII_cpd00067[p]', [], 0, 0, 1000000, 0);

% Add proton exchange between symbiosome and periplasm
nodulatedPlantDisconnected = addReaction(nodulatedPlantDisconnected, 'ZoneIId_ProtonConv', ...
    {'BacteroidIId_cpd00067[p]', 'BacteroidIId_cpd00067[e]'}, ...
    [-1 1], 0, 0, 1000, 0);
nodulatedPlantDisconnected = addReaction(nodulatedPlantDisconnected, 'ZoneIIp_ProtonConv', ...
    {'BacteroidIIp_cpd00067[p]', 'BacteroidIIp_cpd00067[e]'}, ...
    [-1 1], 0, 0, 1000, 0);
nodulatedPlantDisconnected = addReaction(nodulatedPlantDisconnected, 'ZoneIZ_ProtonConv', ...
    {'BacteroidIZ_cpd00067[p]', 'BacteroidIZ_cpd00067[e]'}, ...
    [-1 1], 0, 0, 1000, 0);
nodulatedPlantDisconnected = addReaction(nodulatedPlantDisconnected, 'ZoneIII_ProtonConv', ...
    {'BacteroidIII_cpd00067[p]', 'BacteroidIII_cpd00067[e]'}, ...
    [-1 1], 0, 0, 1000, 0);

% Update symport reactions
transportRxns = findRxnsFromMets(nodulatedPlantDisconnected, 'BacteroidIII_cpd00067[p]');
toKeep = {'BacteroidIII_rxn10122'; 'BacteroidIII_rxn12750'; 'BacteroidIII_rxn10042'; ...
    'BacteroidIII_rxn10121'; 'BacteroidIII_rxn10120'; 'BacteroidIII_rxn05890'; ...
    'BacteroidIII_rxn01806'; 'BacteroidIII_rxn00001a'; 'BacteroidIII_rxn11268'; ...
    'ZoneIII_ProtonConv'};
transportRxns = setdiff(transportRxns, toKeep);
transportRxnsTemp = findRxnsFromMets(nodulatedPlantDisconnected, 'BacteroidIZ_cpd00067[p]');
toKeep = {'BacteroidIZ_rxn10122'; 'BacteroidIZ_rxn12750'; 'BacteroidIZ_rxn10042'; ...
    'BacteroidIZ_rxn10121'; 'BacteroidIZ_rxn10120'; 'BacteroidIZ_rxn05890'; ...
    'BacteroidIZ_rxn01806'; 'BacteroidIZ_rxn00001a'; 'BacteroidIZ_rxn11268'; ...
    'ZoneIII_ProtonConv'};
transportRxnsTemp = setdiff(transportRxnsTemp, toKeep);
transportRxns = vertcat(transportRxns, transportRxnsTemp);
transportRxnsTemp = findRxnsFromMets(nodulatedPlantDisconnected, 'BacteroidIId_cpd00067[p]');
toKeep = {'BacteroidIId_rxn10122'; 'BacteroidIId_rxn12750'; 'BacteroidIId_rxn10042'; ...
    'BacteroidIId_rxn10121'; 'BacteroidIId_rxn10120'; 'BacteroidIId_rxn05890'; ...
    'BacteroidIId_rxn01806'; 'BacteroidIId_rxn00001a'; 'BacteroidIId_rxn11268'; ...
    'ZoneIII_ProtonConv'};
transportRxnsTemp = setdiff(transportRxnsTemp, toKeep);
transportRxns = vertcat(transportRxns, transportRxnsTemp);
transportRxnsTemp = findRxnsFromMets(nodulatedPlantDisconnected, 'BacteroidIIp_cpd00067[p]');
toKeep = {'BacteroidIIp_rxn10122'; 'BacteroidIIp_rxn12750'; 'BacteroidIIp_rxn10042'; ...
    'BacteroidIIp_rxn10121'; 'BacteroidIIp_rxn10120'; 'BacteroidIIp_rxn05890'; ...
    'BacteroidIIp_rxn01806'; 'BacteroidIIp_rxn00001a'; 'BacteroidIIp_rxn11268'; ...
    'ZoneIII_ProtonConv'};
transportRxnsTemp = setdiff(transportRxnsTemp, toKeep);
transportRxns = vertcat(transportRxns, transportRxnsTemp);
transportRxnIDs = findRxnIDs(nodulatedPlantDisconnected, transportRxns);
nodulatedPlantDisconnected.rev(transportRxnIDs) = 0;
nodulatedPlantDisconnected.lb(transportRxnIDs) = 0;

% Add transfer reactions
connectNodule;

%% Set the objective function

% Make reactions for zone specific biomasses
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'ZoneI_Biomass', ...
    {'NoduleI_BiomassRoot[c]', 'BiomassZoneI[c]'}, ...
    [-1 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'ZoneIId_Biomass', ...
    {'NoduleIId_BiomassRoot[c]', 'BacteroidIId_cpd11416[c]', 'BiomassZoneIId[c]'}, ...
    [-0.75 -0.25 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'ZoneIIp_Biomass', ...
    {'NoduleIIp_BiomassRoot[c]', 'BacteroidIIp_cpd11416[c]', 'BiomassZoneIIp[c]'}, ...
    [-0.75 -0.25 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'ZoneIZ_Biomass', ...
    {'NoduleIZ_BiomassRoot[c]', 'BacteroidIZ_cpd11416[c]', 'BiomassZoneIZ[c]'}, ...
    [-0.75 -0.25 1], 0, 0, 1000, 0);

% Make reaction to combine the nodule biomass
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'NoduleBiomass', ...
    {'BiomassZoneI[c]', 'BiomassZoneIId[c]', 'BiomassZoneIIp[c]', 'BiomassZoneIZ[c]', 'BiomassNodule[c]'}, ...
    [-0.05 -0.45 -0.45 -0.05 1], 0, 0, 1000, 0);

% Make reaction to combine plant biomass
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'PlantBiomass', ...
    {'BiomassRoot[c]', 'BiomassShoot[c]', 'BiomassPlant[c]'}, [-0.33333 -0.66667 1], 0, 0, 1000, 0);

% Make new overall Biomass reaction with nodule Biomass
nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'Biomass');
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'Biomass', ...
    {'BiomassPlant[c]', 'BiomassNodule[c]', 'BiomassTotal[c]'}, [-0.98 -0.02 1], 0, 0, 1000, 0);

% Add biomass exchange reaction
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_Biomass', ...
    {'BiomassTotal[c]'}, [-1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_NoduleBiomass', ...
    {'BiomassNodule[c]'}, [-1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'EX_PlantBiomass', ...
    {'BiomassPlant[c]'}, [-1], 0, 0, 1000, 0);

% Set the objective function
nodulatedPlantConnected.c(:) = 0;
nodulatedPlantConnected = changeObjective(nodulatedPlantConnected, 'Biomass');
optimizeCbModel(nodulatedPlantConnected)

%% Add maintenance costs

% Root and shoot maintenance based on GrowthRateComparison (Thomas Pfau)
MaintenanceRequirement = 0.6250 / 180 * 32 * 1000;
ATPMaintenanceShoot = find(ismember(nodulatedPlantConnected.rxns,'Leave_ATPSYN-RXN_C'));
ATPMaintenanceRoot = find(ismember(nodulatedPlantConnected.rxns,'Root_ATPSYN-RXN_C'));
nodulatedPlantConnected.lb(ATPMaintenanceShoot) = MaintenanceRequirement * 2/3;
nodulatedPlantConnected.lb(ATPMaintenanceRoot) = MaintenanceRequirement * 1/3;

% Add nodule maintenance costs
MaintenanceRequirementNodule = MaintenanceRequirement * 0.02;
ATPMaintenanceZI = findRxnIDs(nodulatedPlantConnected, 'NoduleI_ATPSYN-RXN_C');
ATPMaintenanceZIId = findRxnIDs(nodulatedPlantConnected, 'NoduleIId_ATPSYN-RXN_C');
ATPMaintenanceZIIp = findRxnIDs(nodulatedPlantConnected, 'NoduleIIp_ATPSYN-RXN_C');
ATPMaintenanceZIZ = findRxnIDs(nodulatedPlantConnected, 'NoduleIZ_ATPSYN-RXN_C');
ATPMaintenanceZIII = findRxnIDs(nodulatedPlantConnected, 'NoduleIII_ATPSYN-RXN_C');
nodulatedPlantConnected.lb(ATPMaintenanceZI) = MaintenanceRequirementNodule * 0.05;
nodulatedPlantConnected.lb(ATPMaintenanceZIId) = MaintenanceRequirementNodule * 0.45 * 0.75;
nodulatedPlantConnected.lb(ATPMaintenanceZIIp) = MaintenanceRequirementNodule * 0.45 * 0.75;
nodulatedPlantConnected.lb(ATPMaintenanceZIZ) = MaintenanceRequirementNodule * 0.05 * 0.75;
nodulatedPlantConnected.lb(ATPMaintenanceZIII) = MaintenanceRequirementNodule * 0.50 * 0.75;

% Add bacteroid maintenance costs
MaintenanceRequirementBacteroid = 8400 * 0.02 * 0.3;
ATPMaintenanceZIId = findRxnIDs(nodulatedPlantConnected, 'BacteroidIId_rxnNGAM');
ATPMaintenanceZIIp = findRxnIDs(nodulatedPlantConnected, 'BacteroidIIp_rxnNGAM');
ATPMaintenanceZIZ = findRxnIDs(nodulatedPlantConnected, 'BacteroidIZ_rxnNGAM');
ATPMaintenanceZIII = findRxnIDs(nodulatedPlantConnected, 'BacteroidIII_rxnNGAM');
nodulatedPlantConnected.lb(ATPMaintenanceZIId) = MaintenanceRequirementBacteroid * 0.45 * 0.25;
nodulatedPlantConnected.lb(ATPMaintenanceZIIp) = MaintenanceRequirementBacteroid * 0.45 * 0.25;
nodulatedPlantConnected.lb(ATPMaintenanceZIZ) = MaintenanceRequirementBacteroid * 0.05 * 0.25;
nodulatedPlantConnected.lb(ATPMaintenanceZIII) = MaintenanceRequirementBacteroid * 0.50 * 0.25;

% Switch name
nodulatedPlant = nodulatedPlantConnected;

%% Set boundaries to reasonable limit

% Change boundaries
for n = 1:length(nodulatedPlant.rxns)
    if nodulatedPlant.ub(n) >= 1000
        nodulatedPlant.ub(n) = 1000;
    end
    if nodulatedPlant.lb(n) <= -1000
        nodulatedPlant.lb(n) = -1000;
    end
end

%% Fix the rules and grRules
% Just in case

% Change model name
model = nodulatedPlant;

% Fix the grRules string formatting
for n = 1:length(model.grRules)
    if isstring(model.grRules{n})
        model.grRules{n} = str2mat(model.grRules{n});
    else
        model.grRules{n} = model.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(model.rules)
    if isstring(model.rules{n})
        model.rules{n} = str2mat(model.rules{n});
    else
        model.rules{n} = model.rules{n};
    end
end

% Fix the grRules
for n = 1:length(model.grRules)
    if ~isempty(model.grRules{n})
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(strrep(strrep(model.grRules{n}, ' or ', 'AAA'), ' or', ''),'AAA', ' or ');
        model.grRules{n} = strrep(strrep(strrep(model.grRules{n}, ' or ', 'AAA'), 'or ', ''),'AAA', ' or ');
    end
end

% Fix the rules
for n = 1:length(model.rules)
    if ~isempty(model.rules{n})
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(strrep(strrep(model.rules{n}, ' | ', 'AAA'), '| ', ''),'AAA', ' | ');
    end
end

% Fix the grRules string formatting
for n = 1:length(model.grRules)
    if isstring(model.grRules{n})
        model.grRules{n} = str2mat(model.grRules{n});
    else
        model.grRules{n} = model.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(model.rules)
    if isstring(model.rules{n})
        model.rules{n} = str2mat(model.rules{n});
    else
        model.rules{n} = model.rules{n};
    end
end

% Change model name
nodulatedPlant = model;

%% Remove deadends

model = nodulatedPlant;
deadEndMetabolites = detectDeadEnds(model);
if ~isempty(deadEndMetabolites)
    test = 1;
    while test > 0
        deadEndMetabolites = detectDeadEnds(model);
        if isempty(deadEndMetabolites)
            test = 0;
        else
            deadEndMetabolites3 = cell(length(deadEndMetabolites),1);
            for n = 1:length(deadEndMetabolites);
                deadEndMetabolites2 = num2cell(deadEndMetabolites);
                deadEndMetabolites3{n,1} = ...
                    model.mets{deadEndMetabolites2{n,1},1};
            end
            [deadEndReactions] = ...
                findRxnsFromMets(model,deadEndMetabolites3);
            model = tncore_remove_reactions(model,deadEndReactions);
            clear deadEndMetabolites
            clear deadEndMetabolites2
        end
    end
    model = tncore_remove(model);
end
nodulatedPlantOriginal = nodulatedPlant;
nodulatedPlant = model;

%% Set up plant transport boundaries
% This code is mostly taken from GrowthRateComparison of Thomas Pfau

% Turn off nitrogen uptake by the roots
NitrateUptakeRoot = find(ismember(nodulatedPlant.rxns,'Root_TEC_NITRATE'));
AmmoniumUptakeRoot = find(ismember(nodulatedPlant.rxns,'Root_TEC_AMMONIUM'));
nodulatedPlant.ub(NitrateUptakeRoot) = 0;
nodulatedPlant.ub(AmmoniumUptakeRoot) = 0;

% Turn off starch uptake
StarchImport = find(ismember(nodulatedPlant.rxns,'Leave_TEH_Starch'));
nodulatedPlant.ub(StarchImport) = 0;

% Make sure light is turned on
LightImport = find(ismember(nodulatedPlant.rxns,'Leave_TEH_Light'));
nodulatedPlant.ub(LightImport) = 1000;

% Make sure starch is accumulated
wrongBiomass = find(ismember(nodulatedPlant.rxns,'Leave_BiomassShootWithOutStarch'));
nodulatedPlant.ub(wrongBiomass) = 0;
nodulatedPlant.lb(wrongBiomass) = 0;

%% Save and clean workspace

save('allWorkspace.mat');
save('nodulatedPlant.mat', 'nodulatedPlant');
clear;
