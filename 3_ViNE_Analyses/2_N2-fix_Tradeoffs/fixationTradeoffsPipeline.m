%% Initial fitness cost calculations

initialFitnessTests;

%% Calculate effects of rate of N2-fixation

% Load and prepare the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Pareto frontier for nitrogen fixation and growth
output_paretoNitrogenFixation = paretoFrontier(model, 'BacteroidIII_rxnConvFixed', 'Biomass', 100);

% Pareto frontier for nitrogen fixation and growth without oxygen limit
model2 = changeRxnBounds(model, 'TNIII_OXYGEN-MOLECULE', 1000000, 'u');
model2 = changeRxnBounds(model2, 'NoduleIII_EXCT_for_MNXM4_e0', 1000000, 'u');
output_paretoNitrogenFixationOxygen = paretoFrontier(model2, 'BacteroidIII_rxnConvFixed', 'Biomass', 100);

% Clean the workspace
clearvars -except output*

%% Calculate effects of rate of nodulation

% Load and prepare the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / 0.02;

% Prepare outout
output_paretoNodulation = cell(101,4);

% Run the analysis
for n = 0:200
    n
    nodule = n / 2000;
    plant = 1 - nodule;
    modelTemp = model;
    modelTemp.S(findMetIDs(modelTemp, 'BiomassNodule[c]'),findRxnIDs(modelTemp, 'Biomass')) = -1 * nodule;
    modelTemp.S(findMetIDs(modelTemp, 'BiomassPlant[c]'),findRxnIDs(modelTemp, 'Biomass')) = -1 * plant;
    modelTemp.lb(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) = n2fix * nodule;
    modelTemp.ub(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) = n2fix * nodule;
    modelTemp.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * nodule / 0.02;
    modelTemp.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * nodule / 0.02;
    modelTemp.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * nodule / 0.02;
    sol = optimizeCbModel(modelTemp);
    output_paretoNodulation{n+1,1} = nodule;
    if sol.f == 0
        output_paretoNodulation{n+1,2} = 0;
        output_paretoNodulation{n+1,3} = 0;
        output_paretoNodulation{n+1,4} = 0;
    else
        output_paretoNodulation{n+1,2} = sol.x(findRxnIDs(model, 'PlantBiomass'));
        output_paretoNodulation{n+1,3} = sol.x(findRxnIDs(model, 'NoduleBiomass'));
        output_paretoNodulation{n+1,4} = sol.x(findRxnIDs(model, 'Biomass'));
    end
end

% Add headers
headers = {'Percent_nodule', 'PlantBiomass', 'NoduleBiomass', 'TotalBiomass'};
output_paretoNodulation = vertcat(headers, output_paretoNodulation);

% Clean the workspace
clearvars -except output*

%% Examine effect of changing N2-fixation efficiency allowing nodulation to also vary

fixationEfficiencyAnalysis;

%% Clean workspace and save

clearvars -except output*
save('fixationTradeoffsOutput.mat');
clear;
