%% Test effect of changing N2-fixation efficiency - limited O2, no nodule max

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = 0.01 * sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of nitrogen fixation reaction
pos = findRxnIDs(model, 'BacteroidIII_rxnConvFixed');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Determine the standard rate of N2-fixation
n2fixMax = solOrig.x(pos);

% Prepare output variable
output_fixationEfficiency_max0_limO2 = cell(100, 5);

% Vary the rate of N2-fixation
for n = 1:100
    modelTemp = model;
    % Set efficiency
    n2fixNew = n2fix * n / 100;
    nodule = 0.01 * n2fixMax / n2fixNew;
    % Set initial nodule amount to be 10% max
    if nodule > 0.1
        nodule = 0.1;
    end
    % Get original solution
    sol = optimizeCbModel(modelTemp);
    x = 1;
    % Get the initial local maximum
    while x == 1
        % Record old solution
        solOld = sol;
        % Adjust oxygen and nitrogen fixation
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Adjust the nodule accordingly
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
        nodule = 0.01 * sol.x(pos) / n2fixNew;
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        sol = optimizeCbModel(modelTemp);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Stop the loop when local maximum obtained
        if sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
            x = 0;
            nodule = 0.01 * sol.x(pos) / n2fixNew;
        else
            nodule = 0.01 * sol.x(pos) / n2fixNew;
        end
    end
    y = 1;
    % See if increasing nodule amount can increase growth rate, iteratively
    while y == 1
        % Record old values
        noduleOld = nodule;
        solOld2 = sol;
        % Increase the amount of nodule by 1%
        nodule = nodule + 0.001;
        x = 1;
        % Get the local maximum as done abolve
        while x == 1
            solOld = sol;
            modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
            nodule = 0.01 * sol.x(pos) / n2fixNew;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            if sol.f == 0
                x = 0;
            elseif sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
                x = 0;
                nodule = 0.01 * sol.x(pos) / n2fixNew;
            else
                nodule = 0.01 * sol.x(pos) / n2fixNew;
            end
        end
        % If the local maximum is lower than the previous local maximum,
        % then the global maximum has already been found and we can stop
        if sol.f <= solOld2.f
            y = 0;
        end
    end
    % Save the output
    output_fixationEfficiency_max0_limO2{n,1} = n2fixNew;
    output_fixationEfficiency_max0_limO2{n,2} = solOld2.x(findRxnIDs(model, 'PlantBiomass'));
    output_fixationEfficiency_max0_limO2{n,3} = solOld2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_fixationEfficiency_max0_limO2{n,4} = output_fixationEfficiency_max0_limO2{n,3} / ...
        (output_fixationEfficiency_max0_limO2{n,2} + output_fixationEfficiency_max0_limO2{n,3});
    output_fixationEfficiency_max0_limO2{n,5} = solOld2.x(pos);
end

% Add headers
headers = {'Fixation_efficiency', 'Plant_biomass', 'Nodule_biomass', ...
    'Percent_nodule', 'Fixation_total'};
output_fixationEfficiency_max0_limO2 = vertcat(headers, output_fixationEfficiency_max0_limO2);

% Clean workspace
clearvars -except output*

%% Test effect of changing N2-fixation efficiency - unlimited O2, no nodule max

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Remove oxygen limit
model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) = 1000000;

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = 0.01 * sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of nitrogen fixation reaction
pos = findRxnIDs(model, 'BacteroidIII_rxnConvFixed');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Determine the standard rate of N2-fixation
n2fixMax = solOrig.x(pos);

% Prepare output variable
output_fixationEfficiency_max0_unlimO2 = cell(100, 5);

% Vary the rate of N2-fixation
for n = 1:100
    modelTemp = model;
    % Set efficiency
    n2fixNew = n2fix * n / 100;
    nodule = 0.01 * n2fixMax / n2fixNew;
    % Set initial nodule amount to be 10% max
    if nodule > 0.1
        nodule = 0.1;
    end
    % Get original solution
    sol = optimizeCbModel(modelTemp);
    x = 1;
    % Get the initial local maximum
    while x == 1
        % Record old solution
        solOld = sol;
        % Adjust oxygen and nitrogen fixation
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Adjust the nodule accordingly
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
        nodule = 0.01 * sol.x(pos) / n2fixNew;
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        sol = optimizeCbModel(modelTemp);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Stop the loop when local maximum obtained
        if sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
            x = 0;
            nodule = 0.01 * sol.x(pos) / n2fixNew;
        else
            nodule = 0.01 * sol.x(pos) / n2fixNew;
        end
    end
    y = 1;
    % See if increasing nodule amount can increase growth rate, iteratively
    while y == 1
        % Record old values
        noduleOld = nodule;
        solOld2 = sol;
        % Increase the amount of nodule by 1%
        nodule = nodule + 0.001;
        x = 1;
        % Get the local maximum as done abolve
        while x == 1
            solOld = sol;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
            nodule = 0.01 * sol.x(pos) / n2fixNew;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            if sol.f == 0
                x = 0;
            elseif sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
                x = 0;
                nodule = 0.01 * sol.x(pos) / n2fixNew;
            else
                nodule = 0.01 * sol.x(pos) / n2fixNew;
            end
        end
        % If the local maximum is lower than the previous local maximum,
        % then the global maximum has already been found and we can stop
        if sol.f <= solOld2.f
            y = 0;
        end
    end
    % Save the output
    output_fixationEfficiency_max0_unlimO2{n,1} = n2fixNew;
    output_fixationEfficiency_max0_unlimO2{n,2} = solOld2.x(findRxnIDs(model, 'PlantBiomass'));
    output_fixationEfficiency_max0_unlimO2{n,3} = solOld2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_fixationEfficiency_max0_unlimO2{n,4} = output_fixationEfficiency_max0_unlimO2{n,3} / ...
        (output_fixationEfficiency_max0_unlimO2{n,2} + output_fixationEfficiency_max0_unlimO2{n,3});
    output_fixationEfficiency_max0_unlimO2{n,5} = solOld2.x(pos);
end

% Add headers
headers = {'Fixation_efficiency', 'Plant_biomass', 'Nodule_biomass', ...
    'Percent_nodule', 'Fixation_total'};
output_fixationEfficiency_max0_unlimO2 = vertcat(headers, output_fixationEfficiency_max0_unlimO2);

% Clean workspace
clearvars -except output*

%% Test effect of changing N2-fixation efficiency - limited O2, 10% nodule max

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = 0.01 * sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of nitrogen fixation reaction
pos = findRxnIDs(model, 'BacteroidIII_rxnConvFixed');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Determine the standard rate of N2-fixation
n2fixMax = solOrig.x(pos);

% Prepare output variable
output_fixationEfficiency_max10_limO2 = cell(100, 5);

% Vary the rate of N2-fixation
for n = 1:100
    modelTemp = model;
    % Set efficiency
    n2fixNew = n2fix * n / 100;
    nodule = 0.01 * n2fixMax / n2fixNew;
    % Set initial nodule amount to be 10% max
    if nodule > 0.1
        nodule = 0.1;
    end
    % Get original solution
    sol = optimizeCbModel(modelTemp);
    x = 1;
    % Get the initial local maximum
    while x == 1
        % Record old solution
        solOld = sol;
        % Adjust oxygen and nitrogen fixation
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Adjust the nodule accordingly
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
        nodule = 0.01 * sol.x(pos) / n2fixNew;
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        sol = optimizeCbModel(modelTemp);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Stop the loop when local maximum obtained
        if sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
            x = 0;
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
        else
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
        end
    end
    y = 1;
    % See if increasing nodule amount can increase growth rate, iteratively
    while y == 1
        % Record old values
        noduleOld = nodule;
        solOld2 = sol;
        % Increase the amount of nodule by 1%
        nodule = nodule + 0.001;
                if nodule > 0.1
                    nodule = 0.1;
                end
        x = 1;
        % Get the local maximum as done abolve
        while x == 1
            solOld = sol;
            modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
            nodule = 0.01 * sol.x(pos) / n2fixNew;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            if sol.f == 0
                x = 0;
            elseif sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
                x = 0;
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
            else
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
            end
        end
        % If the local maximum is lower than the previous local maximum,
        % then the global maximum has already been found and we can stop
        if sol.f <= solOld2.f
            y = 0;
        end
    end
    % Save the output
    output_fixationEfficiency_max10_limO2{n,1} = n2fixNew;
    output_fixationEfficiency_max10_limO2{n,2} = solOld2.x(findRxnIDs(model, 'PlantBiomass'));
    output_fixationEfficiency_max10_limO2{n,3} = solOld2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_fixationEfficiency_max10_limO2{n,4} = output_fixationEfficiency_max10_limO2{n,3} / ...
        (output_fixationEfficiency_max10_limO2{n,2} + output_fixationEfficiency_max10_limO2{n,3});
    output_fixationEfficiency_max10_limO2{n,5} = solOld2.x(pos);
end

% Add headers
headers = {'Fixation_efficiency', 'Plant_biomass', 'Nodule_biomass', ...
    'Percent_nodule', 'Fixation_total'};
output_fixationEfficiency_max10_limO2 = vertcat(headers, output_fixationEfficiency_max10_limO2);

% Clean workspace
clearvars -except output*

%% Test effect of changing N2-fixation efficiency - unlimited O2, 10% nodule max

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Remove oxygen limit
model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) = 1000000;

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = 0.01 * sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of nitrogen fixation reaction
pos = findRxnIDs(model, 'BacteroidIII_rxnConvFixed');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Determine the standard rate of N2-fixation
n2fixMax = solOrig.x(pos);

% Prepare output variable
output_fixationEfficiency_max10_unlimO2 = cell(100, 5);

% Vary the rate of N2-fixation
for n = 1:100
    modelTemp = model;
    % Set efficiency
    n2fixNew = n2fix * n / 100;
    nodule = 0.01 * n2fixMax / n2fixNew;
    % Set initial nodule amount to be 10% max
    if nodule > 0.1
        nodule = 0.1;
    end
    % Get original solution
    sol = optimizeCbModel(modelTemp);
    x = 1;
    % Get the initial local maximum
    while x == 1
        % Record old solution
        solOld = sol;
        % Adjust oxygen and nitrogen fixation
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Adjust the nodule accordingly
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
        nodule = 0.01 * sol.x(pos) / n2fixNew;
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        sol = optimizeCbModel(modelTemp);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Stop the loop when local maximum obtained
        if sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
            x = 0;
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
        else
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
        end
    end
    y = 1;
    % See if increasing nodule amount can increase growth rate, iteratively
    while y == 1
        % Record old values
        noduleOld = nodule;
        solOld2 = sol;
        % Increase the amount of nodule by 1%
        nodule = nodule + 0.001;
                if nodule > 0.1
                    nodule = 0.1;
                end
        x = 1;
        % Get the local maximum as done abolve
        while x == 1
            solOld = sol;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
            nodule = 0.01 * sol.x(pos) / n2fixNew;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            if sol.f == 0
                x = 0;
            elseif sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
                x = 0;
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
            else
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.1
                    nodule = 0.1;
                end
            end
        end
        % If the local maximum is lower than the previous local maximum,
        % then the global maximum has already been found and we can stop
        if sol.f <= solOld2.f
            y = 0;
        end
    end
    % Save the output
    output_fixationEfficiency_max10_unlimO2{n,1} = n2fixNew;
    output_fixationEfficiency_max10_unlimO2{n,2} = solOld2.x(findRxnIDs(model, 'PlantBiomass'));
    output_fixationEfficiency_max10_unlimO2{n,3} = solOld2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_fixationEfficiency_max10_unlimO2{n,4} = output_fixationEfficiency_max10_unlimO2{n,3} / ...
        (output_fixationEfficiency_max10_unlimO2{n,2} + output_fixationEfficiency_max10_unlimO2{n,3});
    output_fixationEfficiency_max10_unlimO2{n,5} = solOld2.x(pos);
end

% Add headers
headers = {'Fixation_efficiency', 'Plant_biomass', 'Nodule_biomass', ...
    'Percent_nodule', 'Fixation_total'};
output_fixationEfficiency_max10_unlimO2 = vertcat(headers, output_fixationEfficiency_max10_unlimO2);

% Clean workspace
clearvars -except output*

%% Test effect of changing N2-fixation efficiency - limited O2, 5% nodule max

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = 0.01 * sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of nitrogen fixation reaction
pos = findRxnIDs(model, 'BacteroidIII_rxnConvFixed');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Determine the standard rate of N2-fixation
n2fixMax = solOrig.x(pos);

% Prepare output variable
output_fixationEfficiency_max5_limO2 = cell(100, 5);

% Vary the rate of N2-fixation
for n = 1:100
    modelTemp = model;
    % Set efficiency
    n2fixNew = n2fix * n / 100;
    nodule = 0.01 * n2fixMax / n2fixNew;
    % Set initial nodule amount to be 10% max
    if nodule > 0.1
        nodule = 0.1;
    end
    % Get original solution
    sol = optimizeCbModel(modelTemp);
    x = 1;
    % Get the initial local maximum
    while x == 1
        % Record old solution
        solOld = sol;
        % Adjust oxygen and nitrogen fixation
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Adjust the nodule accordingly
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
        nodule = 0.01 * sol.x(pos) / n2fixNew;
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        sol = optimizeCbModel(modelTemp);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Stop the loop when local maximum obtained
        if sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
            x = 0;
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
        else
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
        end
    end
    y = 1;
    % See if increasing nodule amount can increase growth rate, iteratively
    while y == 1
        % Record old values
        noduleOld = nodule;
        solOld2 = sol;
        % Increase the amount of nodule by 1%
        nodule = nodule + 0.001;
                if nodule > 0.05
                    nodule = 0.05;
                end
        x = 1;
        % Get the local maximum as done abolve
        while x == 1
            solOld = sol;
            modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * (nodule / 0.02);
            nodule = 0.01 * sol.x(pos) / n2fixNew;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            if sol.f == 0
                x = 0;
            elseif sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
                x = 0;
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
            else
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
            end
        end
        % If the local maximum is lower than the previous local maximum,
        % then the global maximum has already been found and we can stop
        if sol.f <= solOld2.f
            y = 0;
        end
    end
    % Save the output
    output_fixationEfficiency_max5_limO2{n,1} = n2fixNew;
    output_fixationEfficiency_max5_limO2{n,2} = solOld2.x(findRxnIDs(model, 'PlantBiomass'));
    output_fixationEfficiency_max5_limO2{n,3} = solOld2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_fixationEfficiency_max5_limO2{n,4} = output_fixationEfficiency_max5_limO2{n,3} / ...
        (output_fixationEfficiency_max5_limO2{n,2} + output_fixationEfficiency_max5_limO2{n,3});
    output_fixationEfficiency_max5_limO2{n,5} = solOld2.x(pos);
end

% Add headers
headers = {'Fixation_efficiency', 'Plant_biomass', 'Nodule_biomass', ...
    'Percent_nodule', 'Fixation_total'};
output_fixationEfficiency_max5_limO2 = vertcat(headers, output_fixationEfficiency_max5_limO2);

% Clean workspace
clearvars -except output*

%% Test effect of changing N2-fixation efficiency - unlimited O2, 5% nodule max

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Add ammonium sinks
model = addReaction(model, 'NoduleIII_SINK_MNXM15', {'NoduleIII_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);
model = addReaction(model, 'Root_SINK_MNXM15', {'Root_MNXM15[C]'}, [-1], 0, 0, 1000000, 0);

% Remove oxygen limit
model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) = 1000000;

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = 0.01 * sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of nitrogen fixation reaction
pos = findRxnIDs(model, 'BacteroidIII_rxnConvFixed');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Determine the standard rate of N2-fixation
n2fixMax = solOrig.x(pos);

% Prepare output variable
output_fixationEfficiency_max5_unlimO2 = cell(100, 5);

% Vary the rate of N2-fixation
for n = 1:100
    modelTemp = model;
    % Set efficiency
    n2fixNew = n2fix * n / 100;
    nodule = 0.01 * n2fixMax / n2fixNew;
    % Set initial nodule amount to be 10% max
    if nodule > 0.1
        nodule = 0.1;
    end
    % Get original solution
    sol = optimizeCbModel(modelTemp);
    x = 1;
    % Get the initial local maximum
    while x == 1
        % Record old solution
        solOld = sol;
        % Adjust oxygen and nitrogen fixation
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Adjust the nodule accordingly
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
        nodule = 0.01 * sol.x(pos) / n2fixNew;
        modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
        sol = optimizeCbModel(modelTemp);
        % Get new growth
        sol = optimizeCbModel(modelTemp);
        % Stop the loop when local maximum obtained
        if sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
            x = 0;
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
        else
            nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
        end
    end
    y = 1;
    % See if increasing nodule amount can increase growth rate, iteratively
    while y == 1
        % Record old values
        noduleOld = nodule;
        solOld2 = sol;
        % Increase the amount of nodule by 1%
        nodule = nodule + 0.001;
                if nodule > 0.05
                    nodule = 0.05;
                end
        x = 1;
        % Get the local maximum as done abolve
        while x == 1
            solOld = sol;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = (nodule * sol.x(findRxnIDs(model, 'PlantBiomass'))) / (1 - nodule);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * (nodule / 0.02);
            modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * (nodule / 0.02);
            nodule = 0.01 * sol.x(pos) / n2fixNew;
            modelTemp.ub(pos) = n2fixNew * (nodule / 0.01);
            sol = optimizeCbModel(modelTemp);
            if sol.f == 0
                x = 0;
            elseif sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100) > 0.995 * n2fixNew && sol.x(pos) / ((sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(modelTemp, 'PlantBiomass')))) * 100)  < 1.005 * n2fixNew
                x = 0;
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
            else
                nodule = 0.01 * sol.x(pos) / n2fixNew;
                if nodule > 0.05
                    nodule = 0.05;
                end
            end
        end
        % If the local maximum is lower than the previous local maximum,
        % then the global maximum has already been found and we can stop
        if sol.f <= solOld2.f
            y = 0;
        end
    end
    % Save the output
    output_fixationEfficiency_max5_unlimO2{n,1} = n2fixNew;
    output_fixationEfficiency_max5_unlimO2{n,2} = solOld2.x(findRxnIDs(model, 'PlantBiomass'));
    output_fixationEfficiency_max5_unlimO2{n,3} = solOld2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_fixationEfficiency_max5_unlimO2{n,4} = output_fixationEfficiency_max5_unlimO2{n,3} / ...
        (output_fixationEfficiency_max5_unlimO2{n,2} + output_fixationEfficiency_max5_unlimO2{n,3});
    output_fixationEfficiency_max5_unlimO2{n,5} = solOld2.x(pos);
end

% Add headers
headers = {'Fixation_efficiency', 'Plant_biomass', 'Nodule_biomass', ...
    'Percent_nodule', 'Fixation_total'};
output_fixationEfficiency_max5_unlimO2 = vertcat(headers, output_fixationEfficiency_max5_unlimO2);

% Clean workspace
clearvars -except output*
