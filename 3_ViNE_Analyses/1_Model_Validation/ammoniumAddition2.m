%% Test effect of adding ammonium to soil
% Assume rate of N2-fixation remains constant but nodulation varies

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Check standard flux through nodule biomass reaction
sol = optimizeCbModel(model);
noduleFlux = sol.x(findRxnIDs(model, 'NoduleBiomass'));

% Remove nodule from biomass
model.S(findMetIDs(model, 'BiomassNodule[c]'),findRxnIDs(model, 'Biomass')) = 0;

% Force nodule biomass through a sink
model = addExchangeRxn(model, {'BiomassNodule[c]'}, 0, 1000000);
model.lb(findRxnIDs(model, 'EX_BiomassNodule[c]')) = noduleFlux;

% Determine ammonium uptake at maximal growth
model = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
sol = optimizeCbModel(model);
maxNH4 = sol.x(findRxnIDs(model, 'Root_TEC_AMMONIUM'));
model = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 0, 'u');

% Determine rate of N2-fixation per nodule mass
sol = optimizeCbModel(model);
n2fix = sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / ...
    (sol.x(findRxnIDs(model, 'NoduleBiomass')) / (sol.x(findRxnIDs(model, 'PlantBiomass')) + sol.x(findRxnIDs(model, 'NoduleBiomass'))));

% Find position of ammonium uptake reaction
pos = findRxnIDs(model, 'Root_TEC_AMMONIUM');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Prepare output variable
output_ammoniumNfixation_2 = cell(151, 5);

% Add ammonium
for n = 0:150
    n
    % Set ammonium uptake maximum
    modelTemp = model;
    if n == 0
        modelTemp.lb(pos) = 0;
        modelTemp.ub(pos) = 0;
    else
        modelTemp.lb(pos) = 0;
        modelTemp.ub(pos) = n * (maxNH4 / 100);
    end
    % Solve with new ammonium uptake
    sol = optimizeCbModel(modelTemp);
    % Set the nodule - plant ratio to be proper
    modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = ...
        (sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) * sol.x(findRxnIDs(model, 'PlantBiomass'))) ...
        / (n2fix - sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')));
    % Get solution with new amount of nodule
    sol = optimizeCbModel(modelTemp);
    % Adjust all the NGAM reactions
    modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = ...
        (model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = ...
        (model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = ...
        (model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = ...
        (model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = ...
        (model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = ...
        (model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = ...
        (model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = ...
        (model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = ...
        (model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * ...
        sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
        (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
    % Iteratively adjust model to get correct oxygen and ATPM values
    x = 1;
    while x == 1
        % First solve
        sol = optimizeCbModel(modelTemp);
        % Set the nodule - plant ratio to be proper
        modelTemp.lb(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) = ...
            (sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) * sol.x(findRxnIDs(model, 'PlantBiomass'))) ...
            / (n2fix - sol.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')));
        % Get solution with new amount of nodule
        sol = optimizeCbModel(modelTemp);
        % Adjust all the NGAM reactions
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIId_rxnNGAM')) = ...
            (model.lb(findRxnIDs(model, 'BacteroidIId_rxnNGAM')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIIp_rxnNGAM')) = ...
            (model.lb(findRxnIDs(model, 'BacteroidIIp_rxnNGAM')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIZ_rxnNGAM')) = ...
            (model.lb(findRxnIDs(model, 'BacteroidIZ_rxnNGAM')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'BacteroidIII_rxnNGAM')) = ...
            (model.lb(findRxnIDs(model, 'BacteroidIII_rxnNGAM')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleI_MNXR96136_C')) = ...
            (model.lb(findRxnIDs(model, 'NoduleI_MNXR96136_C')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIId_MNXR96136_C')) = ...
            (model.lb(findRxnIDs(model, 'NoduleIIp_MNXR96136_C')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIIp_MNXR96136_C')) = ...
            (model.lb(findRxnIDs(model, 'NoduleIId_MNXR96136_C')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIZ_MNXR96136_C')) = ...
            (model.lb(findRxnIDs(model, 'NoduleIZ_MNXR96136_C')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        modelTemp.lb(findRxnIDs(modelTemp, 'NoduleIII_MNXR96136_C')) = ...
            (model.lb(findRxnIDs(model, 'NoduleIII_MNXR96136_C')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        % Get the new solution with all the NGAM reactions set
        sol = optimizeCbModel(modelTemp);
        % Adjust the maximum oxygen uptake rate
        modelTemp.ub(findRxnIDs(modelTemp, 'TNIII_OXYGEN-MOLECULE')) = ...
            (model.ub(findRxnIDs(model, 'TNIII_OXYGEN-MOLECULE')) * ...
            sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) / ...
            (sol.x(findRxnIDs(modelTemp, 'EX_BiomassNodule[c]')) + sol.x(findRxnIDs(model, 'PlantBiomass')))) / 0.02;
        % Solve with the proper oxygen uptake rate
        sol2 = optimizeCbModel(modelTemp);
        % Repeat if the N2-fiation rate is not within 1% of the desired value
        if (sol2.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / (sol2.x(findRxnIDs(model, 'NoduleBiomass')) / (sol2.x(findRxnIDs(model, 'NoduleBiomass')) + sol2.x(findRxnIDs(model, 'PlantBiomass'))))) ...
                / (solOrig.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / (solOrig.x(findRxnIDs(model, 'NoduleBiomass')) / (solOrig.x(findRxnIDs(model, 'NoduleBiomass')) + solOrig.x(findRxnIDs(model, 'PlantBiomass'))))) > 0.999 && ...
                (sol2.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / (sol2.x(findRxnIDs(model, 'NoduleBiomass')) / (sol2.x(findRxnIDs(model, 'NoduleBiomass')) + sol2.x(findRxnIDs(model, 'PlantBiomass'))))) ...
                / (solOrig.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')) / (solOrig.x(findRxnIDs(model, 'NoduleBiomass')) / (solOrig.x(findRxnIDs(model, 'NoduleBiomass')) + solOrig.x(findRxnIDs(model, 'PlantBiomass'))))) < 1.001
            x = 0;
        elseif round( sol2.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed')), 3) == 0
            x = 0;
        end
    end
    % Store the output
    output_ammoniumNfixation_2{n+1,1} = modelTemp.ub(findRxnIDs(model, 'Root_TEC_AMMONIUM'));
    output_ammoniumNfixation_2{n+1,2} = sol2.x(findRxnIDs(model, 'PlantBiomass'));
    output_ammoniumNfixation_2{n+1,3} = sol2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_ammoniumNfixation_2{n+1,4} = sol2.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed'));
    output_ammoniumNfixation_2{n+1,5} = sol2.x(findRxnIDs(model, 'Root_TEC_AMMONIUM'));
end

% Add headers
headers = {'Max_NH4_uptake', 'PlantBiomass', 'NoduleBiomass', 'N2_fix', 'NH4_uptake'};
output_ammoniumNfixation_2 = vertcat(headers, output_ammoniumNfixation_2);

% Clear unwanted variables
clearvars -except output_*

