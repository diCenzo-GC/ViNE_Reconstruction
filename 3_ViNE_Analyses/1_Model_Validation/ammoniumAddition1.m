%% Test effect of adding ammonium to soil
% Assume rate of N2-fixation changes but nodulation remains constant

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Determine ammonium uptake at maximal growth
model = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
sol = optimizeCbModel(model);
maxNH4 = sol.x(findRxnIDs(model, 'Root_TEC_AMMONIUM'));
model = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 0, 'u');

% Find position of ammonium uptake reaction
pos = findRxnIDs(model, 'Root_TEC_AMMONIUM');

% Get the original solution of the model
solOrig = optimizeCbModel(model);

% Prepare output variable
output_ammoniumNfixation_1 = cell(151, 5);

% Add ammonium
for n = 0:150
    modelTemp = model;
    % Set ammonium uptake maximum
    modelTemp.lb(pos) = 0;
    modelTemp.ub(pos) = n * (maxNH4 / 100);
    % Optimize growth
    sol2 = optimizeCbModel(modelTemp);
    % Store the output
    output_ammoniumNfixation_1{n+1,1} = modelTemp.ub(findRxnIDs(model, 'Root_TEC_AMMONIUM'));
    output_ammoniumNfixation_1{n+1,2} = sol2.x(findRxnIDs(model, 'PlantBiomass'));
    output_ammoniumNfixation_1{n+1,3} = sol2.x(findRxnIDs(model, 'NoduleBiomass'));
    output_ammoniumNfixation_1{n+1,4} = sol2.x(findRxnIDs(model, 'BacteroidIII_rxnConvFixed'));
    output_ammoniumNfixation_1{n+1,5} = sol2.x(findRxnIDs(model, 'Root_TEC_AMMONIUM'));
end

% Add headers
headers = {'Max_NH4_uptake', 'PlantBiomass', 'NoduleBiomass', 'N2_fix', 'NH4_uptake'};
output_ammoniumNfixation_1 = vertcat(headers, output_ammoniumNfixation_1);

% Clear unwanted variables
clearvars -except output_*

