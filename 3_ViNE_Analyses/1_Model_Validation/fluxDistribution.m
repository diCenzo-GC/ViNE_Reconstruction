%% Flux distribution patterns

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run FBA and FVA
sol = optimizeCbModel(model);
[minFlux, maxFlux] = fluxVariability(model, 99.9);

%% Get reaction lists

% Get nodule and bacteroid reactions
plantRxns = model.rxns(strmatch('Nodule', model.rxns));
plantRxns = setdiff(plantRxns, {'NoduleBiomass'});
bacteroidRxns = model.rxns(strmatch('Bacteroid', model.rxns));

% Get plant reactions
shootRxns = model.rxns(strmatch('Leave', model.rxns));
rootRxns = model.rxns(strmatch('Root', model.rxns));
plantRxns = vertcat(shootRxns, rootRxns);

% Get list of zone-independent nodule reactions
plantRxnsOverall = cell(length(plantRxns), 1);
for n = 1:length(plantRxns)
    temp = strsplit(plantRxns{n}, '_');
    plantRxnsOverall{n} = temp{1,2};
end
plantRxnsOverall = unique(plantRxnsOverall);

% Get list of zone-independent bacteroid reactions
bacteroidRxnsOverall = cell(length(bacteroidRxns), 1);
for n = 1:length(bacteroidRxns)
    temp = strsplit(bacteroidRxns{n}, '_');
    bacteroidRxnsOverall{n} = temp{1,2};
end
bacteroidRxnsOverall = unique(bacteroidRxnsOverall);

%% Parse and store the data

% Store the plant data as a matrix
temp = cell(length(plantRxnsOverall), 15);
for n = 1:length(plantRxns)
    rxnOrig = strsplit(plantRxns{n}, '_');
    pos = strmatch(rxnOrig{1,2}, plantRxnsOverall, 'exact');
    if strmatch('NoduleI', rxnOrig{1,1}, 'exact')
        temp{pos,1} = sol.x(n);
        temp{pos,6} = minFlux(n);
        temp{pos,11} = maxFlux(n);
    elseif strmatch('NoduleIId', rxnOrig{1,1}, 'exact')
        temp{pos,2} = sol.x(n);
        temp{pos,7} = minFlux(n);
        temp{pos,12} = maxFlux(n);
    elseif strmatch('NoduleIIp', rxnOrig{1,1}, 'exact')
        temp{pos,3} = sol.x(n);
        temp{pos,8} = minFlux(n);
        temp{pos,13} = maxFlux(n);
    elseif strmatch('NoduleIZ', rxnOrig{1,1}, 'exact')
        temp{pos,4} = sol.x(n);
        temp{pos,9} = minFlux(n);
        temp{pos,14} = maxFlux(n);
    elseif strmatch('NoduleIII', rxnOrig{1,1}, 'exact')
        temp{pos,5} = sol.x(n);
        temp{pos,10} = minFlux(n);
        temp{pos,15} = maxFlux(n);
    end
end
temp = horzcat(plantRxnsOverall, temp);
output_fluxDistribution_plant = ...
    vertcat({'Reaction','Zone_I_fba','Zone_IId_fba','Zone_IIp_fba',...
    'Zone_IZ_fba','Zone_III_fba','Zone_I_min','Zone_IId_min','Zone_IIp_min',...
    'Zone_IZ_min','Zone_III_min','Zone_I_max','Zone_IId_max','Zone_IIp_max',...
    'Zone_IZ_max','Zone_III_max'}, temp);
isEmpty = cellfun(@isempty, output_fluxDistribution_plant);
output_fluxDistribution_plant(isEmpty) = {0};

% Store the bacteroid data as a matrix
temp = cell(length(bacteroidRxnsOverall), 12);
for n = 1:length(bacteroidRxns)
    rxnOrig = strsplit(bacteroidRxns{n}, '_');
    pos = strmatch(rxnOrig{1,2}, bacteroidRxnsOverall, 'exact');
    rxnID = findRxnIDs(model, bacteroidRxns{n});
    if strmatch('BacteroidIId', rxnOrig{1,1}, 'exact')
        temp{pos,1} = sol.x(rxnID);
        temp{pos,5} = minFlux(rxnID);
        temp{pos,9} = maxFlux(rxnID);
    elseif strmatch('BacteroidIIp', rxnOrig{1,1}, 'exact')
        temp{pos,2} = sol.x(rxnID);
        temp{pos,6} = minFlux(rxnID);
        temp{pos,10} = maxFlux(rxnID);
    elseif strmatch('BacteroidIZ', rxnOrig{1,1}, 'exact')
        temp{pos,3} = sol.x(rxnID);
        temp{pos,7} = minFlux(rxnID);
        temp{pos,11} = maxFlux(rxnID);
    elseif strmatch('BacteroidIII', rxnOrig{1,1}, 'exact')
        temp{pos,4} = sol.x(rxnID);
        temp{pos,8} = minFlux(rxnID);
        temp{pos,12} = maxFlux(rxnID);
    end
end
temp = horzcat(bacteroidRxnsOverall, temp);
output_fluxDistribution_bacteroid = ...
    vertcat({'Reaction','Zone_IId_fba','Zone_IIp_fba',...
    'Zone_IZ_fba','Zone_III_fba','Zone_IId_min','Zone_IIp_min',...
    'Zone_IZ_min','Zone_III_min','Zone_IId_max','Zone_IIp_max',...
    'Zone_IZ_max','Zone_III_max'}, temp);
isEmpty = cellfun(@isempty, output_fluxDistribution_bacteroid);
output_fluxDistribution_bacteroid(isEmpty) = {0};

%% Clean the workspace

% Clear unwanted variables
clearvars -except output_*
