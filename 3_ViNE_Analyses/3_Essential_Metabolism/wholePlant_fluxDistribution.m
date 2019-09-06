%% Flux distribution patterns

% Run FBA during symbiosis
sol = optimizeCbModel(model);

% Run FVA during symbiosis
[minFlux, maxFlux] = fluxVariability(model);

% Set up non-symbiotic model
model2 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model2 = changeObjective(model2, 'PlantBiomass');
model2 = addExchangeRxn(model2, {'BiomassPlant[c]'}, 0, 1000000);
noduleReactions = model.rxns(strmatch('Nodule', model.rxns));
bacteroidReactions = model.rxns(strmatch('Bacteroid', model.rxns));
model2 = changeRxnBounds(model2, noduleReactions, 0, 'b');
model2 = changeRxnBounds(model2, bacteroidReactions, 0, 'b');

% Run FBA during ammonium uptake
sol_amm = optimizeCbModel(model2);

% Run FVA during ammonium uptake
[minFlux_amm, maxFlux_amm] = fluxVariability(model2);

% Get nodule and bacteroid reactions
noduleRxns = model.rxns(strmatch('Nodule', model.rxns));
noduleRxns = setdiff(noduleRxns, {'NoduleBiomass'});
bacteroidRxns = model.rxns(strmatch('Bacteroid', model.rxns));

% Get plant reactions
shootRxns = model.rxns(strmatch('Leave', model.rxns));
rootRxns = model.rxns(strmatch('Root', model.rxns));
plantRxns = vertcat(shootRxns, rootRxns);

% Get list of zone-independent nodule reactions
noduleRxnsOverall = cell(length(noduleRxns), 1);
for n = 1:length(noduleRxns)
    temp = strsplit(noduleRxns{n}, '_');
    noduleRxnsOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            noduleRxnsOverall{n} = strcat(noduleRxnsOverall{n}, '_', temp{m});
        end
    end
end
noduleRxnsOverall = unique(noduleRxnsOverall);

% Get list of zone-independent bacteroid reactions
bacteroidRxnsOverall = cell(length(bacteroidRxns), 1);
for n = 1:length(bacteroidRxns)
    temp = strsplit(bacteroidRxns{n}, '_');
    bacteroidRxnsOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            bacteroidRxnsOverall{n} = strcat(bacteroidRxnsOverall{n}, '_', temp{m});
        end
    end
end
bacteroidRxnsOverall = unique(bacteroidRxnsOverall);

% Get list of tissue-independent plant genes
plantRxnsOverall = cell(length(plantRxns), 1);
for n = 1:length(plantRxns)
    temp = strsplit(plantRxns{n}, '_');
    plantRxnsOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            plantRxnsOverall{n} = strcat(plantRxnsOverall{n}, '_', temp{m});
        end
    end
end
plantRxnsOverall = unique(plantRxnsOverall);

% Store the nodule data as a matrix
temp = cell(length(noduleRxnsOverall), 15);
for n = 1:length(noduleRxns)
    rxnOrig = strsplit(noduleRxns{n}, '_');
    temporary = rxnOrig{1,2};
    if length(rxnOrig) > 2
        for m = 3:length(rxnOrig)
            temporary = strcat(temporary, '_', rxnOrig{m});
        end
    end
    pos = strmatch(temporary, noduleRxnsOverall, 'exact');
    if strmatch('NoduleI', rxnOrig{1,1}, 'exact')
        temp{pos,1} = sol.x(findRxnIDs(model, noduleRxns{n}));
        temp{pos,6} = minFlux(findRxnIDs(model, noduleRxns{n}));
        temp{pos,11} = maxFlux(findRxnIDs(model, noduleRxns{n}));
    elseif strmatch('NoduleIId', rxnOrig{1,1}, 'exact')
        temp{pos,2} = sol.x(findRxnIDs(model, noduleRxns{n}));
        temp{pos,7} = minFlux(findRxnIDs(model, noduleRxns{n}));
        temp{pos,12} = maxFlux(findRxnIDs(model, noduleRxns{n}));
    elseif strmatch('NoduleIIp', rxnOrig{1,1}, 'exact')
        temp{pos,3} = sol.x(findRxnIDs(model, noduleRxns{n}));
        temp{pos,8} = minFlux(findRxnIDs(model, noduleRxns{n}));
        temp{pos,13} = maxFlux(findRxnIDs(model, noduleRxns{n}));
    elseif strmatch('NoduleIZ', rxnOrig{1,1}, 'exact')
        temp{pos,4} = sol.x(findRxnIDs(model, noduleRxns{n}));
        temp{pos,9} = minFlux(findRxnIDs(model, noduleRxns{n}));
        temp{pos,14} = maxFlux(findRxnIDs(model, noduleRxns{n}));
    elseif strmatch('NoduleIII', rxnOrig{1,1}, 'exact')
        temp{pos,5} = sol.x(findRxnIDs(model, noduleRxns{n}));
        temp{pos,10} = minFlux(findRxnIDs(model, noduleRxns{n}));
        temp{pos,15} = maxFlux(findRxnIDs(model, noduleRxns{n}));
    end
end
temp = horzcat(noduleRxnsOverall, temp);
output_fluxDistribution_nodule = ...
    vertcat({'Reaction','Zone_I_fba','Zone_IId_fba','Zone_IIp_fba',...
    'Zone_IZ_fba','Zone_III_fba','Zone_I_min','Zone_IId_min','Zone_IIp_min',...
    'Zone_IZ_min','Zone_III_min','Zone_I_max','Zone_IId_max','Zone_IIp_max',...
    'Zone_IZ_max','Zone_III_max'}, temp);
isEmpty = cellfun(@isempty, output_fluxDistribution_nodule);
output_fluxDistribution_nodule(isEmpty) = {0};

% Store the bacteroid data as a matrix
temp = cell(length(bacteroidRxnsOverall), 12);
for n = 1:length(bacteroidRxns)
    rxnOrig = strsplit(bacteroidRxns{n}, '_');
    temporary = rxnOrig{1,2};
    if length(rxnOrig) > 2
        for m = 3:length(rxnOrig)
            temporary = strcat(temporary, '_', rxnOrig{m});
        end
    end
    pos = strmatch(temporary, bacteroidRxnsOverall, 'exact');
    if strmatch('BacteroidIId', rxnOrig{1,1}, 'exact')
        temp{pos,1} = sol.x(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,5} = minFlux(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,9} = maxFlux(findRxnIDs(model, bacteroidRxns{n}));
    elseif strmatch('BacteroidIIp', rxnOrig{1,1}, 'exact')
        temp{pos,2} = sol.x(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,6} = minFlux(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,10} = maxFlux(findRxnIDs(model, bacteroidRxns{n}));
    elseif strmatch('BacteroidIZ', rxnOrig{1,1}, 'exact')
        temp{pos,3} = sol.x(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,7} = minFlux(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,11} = maxFlux(findRxnIDs(model, bacteroidRxns{n}));
    elseif strmatch('BacteroidIII', rxnOrig{1,1}, 'exact')
        temp{pos,4} = sol.x(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,8} = minFlux(findRxnIDs(model, bacteroidRxns{n}));
        temp{pos,12} = maxFlux(findRxnIDs(model, bacteroidRxns{n}));
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

% Store the plant data as a matrix
temp = cell(length(plantRxnsOverall), 12);
for n = 1:length(plantRxns)
    rxnOrig = strsplit(plantRxns{n}, '_');
    temporary = rxnOrig{1,2};
    if length(rxnOrig) > 2
        for m = 3:length(rxnOrig)
            temporary = strcat(temporary, '_', rxnOrig{m});
        end
    end
    pos = strmatch(temporary, plantRxnsOverall, 'exact');
    if strmatch('Leave', rxnOrig{1,1}, 'exact')
        temp{pos,1} = sol.x(findRxnIDs(model, plantRxns{n}));
        temp{pos,3} = sol_amm.x(findRxnIDs(model, plantRxns{n}));
        temp{pos,5} = minFlux(findRxnIDs(model, plantRxns{n}));
        temp{pos,7} = minFlux_amm(findRxnIDs(model, plantRxns{n}));
        temp{pos,9} = maxFlux(findRxnIDs(model, plantRxns{n}));
        temp{pos,11} = maxFlux_amm(findRxnIDs(model, plantRxns{n}));
    elseif strmatch('Root', rxnOrig{1,1}, 'exact')
        temp{pos,2} = sol.x(findRxnIDs(model, plantRxns{n}));
        temp{pos,4} = sol_amm.x(findRxnIDs(model, plantRxns{n}));
        temp{pos,6} = minFlux(findRxnIDs(model, plantRxns{n}));
        temp{pos,8} = minFlux_amm(findRxnIDs(model, plantRxns{n}));
        temp{pos,10} = maxFlux(findRxnIDs(model, plantRxns{n}));
        temp{pos,12} = maxFlux_amm(findRxnIDs(model, plantRxns{n}));
    end
end
temp = horzcat(plantRxnsOverall, temp);
output_fluxDistribution_plant = ...
    vertcat({'Reaction', 'Shoot_nod_fba', 'Root_nod_fba', 'Shoot_amm_fba', ...
    'Root_amm_fba', 'Shoot_nod_min', 'Root_nod_min', 'Shoot_amm_min', ...
    'Root_amm_min', 'Shoot_nod_max', 'Root_nod_max', 'Shoot_amm_max', 'Root_amm_max'}, temp);
isEmpty = cellfun(@isempty, output_fluxDistribution_plant);
output_fluxDistribution_plant(isEmpty) = {0};
