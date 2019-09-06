%% Zone-independent reaction deletion

% Get plant and bacteroid reactions
bacteroidRxns = model.rxns(strmatch('Bacteroid', model.rxns));
plantRxns = vertcat(model.rxns(strmatch('Nodule', model.rxns)), ...
    model.rxns(strmatch('Leave_', model.rxns)), model.rxns(strmatch('Root_', model.rxns)));
plantRxns = setdiff(plantRxns, {'NoduleBiomass'});

% Get list of zone-independent plant reactions
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

% Get plant reactions pairings
plantRxnPairs = cell(length(plantRxns), 2);
for n = 1:length(plantRxns)
    temp = strsplit(plantRxns{n}, '_');
    plantRxnPairs{n,1} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            plantRxnPairs{n,1} = strcat(plantRxnPairs{n,1}, '_', temp{m});
        end
    end
    plantRxnPairs{n,2} = plantRxns{n};
end

% Get bacteroid reactions pairings
bacteroidRxnPairs = cell(length(bacteroidRxns), 2);
for n = 1:length(bacteroidRxns)
    temp = strsplit(bacteroidRxns{n}, '_');
    bacteroidRxnPairs{n,1} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            bacteroidRxnPairs{n,1} = strcat(bacteroidRxnPairs{n,1}, '_', temp{m});
        end
    end
    bacteroidRxnPairs{n,2} = bacteroidRxns{n};
end

% Get the original solution during symbiosis
solOrig = optimizeCbModel(model);

% Nodule deletion analysis during symbiosis
output_zoneIndependentRxnDeletion_nodule = cell(length(plantRxnsOverall), 3);
for n = 1:length(plantRxnsOverall)
    pos = strmatch(plantRxnsOverall{n}, plantRxnPairs(:,1), 'exact');
    rxnGroup = plantRxnPairs(pos,2);
    modelTemp = tncore_remove_reactions(model, rxnGroup);
    sol = optimizeCbModel(modelTemp);
    output_zoneIndependentRxnDeletion_nodule{n,1} = plantRxnsOverall{n};
    output_zoneIndependentRxnDeletion_nodule{n,2} = round(sol.f / solOrig.f, 2);
end

% Bacteroid deletion analysis during symbiosis
output_zoneIndependentRxnDeletion_bacteriod = cell(length(bacteroidRxnsOverall), 2);
for n = 1:length(bacteroidRxnsOverall)
    pos = strmatch(bacteroidRxnsOverall{n}, bacteroidRxnPairs(:,1), 'exact');
    rxnGroup = bacteroidRxnPairs(pos,2);
    modelTemp = tncore_remove_reactions(model, rxnGroup);
    sol = optimizeCbModel(modelTemp);
    output_zoneIndependentRxnDeletion_bacteriod{n,1} = bacteroidRxnsOverall{n};
    output_zoneIndependentRxnDeletion_bacteriod{n,2} = round(sol.f / solOrig.f, 2);
end

% Get the original solution with ammonium
model2 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model2 = changeObjective(model2, 'PlantBiomass');
model2 = addExchangeRxn(model2, {'BiomassPlant[c]'}, 0, 1000000);
noduleReactions = model.rxns(strmatch('Nodule', model.rxns));
bacteroidReactions = model.rxns(strmatch('Bacteroid', model.rxns));
model2 = changeRxnBounds(model2, noduleReactions, 0, 'b');
model2 = changeRxnBounds(model2, bacteroidReactions, 0, 'b');
solOrig_amm = optimizeCbModel(model2);

% Nodule deletion analysis with ammonium
for n = 1:length(plantRxnsOverall)
    pos = strmatch(plantRxnsOverall{n}, plantRxnPairs(:,1), 'exact');
    rxnGroup = plantRxnPairs(pos,2);
    modelTemp = tncore_remove_reactions(model2, rxnGroup);
    sol = optimizeCbModel(modelTemp);
    output_zoneIndependentRxnDeletion_nodule{n,3} = round(sol.f / solOrig_amm.f, 2);
end