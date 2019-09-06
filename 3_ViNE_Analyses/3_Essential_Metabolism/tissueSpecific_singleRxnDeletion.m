%% Tissue-specific plant reaction deletion

% Get plant reactions
shootRxns = model.rxns(strmatch('Leave', model.rxns));
rootRxns = model.rxns(strmatch('Root', model.rxns));
plantRxns = vertcat(shootRxns, rootRxns);

% Perform gene deletion analyses
grRatioPlant = singleRxnDeletion(model, 'FBA', plantRxns);
grRatioPlant(isnan(grRatioPlant)) = 0;

% Perform gene deletion analysis with ammonium
model2 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model2 = changeObjective(model2, 'PlantBiomass');
model2 = addExchangeRxn(model2, {'BiomassPlant[c]'}, 0, 1000000);
noduleReactions = model.rxns(strmatch('Nodule', model.rxns));
bacteroidReactions = model.rxns(strmatch('Bacteroid', model.rxns));
model2 = changeRxnBounds(model2, noduleReactions, 0, 'b');
model2 = changeRxnBounds(model2, bacteroidReactions, 0, 'b');
grRatioPlant_amm = singleRxnDeletion(model2, 'FBA', plantRxns);
grRatioPlant_amm(isnan(grRatioPlant_amm)) = 0;

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

% Store the data as a matrix
temp = cell(length(plantRxnsOverall), 4);
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
        temp{pos,1} = round(grRatioPlant(n), 2);
        temp{pos,3} = round(grRatioPlant_amm(n), 2);
    elseif strmatch('Root', rxnOrig{1,1}, 'exact')
        temp{pos,2} = round(grRatioPlant(n), 2);
        temp{pos,4} = round(grRatioPlant_amm(n), 2);
    end
end
temp = horzcat(plantRxnsOverall, temp);
output_tissueSpecificRxnDeletion_plant = ...
    vertcat({'Gene','Shoot_nod','Root_nod','Shoot_amm','Root_amm'}, temp);
isEmpty = cellfun(@isempty, output_tissueSpecificRxnDeletion_plant);
output_tissueSpecificRxnDeletion_plant(isEmpty) = {1};