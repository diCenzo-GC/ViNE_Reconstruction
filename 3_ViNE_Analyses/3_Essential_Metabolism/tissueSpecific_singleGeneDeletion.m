%% Tissue-specific plant gene deletion

% Get plant genes
shootGenes = model.genes(strmatch('Leave', model.genes));
rootGenes = model.genes(strmatch('Root', model.genes));
plantGenes = vertcat(shootGenes, rootGenes);

% Perform gene deletion analyses
grRatioPlant = singleGeneDeletion(model, 'FBA', plantGenes);
grRatioPlant(isnan(grRatioPlant)) = 0;

% Perform gene deletion analysis with ammonium
model2 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model2 = changeObjective(model2, 'PlantBiomass');
model2 = addExchangeRxn(model2, {'BiomassPlant[c]'}, 0, 1000000);
noduleReactions = model.rxns(strmatch('Nodule', model.rxns));
bacteroidReactions = model.rxns(strmatch('Bacteroid', model.rxns));
model2 = changeRxnBounds(model2, noduleReactions, 0, 'b');
model2 = changeRxnBounds(model2, bacteroidReactions, 0, 'b');
grRatioPlant_amm = singleGeneDeletion(model2, 'FBA', plantGenes);
grRatioPlant_amm(isnan(grRatioPlant_amm)) = 0;

% Get list of tissue-independent plant genes
plantGenesOverall = cell(length(plantGenes), 1);
for n = 1:length(plantGenes)
    temp = strsplit(plantGenes{n}, '_');
    plantGenesOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            plantGenesOverall{n} = strcat(plantGenesOverall{n}, '_', temp{m});
        end
    end
end
plantGenesOverall = unique(plantGenesOverall);

% Store the data as a matrix
temp = cell(length(plantGenesOverall), 4);
for n = 1:length(plantGenes)
    geneOrig = strsplit(plantGenes{n}, '_');
    temporary = geneOrig{1,2};
    if length(geneOrig) > 2
        for m = 3:length(geneOrig)
            temporary = strcat(temporary, '_', geneOrig{m});
        end
    end
    pos = strmatch(temporary, plantGenesOverall, 'exact');
    if strmatch('Leave', geneOrig{1,1}, 'exact')
        temp{pos,1} = round(grRatioPlant(n), 2);
        temp{pos,3} = round(grRatioPlant_amm(n), 2);
    elseif strmatch('Root', geneOrig{1,1}, 'exact')
        temp{pos,2} = round(grRatioPlant(n), 2);
        temp{pos,4} = round(grRatioPlant_amm(n), 2);
    end
end
temp = horzcat(plantGenesOverall, temp);
output_tissueSpecificGeneDeletion_plant = ...
    vertcat({'Gene','Shoot_nod','Root_nod','Shoot_amm','Root_amm'}, temp);
isEmpty = cellfun(@isempty, output_tissueSpecificGeneDeletion_plant);
output_tissueSpecificGeneDeletion_plant(isEmpty) = {1};