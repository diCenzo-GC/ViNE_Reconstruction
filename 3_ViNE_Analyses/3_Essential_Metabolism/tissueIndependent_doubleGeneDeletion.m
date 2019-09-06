%% Zone-independent synthetic lethal gene pairs

%% Symbiotic

% Get plant and bacteroid genes
isBacteroid = {};
for n = 1:length(model.genes)
    if strmatch('Bacteroid', model.genes{n});
        isBacteroid{n} = 1;
    else
        isBacteroid{n} = 0;
    end
end
plantGenes = model.genes(~logical(cell2mat(isBacteroid)));
bacteroidGenes = model.genes(logical(cell2mat(isBacteroid)));

% Get list of zone-independent plant genes
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

% Get list of zone-independent bacteroid genes
bacteroidGenesOverall = cell(length(bacteroidGenes), 1);
for n = 1:length(bacteroidGenes)
    temp = strsplit(bacteroidGenes{n}, '_');
    bacteroidGenesOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            bacteroidGenesOverall{n} = strcat(bacteroidGenesOverall{n}, '_', temp{m});
        end
    end
end
bacteroidGenesOverall = unique(bacteroidGenesOverall);

% Get plant gene pairings
plantGenePairs = cell(length(plantGenes), 2);
for n = 1:length(plantGenes)
    temp = strsplit(plantGenes{n}, '_');
    plantGenePairs{n,1} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            plantGenePairs{n,1} = strcat(plantGenePairs{n,1}, '_', temp{m});
        end
    end
    plantGenePairs{n,2} = plantGenes{n};
end

% Get bacteroid gene pairings
bacteroidGenePairs = cell(length(bacteroidGenes), 2);
for n = 1:length(bacteroidGenes)
    temp = strsplit(bacteroidGenes{n}, '_');
    bacteroidGenePairs{n,1} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            bacteroidGenePairs{n,1} = strcat(bacteroidGenePairs{n,1}, '_', temp{m});
        end
    end
    bacteroidGenePairs{n,2} = bacteroidGenes{n};
end

% Get the original solution during symbiosis
solOrig = optimizeCbModel(model);

% Plant genes to analyze (i.e., nonlethal)
nonlethalPlantGenes = {};
x = 0;
for n = 1:length(plantGenesOverall)
    pos = strmatch(plantGenesOverall{n}, plantGenePairs(:,1), 'exact');
    geneGroup = plantGenePairs(pos,2);
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, geneGroup);
    if isempty(constrRxns)
        x = x + 1;
        nonlethalPlantGenes{x,1} = plantGenesOverall{n};
        nonlethalPlantGenes{x,2} = 1;
    else
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) > 0.9
            x = x + 1;
            nonlethalPlantGenes{x,1} = plantGenesOverall{n};
            nonlethalPlantGenes{x,2} = sol.f / solOrig.f;
        end
    end
end
nonlethalPlantGenes = sortrows(nonlethalPlantGenes, 1);

% Bacteroid genes to analyze (i.e., nonlethal)
nonlethalBacteroidGenes = {};
x = 0;
for n = 1:length(bacteroidGenesOverall)
    pos = strmatch(bacteroidGenesOverall{n}, bacteroidGenePairs(:,1), 'exact');
    geneGroup = bacteroidGenePairs(pos,2);
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, geneGroup);
    if isempty(constrRxns)
        x = x + 1;
        nonlethalBacteroidGenes{x,1} = bacteroidGenesOverall{n};
        nonlethalBacteroidGenes{x,2} = 1;
    else
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) > 0.9
            x = x + 1;
            nonlethalBacteroidGenes{x,1} = bacteroidGenesOverall{n};
            nonlethalBacteroidGenes{x,2} = sol.f / solOrig.f;
        end
    end
end
nonlethalBacteroidGenes = sortrows(nonlethalBacteroidGenes, 1);

% Output variables
output_doubleGeneDeletion_plantAll = {};
output_doubleGeneDeletion_bacteroidAll = {};
output_doubleGeneDeletion_overall = {};

% Plant double gene deletion analysis
parpool(15);
parfor n = 1:length(nonlethalPlantGenes)-1
    changeCobraSolver('ibm_cplex')
    posA = strmatch(nonlethalPlantGenes{n,1}, plantGenePairs(:,1), 'exact');
    geneGroupA = plantGenePairs(posA,2);   
    for m = n+1:length(nonlethalPlantGenes)
        posB = strmatch(nonlethalPlantGenes{m,1}, plantGenePairs(:,1), 'exact');
        geneGroupB = plantGenePairs(posB,2);   
        genesToDelete = vertcat(geneGroupA, geneGroupB);
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        if ~isempty(constrRxns)
            modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
            sol = optimizeCbModel(modelTemp);
            if round(sol.f / solOrig.f, 3) < 0.01
                output_temp = {};
                output_temp{1,1} = strcat(nonlethalPlantGenes{n,1}, '__', nonlethalPlantGenes{m,1});
                output_temp{1,2} = nonlethalPlantGenes{n,1};
                output_temp{1,3} = nonlethalPlantGenes{m,1};
                output_temp{1,4} = nonlethalPlantGenes{n,2};
                output_temp{1,5} = nonlethalPlantGenes{m,2};
                output_temp{1,6} = sol.f / solOrig.f;
                output_doubleGeneDeletion_plantAll = vertcat(output_doubleGeneDeletion_plantAll, output_temp);
            end
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

% Bacteroid double gene deletion analysis
parpool(15);
parfor n = 1:length(nonlethalBacteroidGenes)-1
    changeCobraSolver('ibm_cplex')
    posA = strmatch(nonlethalBacteroidGenes{n,1}, bacteroidGenePairs(:,1), 'exact');
    geneGroupA = bacteroidGenePairs(posA,2);   
    for m = n+1:length(nonlethalBacteroidGenes)
        posB = strmatch(nonlethalBacteroidGenes{m,1}, bacteroidGenePairs(:,1), 'exact');
        geneGroupB = bacteroidGenePairs(posB,2);   
        genesToDelete = vertcat(geneGroupA, geneGroupB);
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        if ~isempty(constrRxns)
            modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
            sol = optimizeCbModel(modelTemp);
            if round(sol.f / solOrig.f, 3) < 0.01
                output_temp = {};
                output_temp{1,1} = strcat(nonlethalBacteroidGenes{n,1}, '__', nonlethalBacteroidGenes{m,1});
                output_temp{1,2} = nonlethalBacteroidGenes{n,1};
                output_temp{1,3} = nonlethalBacteroidGenes{m,1};
                output_temp{1,4} = nonlethalBacteroidGenes{n,2};
                output_temp{1,5} = nonlethalBacteroidGenes{m,2};
                output_temp{1,6} = sol.f / solOrig.f;
                output_doubleGeneDeletion_bacteroidAll = vertcat(output_doubleGeneDeletion_bacteroidAll, output_temp);
            end
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

% Plant-bacteroid double gene deletion analysis
parpool(15);
parfor n = 1:length(nonlethalPlantGenes)
    changeCobraSolver('ibm_cplex')
    posA = strmatch(nonlethalPlantGenes{n,1}, plantGenePairs(:,1), 'exact');
    geneGroupA = plantGenePairs(posA,2);   
    for m = 1:length(nonlethalBacteroidGenes)
        posB = strmatch(nonlethalBacteroidGenes{m,1}, bacteroidGenePairs(:,1), 'exact');
        geneGroupB = bacteroidGenePairs(posB,2);   
        genesToDelete = vertcat(geneGroupA, geneGroupB);
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        if ~isempty(constrRxns)
            modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
            sol = optimizeCbModel(modelTemp);
            if round(sol.f / solOrig.f, 3) < 0.01
                output_temp = {};
                output_temp{1,1} = strcat(nonlethalPlantGenes{n,1}, '__', nonlethalBacteroidGenes{m,1});
                output_temp{1,2} = nonlethalPlantGenes{n,1};
                output_temp{1,3} = nonlethalBacteroidGenes{m,1};
                output_temp{1,4} = nonlethalPlantGenes{n,2};
                output_temp{1,5} = nonlethalBacteroidGenes{m,2};
                output_temp{1,6} = sol.f / solOrig.f;
                output_doubleGeneDeletion_overall = vertcat(output_doubleGeneDeletion_overall, output_temp);
            end
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

%% Non-symbiotic

model2 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model2 = changeObjective(model2, 'PlantBiomass');
model2 = addExchangeRxn(model2, {'BiomassPlant[c]'}, 0, 1000000);
noduleReactions = model.rxns(strmatch('Nodule', model.rxns));
bacteroidReactions = model.rxns(strmatch('Bacteroid', model.rxns));
model2 = changeRxnBounds(model2, noduleReactions, 0, 'b');
model2 = changeRxnBounds(model2, bacteroidReactions, 0, 'b');

% Get plant and bacteroid genes
isBacteroid = {};
for n = 1:length(model2.genes)
    if strmatch('Bacteroid', model2.genes{n});
        isBacteroid{n} = 1;
    else
        isBacteroid{n} = 0;
    end
end
plantGenes = model2.genes(~logical(cell2mat(isBacteroid)));
bacteroidGenes = model2.genes(logical(cell2mat(isBacteroid)));

% Get list of zone-independent plant genes
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

% Get plant gene pairings
plantGenePairs = cell(length(plantGenes), 2);
for n = 1:length(plantGenes)
    temp = strsplit(plantGenes{n}, '_');
    plantGenePairs{n,1} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            plantGenePairs{n,1} = strcat(plantGenePairs{n,1}, '_', temp{m});
        end
    end
    plantGenePairs{n,2} = plantGenes{n};
end

% Get the original solution during symbiosis
solOrig = optimizeCbModel(model2);

% Plant genes to analyze (i.e., nonlethal)
nonlethalPlantGenes = {};
x = 0;
for n = 1:length(plantGenesOverall)
    pos = strmatch(plantGenesOverall{n}, plantGenePairs(:,1), 'exact');
    geneGroup = plantGenePairs(pos,2);
    [modelTemp, ~, constrRxns] = deleteModelGenes(model2, geneGroup);
    if isempty(constrRxns)
        x = x + 1;
        nonlethalPlantGenes{x,1} = plantGenesOverall{n};
        nonlethalPlantGenes{x,2} = 1;
    else
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) > 0.9
            x = x + 1;
            nonlethalPlantGenes{x,1} = plantGenesOverall{n};
            nonlethalPlantGenes{x,2} = sol.f / solOrig.f;
        end
    end
end
nonlethalPlantGenes = sortrows(nonlethalPlantGenes, 1);

% Output variables
output_doubleGeneDeletionAmm_plantAll = {};

% Plant double gene deletion analysis
parpool(15);
parfor n = 1:length(nonlethalPlantGenes)-1
    changeCobraSolver('ibm_cplex')
    posA = strmatch(nonlethalPlantGenes{n,1}, plantGenePairs(:,1), 'exact');
    geneGroupA = plantGenePairs(posA,2);   
    for m = n+1:length(nonlethalPlantGenes)
        posB = strmatch(nonlethalPlantGenes{m,1}, plantGenePairs(:,1), 'exact');
        geneGroupB = plantGenePairs(posB,2);   
        genesToDelete = vertcat(geneGroupA, geneGroupB);
        [modelTemp, ~, constrRxns] = deleteModelGenes(model2, genesToDelete);
        if ~isempty(constrRxns)
            modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
            sol = optimizeCbModel(modelTemp);
            if round(sol.f / solOrig.f, 3) < 0.01
                output_temp = {};
                output_temp{1,1} = strcat(nonlethalPlantGenes{n,1}, '__', nonlethalPlantGenes{m,1});
                output_temp{1,2} = nonlethalPlantGenes{n,1};
                output_temp{1,3} = nonlethalPlantGenes{m,1};
                output_temp{1,4} = nonlethalPlantGenes{n,2};
                output_temp{1,5} = nonlethalPlantGenes{m,2};
                output_temp{1,6} = sol.f / solOrig.f;
                output_doubleGeneDeletionAmm_plantAll = vertcat(output_doubleGeneDeletionAmm_plantAll, output_temp);
            end
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

