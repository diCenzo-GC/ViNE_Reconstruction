%% Prepare the models

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Make a nodule-free model
model2 = changeRxnBounds(model, 'Root_TEC_AMMONIUM', 1000000, 'u');
model2 = changeObjective(model2, 'PlantBiomass');
model2 = addExchangeRxn(model2, {'BiomassPlant[c]'}, 0, 1000000);
noduleReactions = model.rxns(strmatch('Nodule', model.rxns));
bacteroidReactions = model.rxns(strmatch('Bacteroid', model.rxns));
model2 = changeRxnBounds(model2, noduleReactions, 0, 'b');
model2 = changeRxnBounds(model2, bacteroidReactions, 0, 'b');

%% Prepare gene lists

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

% Get list of zone-independent nodule genes
plantGenesOverall = cell(length(plantGenes), 1);
for n = 1:length(plantGenes)
    temp = strsplit(plantGenes{n}, '_');
    plantGenesOverall{n} = temp{1,2};
end
plantGenesOverall = unique(plantGenesOverall);

% Get list of zone-independent bacteroid genes
bacteroidGenesOverall = cell(length(bacteroidGenes), 1);
for n = 1:length(bacteroidGenes)
    temp = strsplit(bacteroidGenes{n}, '_');
    bacteroidGenesOverall{n} = temp{1,2};
end
bacteroidGenesOverall = unique(bacteroidGenesOverall);
bacteroidGenesOverall = setdiff(bacteroidGenesOverall, {'BacteroidIId'; ...
    'BacteroidIIp'; 'BacteroidIZ'; 'BacteroidIII'; 'Spontaneous'; 'Unknown'});

% Get plant gene pairings
plantGenePairs = cell(length(plantGenes), 2);
for n = 1:length(plantGenes)
    temp = strsplit(plantGenes{n}, '_');
    plantGenePairs{n,1} = temp{1,2};
    plantGenePairs{n,2} = plantGenes{n};
end

% Get bacteroid gene pairings
bacteroidGenePairs = cell(length(bacteroidGenes), 2);
for n = 1:length(bacteroidGenes)
    temp = strsplit(bacteroidGenes{n}, '_');
    bacteroidGenePairs{n,1} = temp{1,2};
    bacteroidGenePairs{n,2} = bacteroidGenes{n};
end

%% Determine plant gene phenotypes with and without nodule

% Get the original solution during symbiosis
solOrig = optimizeCbModel(model);

% Plant deletion analysis during symbiosis
output_geneDeletion_plant = cell(length(plantGenesOverall), 3);
for n = 1:length(plantGenesOverall)
    pos = strmatch(plantGenesOverall{n}, plantGenePairs(:,1), 'exact');
    geneGroup = plantGenePairs(pos,2);
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, geneGroup);
    if isempty(constrRxns)
        output_geneDeletion_plant{n,1} = plantGenesOverall{n};
        output_geneDeletion_plant{n,2} = 1;
    else
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        output_geneDeletion_plant{n,1} = plantGenesOverall{n};
        if sol.f == 0
            output_geneDeletion_plant{n,2} = 0;
        else
            output_geneDeletion_plant{n,2} = round(sol.x(findRxnIDs(modelTemp, 'PlantBiomass')) / ...
                solOrig.x(findRxnIDs(model, 'PlantBiomass')), 2);
        end
    end
end

% Get the original solution with ammonium
solOrig_amm = optimizeCbModel(model2);

% Plant deletion analysis with ammonium
for n = 1:length(plantGenesOverall)
    pos = strmatch(plantGenesOverall{n}, plantGenePairs(:,1), 'exact');
    geneGroup = plantGenePairs(pos,2);
    [modelTemp, ~, constrRxns] = deleteModelGenes(model2, geneGroup);
    if isempty(constrRxns)
        output_geneDeletion_plant{n,1} = plantGenesOverall{n};
        output_geneDeletion_plant{n,3} = 1;
    else
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if sol.f == 0
            output_geneDeletion_plant{n,3} = 0;
        else
            output_geneDeletion_plant{n,3} = round(sol.x(findRxnIDs(modelTemp, 'PlantBiomass')) / ...
                solOrig_amm.x(findRxnIDs(model, 'PlantBiomass')), 2);
        end
    end
end

% Add headers
headers = {'Gene', 'Symbiosis', 'Ammonium'};
output_geneDeletion_plant = vertcat(headers, output_geneDeletion_plant);

%% Determine bacteroid gene phenotypes with nodule

% Get the original solution during symbiosis
solOrig = optimizeCbModel(model);

% Bacteroid deletion analysis during symbiosis
output_geneDeletion_bacteriod = cell(length(bacteroidGenesOverall), 2);
for n = 1:length(bacteroidGenesOverall)
    pos = strmatch(bacteroidGenesOverall{n}, bacteroidGenePairs(:,1), 'exact');
    geneGroup = bacteroidGenePairs(pos,2);
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, geneGroup);
    if isempty(constrRxns)
        output_geneDeletion_bacteriod{n,1} = bacteroidGenesOverall{n};
        output_geneDeletion_bacteriod{n,2} = 1;
    else
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        output_geneDeletion_bacteriod{n,1} = bacteroidGenesOverall{n};
        if sol.f == 0
            output_geneDeletion_bacteriod{n,2} = 0;
        else
            output_geneDeletion_bacteriod{n,2} = round(sol.x(findRxnIDs(modelTemp, 'PlantBiomass')) / ...
                solOrig.x(findRxnIDs(model, 'PlantBiomass')), 2);
        end
    end
end

% Add headers
headers = {'Gene', 'Symbiosis'};
output_geneDeletion_bacteriod = vertcat(headers, output_geneDeletion_bacteriod);

%% Clean the workspace

% Clear unwanted variables
clearvars -except output_*
