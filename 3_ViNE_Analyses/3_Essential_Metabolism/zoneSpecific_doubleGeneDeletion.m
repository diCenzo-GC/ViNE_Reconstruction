%% Zone-specific synthetic lethal gene pairs

%% Zone I

% Get zone I genes
noduleGenes = model.genes(strmatch('NoduleI_', model.genes));
zoneGenes = vertcat(noduleGenes);

% Original solution
solOrig = optimizeCbModel(model);

% Genes to analyze (i.e., nonlethal)
nonlethalGenes = {};
x = 0;
for n = 1:length(zoneGenes)
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, zoneGenes{n});
    modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
    sol = optimizeCbModel(modelTemp);
    if round(sol.f / solOrig.f, 3) > 0.9
        x = x + 1;
        nonlethalGenes{x,1} = zoneGenes{n};
        nonlethalGenes{x,2} = sol.f / solOrig.f;
    end
end
nonlethalGenes = sortrows(nonlethalGenes, 1);

% Output variable
output_doubleGeneDeletion_zoneI = {};

% Perform deletion analysis
parpool(15);
parfor n = 1:length(nonlethalGenes)-1
    changeCobraSolver('ibm_cplex')
    for m = n+1:length(nonlethalGenes)
        genesToDelete = vertcat(nonlethalGenes(n,1), nonlethalGenes(m,1));
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) < 0.01
            output_temp = {};
            output_temp{1,1} = strcat(nonlethalGenes{n,1}, '__', nonlethalGenes{m,1});
            output_temp{1,2} = nonlethalGenes{n,1};
            output_temp{1,3} = nonlethalGenes{m,1};
            output_temp{1,4} = nonlethalGenes{n,2};
            output_temp{1,5} = nonlethalGenes{m,2};
            output_temp{1,6} = sol.f / solOrig.f;
            output_doubleGeneDeletion_zoneI = vertcat(output_doubleGeneDeletion_zoneI, output_temp);
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

%% Zone IId

% Get zone IId genes
noduleGenes = model.genes(strmatch('NoduleIId_', model.genes));
bacteroidGenes = model.genes(strmatch('BacteroidIId_', model.genes));
zoneGenes = vertcat(noduleGenes, bacteroidGenes);

% Original solution
solOrig = optimizeCbModel(model);

% Genes to analyze (i.e., nonlethal)
nonlethalGenes = {};
x = 0;
for n = 1:length(zoneGenes)
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, zoneGenes{n});
    modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
    sol = optimizeCbModel(modelTemp);
    if round(sol.f / solOrig.f, 3) > 0.9
        x = x + 1;
        nonlethalGenes{x,1} = zoneGenes{n};
        nonlethalGenes{x,2} = sol.f / solOrig.f;
    end
end
nonlethalGenes = sortrows(nonlethalGenes, 1);

% Output variable
output_doubleGeneDeletion_zoneIId = {};

% Perform deletion analysis
parpool(15);
parfor n = 1:length(nonlethalGenes)-1
    changeCobraSolver('ibm_cplex')
    for m = n+1:length(nonlethalGenes)
        genesToDelete = vertcat(nonlethalGenes(n,1), nonlethalGenes(m,1));
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) < 0.01
            output_temp = {};
            output_temp{1,1} = strcat(nonlethalGenes{n,1}, '__', nonlethalGenes{m,1});
            output_temp{1,2} = nonlethalGenes{n,1};
            output_temp{1,3} = nonlethalGenes{m,1};
            output_temp{1,4} = nonlethalGenes{n,2};
            output_temp{1,5} = nonlethalGenes{m,2};
            output_temp{1,6} = sol.f / solOrig.f;
            output_doubleGeneDeletion_zoneIId = vertcat(output_doubleGeneDeletion_zoneIId, output_temp);
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

%% Zone IIp

% Get zone IIp genes
noduleGenes = model.genes(strmatch('NoduleIIp_', model.genes));
bacteroidGenes = model.genes(strmatch('BacteroidIIp_', model.genes));
zoneGenes = vertcat(noduleGenes, bacteroidGenes);

% Original solution
solOrig = optimizeCbModel(model);

% Genes to analyze (i.e., nonlethal)
nonlethalGenes = {};
x = 0;
for n = 1:length(zoneGenes)
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, zoneGenes{n});
    modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
    sol = optimizeCbModel(modelTemp);
    if round(sol.f / solOrig.f, 3) > 0.9
        x = x + 1;
        nonlethalGenes{x,1} = zoneGenes{n};
        nonlethalGenes{x,2} = sol.f / solOrig.f;
    end
end
nonlethalGenes = sortrows(nonlethalGenes, 1);

% Output variable
output_doubleGeneDeletion_zoneIIp = {};

% Perform deletion analysis
parpool(15);
parfor n = 1:length(nonlethalGenes)-1
    changeCobraSolver('ibm_cplex')
    for m = n+1:length(nonlethalGenes)
        genesToDelete = vertcat(nonlethalGenes(n,1), nonlethalGenes(m,1));
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) < 0.01
            output_temp = {};
            output_temp{1,1} = strcat(nonlethalGenes{n,1}, '__', nonlethalGenes{m,1});
            output_temp{1,2} = nonlethalGenes{n,1};
            output_temp{1,3} = nonlethalGenes{m,1};
            output_temp{1,4} = nonlethalGenes{n,2};
            output_temp{1,5} = nonlethalGenes{m,2};
            output_temp{1,6} = sol.f / solOrig.f;
            output_doubleGeneDeletion_zoneIIp = vertcat(output_doubleGeneDeletion_zoneIIp, output_temp);
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

%% Zone IZ

% Get zone IZ genes
noduleGenes = model.genes(strmatch('NoduleIZ_', model.genes));
bacteroidGenes = model.genes(strmatch('BacteroidIZ_', model.genes));
zoneGenes = vertcat(noduleGenes, bacteroidGenes);

% Original solution
solOrig = optimizeCbModel(model);

% Genes to analyze (i.e., nonlethal)
nonlethalGenes = {};
x = 0;
for n = 1:length(zoneGenes)
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, zoneGenes{n});
    modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
    sol = optimizeCbModel(modelTemp);
    if round(sol.f / solOrig.f, 3) > 0.9
        x = x + 1;
        nonlethalGenes{x,1} = zoneGenes{n};
        nonlethalGenes{x,2} = sol.f / solOrig.f;
    end
end
nonlethalGenes = sortrows(nonlethalGenes, 1);

% Output variable
output_doubleGeneDeletion_zoneIZ = {};

% Perform deletion analysis
parpool(15);
parfor n = 1:length(nonlethalGenes)-1
    changeCobraSolver('ibm_cplex')
    for m = n+1:length(nonlethalGenes)
        genesToDelete = vertcat(nonlethalGenes(n,1), nonlethalGenes(m,1));
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) < 0.01
            output_temp = {};
            output_temp{1,1} = strcat(nonlethalGenes{n,1}, '__', nonlethalGenes{m,1});
            output_temp{1,2} = nonlethalGenes{n,1};
            output_temp{1,3} = nonlethalGenes{m,1};
            output_temp{1,4} = nonlethalGenes{n,2};
            output_temp{1,5} = nonlethalGenes{m,2};
            output_temp{1,6} = sol.f / solOrig.f;
            output_doubleGeneDeletion_zoneIZ = vertcat(output_doubleGeneDeletion_zoneIZ, output_temp);
        end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

%% Zone III

% Get zone III genes
noduleGenes = model.genes(strmatch('NoduleIII_', model.genes));
bacteroidGenes = model.genes(strmatch('BacteroidIII_', model.genes));
zoneGenes = vertcat(noduleGenes, bacteroidGenes);

% Original solution
solOrig = optimizeCbModel(model);

% Genes to analyze (i.e., nonlethal)
nonlethalGenes = {};
x = 0;
for n = 1:length(zoneGenes)
    [modelTemp, ~, constrRxns] = deleteModelGenes(model, zoneGenes{n});
    modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
    sol = optimizeCbModel(modelTemp);
    if round(sol.f / solOrig.f, 3) > 0.9
        x = x + 1;
        nonlethalGenes{x,1} = zoneGenes{n};
        nonlethalGenes{x,2} = sol.f / solOrig.f;
    end
end
nonlethalGenes = sortrows(nonlethalGenes, 1);

% Output variable
output_doubleGeneDeletion_zoneIII = {};

% Perform deletion analysis
parpool(15);
parfor n = 1:length(nonlethalGenes)-1
    changeCobraSolver('ibm_cplex')
    for m = n+1:length(nonlethalGenes)
        genesToDelete = vertcat(nonlethalGenes(n,1), nonlethalGenes(m,1));
        [modelTemp, ~, constrRxns] = deleteModelGenes(model, genesToDelete);
        modelTemp = tncore_remove_reactions(modelTemp, constrRxns);
        sol = optimizeCbModel(modelTemp);
        if round(sol.f / solOrig.f, 3) < 0.01
            output_temp = {};
            output_temp{1,1} = strcat(nonlethalGenes{n,1}, '__', nonlethalGenes{m,1});
            output_temp{1,2} = nonlethalGenes{n,1};
            output_temp{1,3} = nonlethalGenes{m,1};
            output_temp{1,4} = nonlethalGenes{n,2};
            output_temp{1,5} = nonlethalGenes{m,2};
            output_temp{1,6} = sol.f / solOrig.f;
            output_doubleGeneDeletion_zoneIII = vertcat(output_doubleGeneDeletion_zoneIII, output_temp);
       end
    end
end
delete(gcp('nocreate'))
save('temp_doubleDeletion2.mat');

