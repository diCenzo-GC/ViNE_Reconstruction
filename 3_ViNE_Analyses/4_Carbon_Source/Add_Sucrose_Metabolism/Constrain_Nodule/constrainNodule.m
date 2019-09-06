%% Load data

load('nodulatedPlant.mat');
rna_medic = table2cell(readtable('Medicago_TPM_values_average.txt'));
rna_sino = table2cell(readtable('Sinorhizobium_TPM_values_average.txt'));
rna_medic_full = table2cell(readtable('Medicago_TPM_values.txt'));
rna_sino_full = table2cell(readtable('Sinorhizobium_TPM_values.txt'));
rna_medic_original = rna_medic;
rna_sino_original = rna_sino;

%% Increase flux

nodulatedPlant.ub = nodulatedPlant.ub * 1000;
nodulatedPlant.lb = nodulatedPlant.lb * 1000;
optimizeCbModel(nodulatedPlant)

%% Set model constraints

% Determine O2 consumption limit of nodule zone III
nodulatedPlant.ub(findRxnIDs(nodulatedPlant, 'TRN_OXYGEN-MOLECULE')) = 1000 * 649 * 0.02;
sol = optimizeCbModel(nodulatedPlant);
nodulatedPlant.ub(findRxnIDs(nodulatedPlant, 'TNIII_OXYGEN-MOLECULE')) = ...
    sol.x(findRxnIDs(nodulatedPlant, 'TNIII_OXYGEN-MOLECULE'));
nodulatedPlant.ub(findRxnIDs(nodulatedPlant, 'TRN_OXYGEN-MOLECULE')) = 1000000;
optimizeCbModel(nodulatedPlant)

% Remove unwanted NoduleIII_EXCT reactions
toRemove = {'NoduleIII_EXCT_for_MNXM621_e0';...
    'NoduleIII_EXCT_for_MNXM1503_e0';'NoduleIII_EXCT_for_MNXM198_e0';...
    'NoduleIII_EXCT_for_MNXM165_e0';'NoduleIII_EXCT_for_MNXM468_e0';...
    'NoduleIII_EXCT_for_MNXM615_e0'};
nodulatedPlant = tncore_remove_reactions(nodulatedPlant, toRemove);
sol = optimizeCbModel(nodulatedPlant);
reactions = nodulatedPlant.rxns(strmatch('NoduleIII_EXCT_for', nodulatedPlant.rxns));
toRemove = {};
for n = 1:length(reactions)
    pos = findRxnIDs(nodulatedPlant, reactions{n});
    if abs(sol.x(pos)) < 0.00001
        toRemove = vertcat(toRemove, reactions{n});
    end
end
toRemove = setdiff(toRemove, {'NoduleIII_EXCT_for_MNXM167_e0'});
nodulatedPlant = tncore_remove_reactions(nodulatedPlant, toRemove);

%% Set growth rate threshold

sol = optimizeCbModel(nodulatedPlant);
growthThresh = 0.99 * sol.f;

%% Prepare the RNA-seq data

processRNA;

%% Constrain the nodule

% Turn on tiger
start_tiger('cplex');

% Produce tiger model
nodulatedPlant_tiger = cobra_to_tiger(nodulatedPlant);

% Split the rnaseq files
medicago_rnaseq = cell2mat(medicago_rnaseq_data_model(:,2));
meliloti_rnaseq = cell2mat(meliloti_rnaseq_data_model(:,2));
medicago_genes = medicago_rnaseq_data_model(:,1);
meliloti_genes = meliloti_rnaseq_data_model(:,1);

% Force FixNOQP to be on
meliloti_rnaseq(strmatch('BacteroidIII_sma0767', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('BacteroidIII_sma0767', meliloti_genes)) * 10;
meliloti_rnaseq(strmatch('BacteroidIII_sma1214', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('BacteroidIII_sma1214', meliloti_genes)) * 10;
meliloti_rnaseq(strmatch('BacteroidIZ_sma0767', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('BacteroidIZ_sma0767', meliloti_genes)) * 10;
meliloti_rnaseq(strmatch('BacteroidIZ_sma1214', meliloti_genes)) = ...
    meliloti_rnaseq(strmatch('BacteroidIZ_sma1214', meliloti_genes)) * 10;

% Prepare the input for tncore_multiGIMME
fields = {'Medicago'; 'Meliloti'};
expressStruct = struct();
genesStruct = struct();
threshStruct = struct();
expressStruct.Medicago = medicago_rnaseq;
genesStruct.Medicago = medicago_genes;
threshStruct.Medicago = medic_thresh;
expressStruct.Meliloti = meliloti_rnaseq;
genesStruct.Meliloti = meliloti_genes;
threshStruct.Meliloti = sino_thresh;

save('temp_beforeGIMME.mat');

% Run the modified version of gimme
[geneStates, genesOut, sol, tiger, weightStruct] = tncore_multi_gimme(nodulatedPlant_tiger, ...
    fields, expressStruct, genesStruct, threshStruct, 0.99);

save('temp_afterGIMME.mat');

%% Build context specific model in COBRA based on active rxns

buildContextModel;

%% Fix the rules and grRules
% Just in case, may not be necessary

% Change model name
model = finalNodulatedPlant;

% Fix the grRules string formatting
for n = 1:length(model.grRules)
    if isstring(model.grRules{n})
        model.grRules{n} = str2mat(model.grRules{n});
    else
        model.grRules{n} = model.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(model.rules)
    if isstring(model.rules{n})
        model.rules{n} = str2mat(model.rules{n});
    else
        model.rules{n} = model.rules{n};
    end
end

% Fix the grRules
for n = 1:length(model.grRules)
    if ~isempty(model.grRules{n})
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, ' or or', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(model.grRules{n}, 'or or ', '');
        model.grRules{n} = strrep(strrep(strrep(model.grRules{n}, ' or ', 'AAA'), ' or', ''),'AAA', ' or ');
        model.grRules{n} = strrep(strrep(strrep(model.grRules{n}, ' or ', 'AAA'), 'or ', ''),'AAA', ' or ');
    end
end

% Fix the rules
for n = 1:length(model.rules)
    if ~isempty(model.rules{n})
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | | ', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, ' | |', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(model.rules{n}, '| | ', '');
        model.rules{n} = strrep(strrep(strrep(model.rules{n}, ' | ', 'AAA'), '| ', ''),'AAA', ' | ');
    end
end

% Fix the grRules string formatting
for n = 1:length(model.grRules)
    if isstring(model.grRules{n})
        model.grRules{n} = str2mat(model.grRules{n});
    else
        model.grRules{n} = model.grRules{n};
    end
end

% Fix the rules string formatting
for n = 1:length(model.rules)
    if isstring(model.rules{n})
        model.rules{n} = str2mat(model.rules{n});
    else
        model.rules{n} = model.rules{n};
    end
end

% Change model name
constrainedNodule = model;

%% Save and clean workspace

save('allWorkspace.mat');
save('constrainedNodule_sucrose.mat', 'constrainedNodule');
clear;

