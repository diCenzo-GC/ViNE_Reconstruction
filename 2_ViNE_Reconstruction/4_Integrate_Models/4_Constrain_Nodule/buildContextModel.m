%% Build initial core model based on GIMME

% Save a copy of the model before further modification
originalModelAfterGimme = nodulatedPlant;

% Get growth threshold
newSol = optimizeCbModel(nodulatedPlant)
growthThresh = newSol.f * 0.99;

% Adjust RNAseq values
medicago_rnaseq = medicago_rnaseq / weightStruct.Medicago;
meliloti_rnaseq = meliloti_rnaseq / weightStruct.Meliloti;
expressThresh = medic_thresh;

% Make combined RNAseq variable and link to genes
allRnaSeq = vertcat(num2cell(medicago_rnaseq), num2cell(meliloti_rnaseq));
allGenes = vertcat(medicago_genes, meliloti_genes);
modelRNAseq = cell(length(nodulatedPlant.genes), 1);
for n = 1:length(nodulatedPlant.genes)
    pos = strmatch(nodulatedPlant.genes{n}, allGenes);
    modelRNAseq{n} = allRnaSeq{pos};
end
for n = 1:length(nodulatedPlant.genes)
    geneID = findGeneIDs(nodulatedPlant, nodulatedPlant.genes{n,1});
    genesRNAseq{n,1} = nodulatedPlant.genes{n,1};
    genesRNAseq{n,2} = modelRNAseq{geneID,1};
end
genesRNAseq = sortrows(genesRNAseq,2);

% Identify genes off in GIMME and delete from the reduced model
isPresent = logical(geneStates);
genesToDelete = genesOut(~isPresent);
[testModel, ~, ~] = deleteModelGenes(nodulatedPlant, genesToDelete);

% Test if the resulting model grows
testSol = optimizeCbModel(testModel)

% Identify reactions active in GIMME
isActive = {};
for n = 1:length(nodulatedPlant.rxns)
    if sol.x(n) > 1e-6 || sol.x(n) < -1e-6
        isActive{n} = 1;
    else
        isActive{n} = 0;
    end
end
isActive = logical(cell2mat(isActive));

% Find a minimal consistent rxn set based on the active rxns in GIMME
C = findRxnIDs(nodulatedPlant, nodulatedPlant.rxns(isActive));
nodulatedPlant.lb(findRxnIDs(nodulatedPlant, 'Biomass')) = growthThresh;
A = fastcore(C, nodulatedPlant, 1.01e-6);
nodulatedPlant.lb(findRxnIDs(nodulatedPlant, 'Biomass')) = 0;
coreRxns = nodulatedPlant.rxns(A);

% Get rxns associated with active genes and that do not require inactive genes
contextGenes = genesOut(logical(geneStates));
contextGenesInverse = genesOut(~logical(geneStates));
[~, contextRxns] = findRxnsFromGenes(nodulatedPlant, contextGenes, [], 1);
contextRxns = unique(contextRxns(:,1));
[~, ~, contextRxnsInverse] = deleteModelGenes(nodulatedPlant, contextGenesInverse);
contextRxns = setdiff(contextRxns, contextRxnsInverse);

% Combine the context and core reactions
allRxns = unique(vertcat(coreRxns, contextRxns));

% Add to this all the root and shoot
rootRxns = nodulatedPlant.rxns(strmatch('Root_', nodulatedPlant.rxns));
shootRxns = nodulatedPlant.rxns(strmatch('Leave_', nodulatedPlant.rxns));
transfer1Rxns = nodulatedPlant.rxns(strmatch('TRS_', nodulatedPlant.rxns));
transfer2Rxns = nodulatedPlant.rxns(strmatch('TSR_', nodulatedPlant.rxns));
transfer3Rxns = nodulatedPlant.rxns(strmatch('TRN_', nodulatedPlant.rxns));
transfer4Rxns = nodulatedPlant.rxns(strmatch('TNI_', nodulatedPlant.rxns));
allRxns = unique(vertcat(allRxns, rootRxns));
allRxns = unique(vertcat(allRxns, shootRxns));
allRxns = unique(vertcat(allRxns, transfer1Rxns));
allRxns = unique(vertcat(allRxns, transfer2Rxns));
allRxns = unique(vertcat(allRxns, transfer3Rxns));
allRxns = unique(vertcat(allRxns, transfer4Rxns));

% Add to this all peribacteroid transport reactions
exctRxnsI = nodulatedPlant.rxns(strmatch('NoduleI_EXCT_', nodulatedPlant.rxns));
exctRxnsIId = nodulatedPlant.rxns(strmatch('NoduleIId_EXCT_', nodulatedPlant.rxns));
exctRxnsIIp = nodulatedPlant.rxns(strmatch('NoduleIIp_EXCT_', nodulatedPlant.rxns));
exctRxnsIZ = nodulatedPlant.rxns(strmatch('NoduleIZ_EXCT_', nodulatedPlant.rxns));
exctRxnsIII = nodulatedPlant.rxns(strmatch('NoduleIII_EXCT_', nodulatedPlant.rxns));
allRxns = unique(vertcat(allRxns, exctRxnsI));
allRxns = unique(vertcat(allRxns, exctRxnsIId));
allRxns = unique(vertcat(allRxns, exctRxnsIIp));
allRxns = unique(vertcat(allRxns, exctRxnsIZ));
allRxns = unique(vertcat(allRxns, exctRxnsIII));

% Ensure all proton exporters remain in the nodule
protonExport = {'NoduleI_TCE_PROTON'; 'NoduleIId_TCE_PROTON'; 'NoduleIIp_TCE_PROTON';...
    'NoduleIZ_TCE_PROTON'; 'NoduleIII_TCE_PROTON'};
allRxns = unique(vertcat(allRxns, protonExport));

% Prepare the reduced model
nodulatedPlant_reduced = tncore_remove_reactions(nodulatedPlant, setdiff(nodulatedPlant.rxns, allRxns));
nodulatedPlant_reduced = tncore_remove(nodulatedPlant_reduced);
model = nodulatedPlant_reduced;

%% Remove unnecessary genes based on expression

% Force 'on' genes to be kept
offGenes = genesOut(~logical(geneStates));
offGenes = intersect(offGenes, nodulatedPlant_reduced.genes);
onGenes = setdiff(nodulatedPlant_reduced.genes, offGenes);
for n = 1:length(genesRNAseq)
    if strmatch(genesRNAseq{n,1}, onGenes, 'exact')
        genesRNAseq{n,2} = 2 * expressThresh;
    end
end

% Identify unnecessary 'off' to remove
for n = 1:length(model.rxns)
    genesToKeep = {};
    genesToKeepTemp = {};
    genesToRemove = {};
    genesToKeep_inv = {};
    genesToKeep_invNew = {};
    if strfind(model.grRules{n}, ' or ')
        if strfind(model.grRules{n}, ' and ')
            grRulesSplit = transpose(strsplit(model.grRules{n}, ' and '));
            grRulesSplitTemp = {};
            for m = 1:length(grRulesSplit)
                if strfind(grRulesSplit{m}, ') or (')
                else
                    if strfind(grRulesSplit{m}, ' or ')
                        grRulesSplitTemp = vertcat(grRulesSplitTemp, grRulesSplit{m});
                    end
                end
            end
            grRulesSplit = grRulesSplitTemp;
            for m = 1:length(grRulesSplit)
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, '(', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ')', '');
            end
            for m = 1:length(grRulesSplit)
                if strfind(grRulesSplit{m}, ' or ')
                    grRulesSplit2 = transpose(strsplit(grRulesSplit{m}, ' or '));
                    for i = 1:length(grRulesSplit2)
                        grRulesSplit2{i,2} = genesRNAseq{...
                            strmatch(grRulesSplit2{i,1}, genesRNAseq(:,1), 'exact'), 2};
                    end
                end
                grRulesSplit2 = sortrows(grRulesSplit2, -2);
                if max(cell2mat(grRulesSplit2(:,2))) >= expressThresh
                    for i = 1:length(grRulesSplit2)
                        if grRulesSplit2{i,2} >= expressThresh
                        else
                            genesToKeep_invNew = vertcat(genesToKeep_invNew, grRulesSplit2{i,1});
                        end
                    end
                else
                    genesToKeep_invNew = vertcat(genesToKeep_invNew, grRulesSplit2(2:end,1));
                end
            end
            tempGrRules = strrep(model.grRules{n}, ' and ', ' or ');
            grRulesSplit = transpose(strsplit(tempGrRules, ' or '));
            for m = 1:length(grRulesSplit)
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, '(', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ')', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ' ', '');
                grRulesSplit{m,2} = genesRNAseq{...
                    strmatch(grRulesSplit{m,1}, genesRNAseq(:,1), 'exact'), 2};
            end
            grRulesSplit = sortrows(grRulesSplit, -2);
            genesToKeep = {};
            genesToKeep_inv = {};
            for m = 1:length(grRulesSplit)
                if grRulesSplit{m,2} >= expressThresh
                    genesToKeep = vertcat(genesToKeep, grRulesSplit{m,1});
                else
                    genesToKeep_inv = vertcat(genesToKeep_inv, grRulesSplit{m,1});
                end
            end
            genesToKeep_inv = setdiff(genesToKeep_inv, genesToKeep_invNew);
            [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_inv);
            if ~isempty(constrRxns)
                if strmatch(model.rxns{n}, constrRxns, 'exact');
                    genesToKeepTemp = genesToKeep;
                    for m = length(genesToKeep)+1:length(grRulesSplit)
                        if strmatch(grRulesSplit{m,1}, genesToKeep_invNew, 'exact')
                        else
                            genesToKeep_invTemp = setdiff(grRulesSplit(:,1), genesToKeepTemp);
                            genesToKeep_invTemp = vertcat(genesToKeep_invTemp, genesToKeep_invNew);
                            [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_invTemp);
                            if strmatch(model.rxns{n}, constrRxns, 'exact')
                                genesToKeepTemp = vertcat(genesToKeepTemp, grRulesSplit{m,1});
                                genesToKeep_invTemp = setdiff(grRulesSplit(:,1), genesToKeepTemp);
                                [~,~,constrRxns] = deleteModelGenes(model, genesToKeep_invTemp);
                                if isempty(constrRxns)
                                    if ~isempty(genesToKeep)
                                        tempGenes = vertcat(genesToKeep, genesToKeep_invTemp);
                                        [~,~,constrRxns] = deleteModelGenes(model, tempGenes);
                                        if isempty(constrRxns)
                                            genesToKeepTemp = genesToKeepTemp(1:end-1);
                                        elseif ~strmatch(model.rxns{n}, constrRxns, 'exact')
                                            genesToKeepTemp = genesToKeepTemp(1:end-1);
                                        end
                                    end
                                elseif ~strmatch(model.rxns{n}, constrRxns, 'exact')
                                    if ~isempty(genesToKeep)
                                        tempGenes = vertcat(genesToKeep, genesToKeep_invTemp);
                                        [~,~,constrRxns] = deleteModelGenes(model, tempGenes);
                                        if isempty(constrRxns)
                                            genesToKeepTemp = genesToKeepTemp(1:end-1);
                                        elseif ~strmach(model.rxns{n}, constrRxns, 'exact')
                                            genesToKeepTemp = genesToKeepTemp(1:end-1);
                                        end
                                    end
                                end
                            end
                        end
                    end
                    genesToKeep = genesToKeepTemp;
                    genesToRemove = setdiff(grRulesSplit(:,1), genesToKeep);
                end
                genesToKeep2 = cell(length(genesToKeep), 2);
                for m = 1:length(genesToKeep);
                    pos = strmatch(genesToKeep{m,1}, grRulesSplit(:,1), 'exact');
                    genesToKeep2{m,1} = genesToKeep{m,1};
                    genesToKeep2{m,2} = grRulesSplit{pos, 2};
                end
                genesToKeep = sortrows(genesToKeep2, 2);
                genesToRemoveTemp = genesToRemove;
                for m = 1:size(genesToKeep,1)
                    genesToRemoveTemp = vertcat(genesToRemoveTemp, genesToKeep{m,1});
                    [~,~,constrRxns] = deleteModelGenes(model, genesToRemoveTemp);
                    if strmatch(model.rxns{n}, constrRxns, 'exact')
                        genesToRemoveTemp = genesToRemoveTemp(1:end-1);
                    end
                end
                genesToRemove = genesToRemoveTemp;
            else
                genesToRemove = genesToKeep_inv;
            end
        else
            grRulesSplit = transpose(strsplit(model.grRules{n}, ' or '));
            for m = 1:length(grRulesSplit)
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, '(', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ')', '');
                grRulesSplit{m,1} = strrep(grRulesSplit{m,1}, ' ', '');
                grRulesSplit{m,2} = genesRNAseq{...
                    strmatch(grRulesSplit{m,1}, genesRNAseq(:,1), 'exact'), 2};
            end
            grRulesSplit = sortrows(grRulesSplit, -2);
            genesToKeep = {};
            genesToRemove = {};
            if max(cell2mat(grRulesSplit(:,2))) >= expressThresh
                for m = 1:length(grRulesSplit)
                    if grRulesSplit{m,2} >= expressThresh
                        genesToKeep = vertcat(genesToKeep, grRulesSplit{m,1});
                    else
                        genesToRemove = vertcat(genesToRemove, grRulesSplit{m,1});
                    end
                end
            else
                genesToKeep = grRulesSplit{1,1};
                genesToRemove = grRulesSplit(2:end,1);
            end
        end
        if ~isempty(genesToRemove)
            model.rxnGeneMat = full(model.rxnGeneMat);
            for m = 1:length(genesToRemove)
                model.grRules{n} = strrep(model.grRules{n}, ...
                    genesToRemove{m}, 'toDelete');
                genePos = strmatch(genesToRemove{m}, model.genes, 'exact');
                ruleToRemove = ['x(' num2str(genePos) ')'];
                if strmatch('toDelete', model.genes, 'exact')
                    ruleToAdd = ['x(' num2str(length(model.genes)) ')'];
                else
                    ruleToAdd = ['x(' num2str(length(model.genes)+1) ')'];
                end
                model.rules{n} = strrep(model.rules{n}, ...
                    ruleToRemove, ruleToAdd);
                model.rxnGeneMat(n,genePos) = 0;
            end
            if strmatch('toDelete', model.genes, 'exact')
                model.rxnGeneMat(n,length(model.genes)) = 1;
            else
                model.genes{end+1,1} = 'toDelete';
                model.rxnGeneMat(n,length(model.genes)) = 1;
            end
        end
    end
end
model.rxnGeneMat = num2cell(model.rxnGeneMat);
temp = cellfun('isempty',model.rxnGeneMat);
model.rxnGeneMat(temp) = {0};
model.rxnGeneMat = cell2mat(model.rxnGeneMat);
model.rxnGeneMat = sparse(double(model.rxnGeneMat));

%% Prepare the context-specific model

% Prepare the reduced model
model = tncore_remove(model);
[model, ~, constrRxns] = deleteModelGenes(model, {'toDelete'});
model = tncore_delete(model);

% Remove reactions producing deadend metabolites except in shoot and root
tempModel = tncore_deadends(model);
rxnsToRemove = setdiff(model.rxns, tempModel.rxns);
rxnsToRemoveNew = {};
for n = 1:length(rxnsToRemove)
    if strmatch('Root_', rxnsToRemove{n})
    elseif strmatch('Leave_', rxnsToRemove{n})
    else
        rxnsToRemoveNew = vertcat(rxnsToRemoveNew, rxnsToRemove{n});
    end
end
model2 = tncore_remove_reactions(model, rxnsToRemoveNew);
finalNodulatedPlant = tncore_remove(model2);
