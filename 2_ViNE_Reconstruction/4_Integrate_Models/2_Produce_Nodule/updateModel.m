%% Load data

load('../1_Produce_Plant/Data/finalConversions.mat');

%% Update the model

% Save the input model
modelOrigSaved = model;

% Make new gene list
model.genes = unique(finalConversions(:,3));

% Remake the grRules
for n = 1:length(model.grRules)
    if ~isempty(model.grRules{n});
        newGrRules = 'toDelete';
        grRules = strrep(model.grRules{n}, ' or ', ' ');
        rxnGenes = strsplit(grRules, ' ');
        x = 0;
        for m = 1:length(rxnGenes)
            pos = strmatch(rxnGenes{m}, finalConversions(:,1), 'exact');
            if ~isempty(pos)
                x = x + 1;
                for i = 1:length(pos)
                    if x == 1 && i == 1
                        newGrRules = finalConversions{pos(i), 3};
                    else
                        if strfind(newGrRules, finalConversions{pos(i), 3})
                        else
                            newGrRules = strcat(newGrRules, ...
                                [' or ' finalConversions{pos(i), 3}]);
                        end
                    end
                end
            end
        end
        model.grRules{n} = newGrRules;
    end
end

% Add back in the essential genes
pos = findRxnIDs(model, 'RXN__45__9944_H');
model.grRules{pos} = 'Medtr6g008180.1';
model.genes{end+1} = 'Medtr6g008180.1';
pos = findRxnIDs(model, 'RXN__45__7674_H');
model.grRules{pos} = 'Medtr7g069400.1';
model.genes{end+1} = 'Medtr7g069400.1';
pos = findRxnIDs(model, 'PLASTOQUINOL__45____45__PLASTOCYANIN__45__REDUCTASE__45__RXN_H');
model.grRules{pos} = 'Medtr2g021500.1';
model.genes{end+1} = 'Medtr2g021500.1';
pos = findRxnIDs(model, 'RXN0__45__6555_C');
model.grRules{pos} = 'Medtr7g037970.1';
model.genes{end+1} = 'Medtr7g037970.1';

% Update the rules
for n = 1:length(model.rxns)
    if isempty(model.grRules{n});
        model.rules{n} = model.grRules{n};
    elseif strmatch('toDelete', model.grRules{n}, 'exact')
        model.rules{n} = '';
    else
        model.rules{n} = model.grRules{n};
        grRules = strrep(model.grRules{n}, ' or ', ' ');
        rxnGenes = strsplit(grRules, ' ');
        for m = 1:length(rxnGenes)
            geneID = strmatch(rxnGenes{m}, model.genes, 'exact');
            model.rules{n} = strrep(model.rules{n}, rxnGenes{m}, ['x(' mat2str(geneID) ')']);
        end
        model.rules{n} = strrep(model.rules{n}, ' or ', ' | ');
    end
end

% Rebuild the rxnGeneMat
model.rxnGeneMat = cell(length(model.rxns),length(model.genes));
for n = 1:length(model.rxns)
    if ~isempty(model.rules{n})
        rules = strsplit(model.rules{n}, ' | ');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            model.rxnGeneMat{n,str2num(rules{m})} = 1;
        end
    end
end
isEmpty = cellfun('isempty',model.rxnGeneMat);
model.rxnGeneMat(isEmpty) = {0};
model.rxnGeneMat = cell2mat(model.rxnGeneMat);
model.rxnGeneMat = sparse(double(model.rxnGeneMat));

% Delete unwanted reactions
rxnsToDelete = {};
for n = 1:length(model.grRules)
    if strmatch('toDelete', model.grRules{n}, 'exact')
        rxnsToDelete = vertcat(rxnsToDelete, model.rxns{n});
    end
end
updatedModel = tncore_remove_reactions(model, rxnsToDelete);

