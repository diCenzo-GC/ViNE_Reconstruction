%% Load data

% Load the 3.5 to 4.0 conversion table
conversionTable = table2cell(readtable('Data/Mt3.5-Mt4.0v1_conversion_table.txt'));

% Load the 4.0 to 5.0 conversion table
conversionTable2 = table2cell(readtable('Data/MtrunA17r5.0-ANR-EGN-r1.6.gene-repeat_region.vs.JCVI-Mt4.0-gene.kgb.synonymy.txt'));

%% Extract data for model genes

% Extract 3.5 to 4.0 data for model genes
conversionTableReduced = {};
for n = 1:length(model.genes)
    pos = strmatch(model.genes{n}, conversionTable(:,1), 'exact');
    for m = 1:length(pos)
        conversionTableReduced = vertcat(conversionTableReduced, conversionTable(pos(m),:));
    end
end
conversionTableReduced = conversionTableReduced(:,1:2);
for n = 1:length(conversionTableReduced)
    parts = splitstr(conversionTableReduced{n,2},'\.');
    conversionTableReduced{n,2} = parts{1,1};
end

% Extract 4.0 to 5.0 data for model genes
conversionTable2Reduced = {};
for n = 1:length(conversionTableReduced)
    pos = strmatch(conversionTableReduced{n,2}, conversionTable2(:,7), 'exact');
    for m = 1:length(pos)
        conversionTable2Reduced = vertcat(conversionTable2Reduced, conversionTable2(pos(m),:));
    end
end
conversionTable2Reduced_temp = conversionTable2Reduced(:,1);
conversionTable2Reduced = horzcat(conversionTable2Reduced_temp, conversionTable2Reduced(:,7));

%% Prepare final conversion table

% Prepare 3.5 to 5.0 conversion table
finalConversions = {};
for n = 1:length(conversionTableReduced)
    pos = strmatch(conversionTableReduced{n,2}, conversionTable2Reduced(:,2), 'exact');
    for m = 1:length(pos)
        finalConversions{end+1,1} = conversionTableReduced{n,1};
        finalConversions{end,2} = conversionTableReduced{n,2};
        finalConversions{end,3} = conversionTable2Reduced{pos(m),1};
    end
end
finalConversions = sortrows(finalConversions, [1 3]);
finalConversions_temp = {};
for n = 1:length(finalConversions)
    if n == 1
        finalConversions_temp(1,:) = finalConversions(1,:);
    else
        if strmatch(finalConversions{n,1}, finalConversions{n-1,1}, 'exact')
            if strmatch(finalConversions{n,3}, finalConversions{n-1,3}, 'exact')
            else
                finalConversions_temp(end+1,:) = finalConversions(n,:);
            end
        else
            finalConversions_temp(end+1,:) = finalConversions(n,:);
        end
    end
end
finalConversions = finalConversions_temp;
save('Data/finalConversions.mat', 'finalConversions');

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
pos = findRxnIDs(model, 'RXN-9944_H');
model.grRules{pos} = 'Medtr6g008180.1';
model.genes{end+1} = 'Medtr6g008180.1';
pos = findRxnIDs(model, 'RXN-7674_H');
model.grRules{pos} = 'Medtr7g069400.1';
model.genes{end+1} = 'Medtr7g069400.1';
pos = findRxnIDs(model, 'PLASTOQUINOL--PLASTOCYANIN-REDUCTASE-RXN_H');
model.grRules{pos} = 'Medtr2g021500.1';
model.genes{end+1} = 'Medtr2g021500.1';
pos = findRxnIDs(model, 'RXN0-6555_C');
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

