%% Load the data

load('constrainedNodule_sucrose.mat');
model = constrainedNodule;

% Import METANETX reaction database
metanetxRxns = table2cell(readtable('reac_xref.txt', 'ReadVariableNames', false));

% Import METANETX compound database
metanetxChem = table2cell(readtable('chem_xref.txt', 'ReadVariableNames', false));

%% Change reaction names

% Extract just the SEED reaction names
seedRxns = {};
x = 0;
for n = 1:length(metanetxRxns)
    if strmatch('seed', metanetxRxns{n,1})
        x = x + 1;
        seedRxns(x,:) = metanetxRxns(n,:);
    end
end

% Extract just the METACYC reaction names
metaRxns = {};
x = 0;
for n = 1:length(metanetxRxns)
    if strmatch('metacyc', metanetxRxns{n,1})
        x = x + 1;
        metaRxns(x,:) = metanetxRxns(n,:);
    end
end

% Extract model reaction names
reactionNames = cell(length(model.rxns), 3);
for n = 1:length(model.rxns)
    if strmatch('Leave_', model.rxns{n})
        reactionNames{n,1} = 'Leave_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Leave_', '');
    elseif strmatch('Root_', model.rxns{n})
        reactionNames{n,1} = 'Root_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'Root_', '');
    elseif strmatch('NoduleI_', model.rxns{n})
        reactionNames{n,1} = 'NoduleI_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'NoduleI_', '');
    elseif strmatch('NoduleIId_', model.rxns{n})
        reactionNames{n,1} = 'NoduleIId_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'NoduleIId_', '');
    elseif strmatch('NoduleIIp_', model.rxns{n})
        reactionNames{n,1} = 'NoduleIIp_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'NoduleIIp_', '');
    elseif strmatch('NoduleIZ_', model.rxns{n})
        reactionNames{n,1} = 'NoduleIZ_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'NoduleIZ_', '');
    elseif strmatch('NoduleIII_', model.rxns{n})
        reactionNames{n,1} = 'NoduleIII_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'NoduleIII_', '');
    elseif strmatch('BacteroidIId_', model.rxns{n})
        reactionNames{n,1} = 'BacteroidIId_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'BacteroidIId_', '');
    elseif strmatch('BacteroidIIp_', model.rxns{n})
        reactionNames{n,1} = 'BacteroidIIp_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'BacteroidIIp_', '');
    elseif strmatch('BacteroidIZ_', model.rxns{n})
        reactionNames{n,1} = 'BacteroidIZ_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'BacteroidIZ_', '');
    elseif strmatch('BacteroidIII_', model.rxns{n})
        reactionNames{n,1} = 'BacteroidIII_';
        reactionNames{n,2} = strrep(model.rxns{n}, 'BacteroidIII_', '');
    else
        reactionNames{n,2} = model.rxns{n};
    end
    splitString = strsplit(reactionNames{n,2}, '_');
    if length(splitString) == 2
        reactionNames{n,2} = splitString{1};
        reactionNames{n,3} = splitString{2};
    end
end

% Change reaction names, where possible
for n = 1:length(model.rxns)
    if strmatch('Bacteroid', reactionNames{n,1})
        rxn = ['seed:' reactionNames{n,2}];
        pos = strmatch(rxn, seedRxns(:,1), 'exact');
        if ~isempty(pos)
            reactionNames{n,2} = seedRxns{pos,2};
        end
    else
        rxn = ['metacyc:' reactionNames{n,2}];
        pos = strmatch(rxn, metaRxns(:,1), 'exact');
        if ~isempty(pos)
            reactionNames{n,2} = metaRxns{pos,2};
        end
    end
end
for n = 1:length(model.rxns)
    rxn = strcat(reactionNames{n,1}, reactionNames{n,2});
    if ~isempty(reactionNames{n,3})
        rxnTemp = strcat(rxn, '_', reactionNames{n,3});
        if strmatch('MNXR', reactionNames{n,2})
            if strmatch(rxnTemp, model.rxns, 'exact')
                model.rxns{n,1} = strcat(rxn, 'b', '_', reactionNames{n,3});
            else
                model.rxns{n,1} = strcat(rxn, '_', reactionNames{n,3});
            end
        else
            model.rxns{n,1} = strcat(rxn, '_', reactionNames{n,3});
        end
    else
        if strmatch('MNXR', reactionNames{n,2})
            if strmatch(rxn, model.rxns, 'exact')
                model.rxns{n,1} = strcat(rxn, 'b');
            else
                model.rxns{n,1} = rxn;
            end
        else
            model.rxns{n,1} = rxn;
        end
    end
end

%% Change metabolite names

% Extract just the SEED compound names
seedMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('seed', metanetxChem{n,1})
        x = x + 1;
        seedMets(x,:) = metanetxChem(n,:);
    end
end

% Extract just the METACYC compound names
metaMets = {};
x = 0;
for n = 1:length(metanetxChem)
    if strmatch('metacyc', metanetxChem{n,1})
        x = x + 1;
        metaMets(x,:) = metanetxChem(n,:);
    end
end

% Extract model compound names
cpdNames = cell(length(model.mets), 3);
for n = 1:length(model.mets)
    if strmatch('Leave_', model.mets{n})
        cpdNames{n,1} = 'Leave_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Leave_', '');
    elseif strmatch('Root_', model.mets{n})
        cpdNames{n,1} = 'Root_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Root_', '');
    elseif strmatch('Nodule_', model.mets{n})
        cpdNames{n,1} = 'Nodule_';
        cpdNames{n,2} = strrep(model.mets{n}, 'Nodule_', '');
    elseif strmatch('NoduleI_', model.mets{n})
        cpdNames{n,1} = 'NoduleI_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleI_', '');
    elseif strmatch('NoduleIId_', model.mets{n})
        cpdNames{n,1} = 'NoduleIId_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIId_', '');
    elseif strmatch('NoduleIIp_', model.mets{n})
        cpdNames{n,1} = 'NoduleIIp_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIIp_', '');
    elseif strmatch('NoduleIZ_', model.mets{n})
        cpdNames{n,1} = 'NoduleIZ_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIZ_', '');
    elseif strmatch('NoduleIII_', model.mets{n})
        cpdNames{n,1} = 'NoduleIII_';
        cpdNames{n,2} = strrep(model.mets{n}, 'NoduleIII_', '');
    elseif strmatch('BacteroidIId_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIId_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIId_', '');
    elseif strmatch('BacteroidIIp_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIIp_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIIp_', '');
    elseif strmatch('BacteroidIZ_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIZ_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIZ_', '');
    elseif strmatch('BacteroidIII_', model.mets{n})
        cpdNames{n,1} = 'BacteroidIII_';
        cpdNames{n,2} = strrep(model.mets{n}, 'BacteroidIII_', '');
    else
        cpdNames{n,2} = model.mets{n};
    end
    splitString = strsplit(cpdNames{n,2}, '[');
    if length(splitString) == 2
        cpdNames{n,2} = splitString{1};
        cpdNames{n,3} = splitString{2};
    end
end

% Change compound names, where possible
for n = 1:length(model.mets)
    if strmatch('Bacteroid', cpdNames{n,1})
        cpd = ['seed:' cpdNames{n,2}];
        pos = strmatch(cpd, seedMets(:,1), 'exact');
        if ~isempty(pos)
            cpdNames{n,2} = seedMets{pos,2};
        end
    else
        cpd = ['metacyc:' cpdNames{n,2}];
        pos = strmatch(cpd, metaMets(:,1), 'exact');
        if ~isempty(pos)
            cpdNames{n,2} = metaMets{pos,2};
        end
    end
end
for n = 1:length(model.mets)
    cpd = strcat(cpdNames{n,1}, cpdNames{n,2});
    if ~isempty(cpdNames{n,3})
        cpdTemp = strcat(cpd, '_', reactionNames{n,3});
        if strmatch(cpdTemp, model.mets, 'exact')
            model.mets{n,1} = strcat(cpd, 'b', '[', cpdNames{n,3});
        else
            model.mets{n,1} = strcat(cpd, '[', cpdNames{n,3});
        end        
    else
        if strmatch(cpdTemp, model.mets, 'exact')
            model.mets{n,1} = strcat(cpd, 'b');
        else
            model.mets{n,1} = cpd;
        end        
    end
end

%% Remove duplicates

% Get rid of duplicates
[A,B,C] = checkDuplicateRxn(model, 'S');
C = model.rxns(C);
model = tncore_remove_reactions(model, C);

% Remove any unused genes
model = tncore_remove(model);
optimizeCbModel(model)

%% Finish

finalNodulatedPlant = model;
save('allWorkspace.mat');
save('finalNodulatedPlant_sucrose.mat', 'finalNodulatedPlant');
clear

