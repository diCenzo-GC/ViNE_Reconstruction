%% Zone-specific gene deletion

% Get nodule and bacteroid genes
noduleGenes = model.genes(strmatch('Nodule', model.genes));
bacteroidGenes = model.genes(strmatch('Bacteroid', model.genes));

% Perform gene deletion analyses
grRatioNodule = singleGeneDeletion(model, 'FBA', noduleGenes);
grRatioBacteroid = singleGeneDeletion(model, 'FBA', bacteroidGenes);
grRatioNodule(isnan(grRatioNodule)) = 0;
grRatioBacteroid(isnan(grRatioBacteroid)) = 0;

% Get list of zone-independent nodule genes
noduleGenesOverall = cell(length(noduleGenes), 1);
for n = 1:length(noduleGenes)
    temp = strsplit(noduleGenes{n}, '_');
    noduleGenesOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            noduleGenesOverall{n} = strcat(noduleGenesOverall{n}, '_', temp{m});
        end
    end
end
noduleGenesOverall = unique(noduleGenesOverall);

% Get list of zone-independent bacteroid genes
bacteroidGenesOverall = cell(length(bacteroidGenes), 1);
for n = 1:length(bacteroidGenes)
    temp = strsplit(bacteroidGenes{n}, '_');
    bacteroidGenesOverall{n} = temp{1,2};
    if length(bacteroidGenesOverall) > 2
        for m = 3:length(temp)
            bacteroidGenesOverall{n} = strcat(bacteroidGenesOverall{n}, '_', temp{m});
        end
    end
end
bacteroidGenesOverall = unique(bacteroidGenesOverall);

% Store the nodule data as a matrix
temp = cell(length(noduleGenesOverall), 5);
for n = 1:length(noduleGenes)
    geneOrig = strsplit(noduleGenes{n}, '_');
    temporary = geneOrig{1,2};
    if length(geneOrig) > 2
        for m = 3:length(geneOrig)
            temporary = strcat(temporary, '_', geneOrig{m});
        end
    end
    pos = strmatch(temporary, noduleGenesOverall, 'exact');
    if strmatch('NoduleI', geneOrig{1,1}, 'exact')
        temp{pos,1} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIId', geneOrig{1,1}, 'exact')
        temp{pos,2} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIIp', geneOrig{1,1}, 'exact')
        temp{pos,3} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIZ', geneOrig{1,1}, 'exact')
        temp{pos,4} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIII', geneOrig{1,1}, 'exact')
        temp{pos,5} = round(grRatioNodule(n), 2);
    end
end
temp = horzcat(noduleGenesOverall, temp);
output_zoneSpecificGeneDeletion_nodule = ...
    vertcat({'Gene','Zone_I','Zone_IId','Zone_IIp','Zone_IZ','Zone_III'}, temp);
isEmpty = cellfun(@isempty, output_zoneSpecificGeneDeletion_nodule);
output_zoneSpecificGeneDeletion_nodule(isEmpty) = {1};

% Store the bacteroid data as a matrix
temp = cell(length(bacteroidGenesOverall), 4);
for n = 1:length(bacteroidGenes)
    geneOrig = strsplit(bacteroidGenes{n}, '_');
    temporary = geneOrig{1,2};
    if length(geneOrig) > 2
        for m = 3:length(geneOrig)
            temporary = strcat(temporary, '_', geneOrig{m});
        end
    end
    pos = strmatch(temporary, bacteroidGenesOverall, 'exact');
    if strmatch('BacteroidIId', geneOrig{1,1}, 'exact')
        temp{pos,1} = round(grRatioBacteroid(n), 2);
    elseif strmatch('BacteroidIIp', geneOrig{1,1}, 'exact')
        temp{pos,2} = round(grRatioBacteroid(n), 2);
    elseif strmatch('BacteroidIZ', geneOrig{1,1}, 'exact')
        temp{pos,3} = round(grRatioBacteroid(n), 2);
    elseif strmatch('BacteroidIII', geneOrig{1,1}, 'exact')
        temp{pos,4} = round(grRatioBacteroid(n), 2);
    end
end
temp = horzcat(bacteroidGenesOverall, temp);
output_zoneSpecificGeneDeletion_bacteroid = ...
    vertcat({'Gene','Zone_IId','Zone_IIp','Zone_IZ','Zone_III'}, temp);
isEmpty = cellfun(@isempty, output_zoneSpecificGeneDeletion_bacteroid);
output_zoneSpecificGeneDeletion_bacteroid(isEmpty) = {1};