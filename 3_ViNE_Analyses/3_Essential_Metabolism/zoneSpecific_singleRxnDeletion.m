%% Zone-specific reaction deletion

% Get nodule and bacteroid reactions
noduleRxns = model.rxns(strmatch('Nodule', model.rxns));
noduleRxns = setdiff(noduleRxns, {'NoduleBiomass'});
bacteroidRxns = model.rxns(strmatch('Bacteroid', model.rxns));

% Perform reaction deletion analyses
grRatioNodule = singleRxnDeletion(model, 'FBA', noduleRxns);
grRatioBacteroid = singleRxnDeletion(model, 'FBA', bacteroidRxns);
grRatioNodule(isnan(grRatioNodule)) = 0;
grRatioBacteroid(isnan(grRatioBacteroid)) = 0;

% Get list of zone-independent nodule reactions
noduleRxnsOverall = cell(length(noduleRxns), 1);
for n = 1:length(noduleRxns)
    temp = strsplit(noduleRxns{n}, '_');
    noduleRxnsOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            noduleRxnsOverall{n} = strcat(noduleRxnsOverall{n}, '_', temp{m});
        end
    end
end
noduleRxnsOverall = unique(noduleRxnsOverall);

% Get list of zone-independent bacteroid reactions
bacteroidRxnsOverall = cell(length(bacteroidRxns), 1);
for n = 1:length(bacteroidRxns)
    temp = strsplit(bacteroidRxns{n}, '_');
    bacteroidRxnsOverall{n} = temp{1,2};
    if length(temp) > 2
        for m = 3:length(temp)
            bacteroidRxnsOverall{n} = strcat(bacteroidRxnsOverall{n}, '_', temp{m});
        end
    end
end
bacteroidRxnsOverall = unique(bacteroidRxnsOverall);

% Store the nodule data as a matrix
temp = cell(length(noduleRxnsOverall), 5);
for n = 1:length(noduleRxns)
    rxnOrig = strsplit(noduleRxns{n}, '_');
    temporary = rxnOrig{1,2};
    if length(rxnOrig) > 2
        for m = 3:length(rxnOrig)
            temporary = strcat(temporary, '_', rxnOrig{m});
        end
    end
    pos = strmatch(temporary, noduleRxnsOverall, 'exact');
    if strmatch('NoduleI', rxnOrig{1,1}, 'exact')
        temp{pos,1} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIId', rxnOrig{1,1}, 'exact')
        temp{pos,2} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIIp', rxnOrig{1,1}, 'exact')
        temp{pos,3} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIZ', rxnOrig{1,1}, 'exact')
        temp{pos,4} = round(grRatioNodule(n), 2);
    elseif strmatch('NoduleIII', rxnOrig{1,1}, 'exact')
        temp{pos,5} = round(grRatioNodule(n), 2);
    end
end
temp = horzcat(noduleRxnsOverall, temp);
output_zoneSpecificRxnDeletion_nodule = ...
    vertcat({'Reaction','Zone_I','Zone_IId','Zone_IIp','Zone_IZ','Zone_III'}, temp);
isEmpty = cellfun(@isempty, output_zoneSpecificRxnDeletion_nodule);
output_zoneSpecificRxnDeletion_nodule(isEmpty) = {1};

% Store the bacteroid data as a matrix
temp = cell(length(bacteroidRxnsOverall), 4);
for n = 1:length(bacteroidRxns)
    rxnOrig = strsplit(bacteroidRxns{n}, '_');
    temporary = rxnOrig{1,2};
    if length(rxnOrig) > 2
        for m = 3:length(rxnOrig)
            temporary = strcat(temporary, '_', rxnOrig{m});
        end
    end
    pos = strmatch(temporary, bacteroidRxnsOverall, 'exact');
    if strmatch('BacteroidIId', rxnOrig{1,1}, 'exact')
        temp{pos,1} = round(grRatioBacteroid(n), 2);
    elseif strmatch('BacteroidIIp', rxnOrig{1,1}, 'exact')
        temp{pos,2} = round(grRatioBacteroid(n), 2);
    elseif strmatch('BacteroidIZ', rxnOrig{1,1}, 'exact')
        temp{pos,3} = round(grRatioBacteroid(n), 2);
    elseif strmatch('BacteroidIII', rxnOrig{1,1}, 'exact')
        temp{pos,4} = round(grRatioBacteroid(n), 2);
    end
end
temp = horzcat(bacteroidRxnsOverall, temp);
output_zoneSpecificRxnDeletion_bacteroid = ...
    vertcat({'Reaction','Zone_IId','Zone_IIp','Zone_IZ','Zone_III'}, temp);
isEmpty = cellfun(@isempty, output_zoneSpecificRxnDeletion_bacteroid);
output_zoneSpecificRxnDeletion_bacteroid(isEmpty) = {1};