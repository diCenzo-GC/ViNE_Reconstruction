%% Find reactions to connect nodule to root metabolism and the soil

% Identify reactions for import and export between root and soil
importLogical = cell(length(plantModel.rxns),1);
exportLogical = cell(length(plantModel.rxns),1);
for n = 1:length(plantModel.rxns)
    if strmatch('Root_TEC_', plantModel.rxns{n})
        importLogical{n,1} = 1;
    elseif strmatch('Root_TEH_', plantModel.rxns{n})
        importLogical{n,1} = 1;
    else
        importLogical{n,1} = 0;
    end
    if strmatch('Root_TCE_', plantModel.rxns{n})
        exportLogical{n,1} = 1;
    elseif strmatch('Root_THE_', plantModel.rxns{n})
        exportLogical{n,1} = 1;
    elseif strmatch('Root_TGE_', plantModel.rxns{n})
        exportLogical{n,1} = 1;
    else
        exportLogical{n,1} = 0;
    end
end
importLogical = logical(cell2mat(importLogical));
importRxns = plantModel.rxns(importLogical);
importIDs = findRxnIDs(plantModel, importRxns);
exportLogical = logical(cell2mat(exportLogical));
exportRxns = plantModel.rxns(exportLogical);
exportIDs = findRxnIDs(plantModel, exportRxns);

% Identify reactions for exchange between root and shoot
TRS_reactions = plantModel.rxns(strmatch('TRS_', plantModel.rxns));
TSR_reactions = plantModel.rxns(strmatch('TSR_', plantModel.rxns));
transferRxns = vertcat(TRS_reactions, TSR_reactions);
transferIDs = findRxnIDs(plantModel, transferRxns);

%% Add export reactions to the nodule

% Set new model name
nodulatedPlantDisconnected_export = nodulatedPlantDisconnected;

% Add export reactions for zone I
for n = 1:length(exportRxns)
    rxnAbr = strrep(exportRxns{n}, 'Root_', 'NoduleI_');
    rxnFormula = printRxnFormula(plantModel, exportRxns{n});
    rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleI_');
    nodulatedPlantDisconnected_export = addReaction(nodulatedPlantDisconnected_export, ...
        rxnAbr, cell2mat(rxnFormula), [], 0, 0, 1000, 0);
end

% Add export reactions for zone IId
for n = 1:length(exportRxns)
    rxnAbr = strrep(exportRxns{n}, 'Root_', 'NoduleIId_');
    rxnFormula = printRxnFormula(plantModel, exportRxns{n});
    rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIId_');
    nodulatedPlantDisconnected_export = addReaction(nodulatedPlantDisconnected_export, ...
        rxnAbr, cell2mat(rxnFormula), [], 0, 0, 1000, 0);
end

% Add export reactions for zone IIp
for n = 1:length(exportRxns)
    rxnAbr = strrep(exportRxns{n}, 'Root_', 'NoduleIIp_');
    rxnFormula = printRxnFormula(plantModel, exportRxns{n});
    rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIIp_');
    nodulatedPlantDisconnected_export = addReaction(nodulatedPlantDisconnected_export, ...
        rxnAbr, cell2mat(rxnFormula), [], 0, 0, 1000, 0);
end

% Add export reactions for zone IZ
for n = 1:length(exportRxns)
    rxnAbr = strrep(exportRxns{n}, 'Root_', 'NoduleIZ_');
    rxnFormula = printRxnFormula(plantModel, exportRxns{n});
    rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIZ_');
    nodulatedPlantDisconnected_export = addReaction(nodulatedPlantDisconnected_export, ...
        rxnAbr, cell2mat(rxnFormula), [], 0, 0, 1000, 0);
end

% Add export reactions for zone III
for n = 1:length(exportRxns)
    rxnAbr = strrep(exportRxns{n}, 'Root_', 'NoduleIII_');
    rxnFormula = printRxnFormula(plantModel, exportRxns{n});
    rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIII_');
    nodulatedPlantDisconnected_export = addReaction(nodulatedPlantDisconnected_export, ...
        rxnAbr, cell2mat(rxnFormula), [], 0, 0, 1000, 0);
end

%% Allow compounds imported by root to be imported by the nodule

% Set new model name
nodulatedPlantDisconnected_import = nodulatedPlantDisconnected_export;

% Add import to general nodule from soil
for n = 1:length(importRxns)
    rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TRN_');
    rxnFormula = printRxnFormula(plantModel, importRxns{n});
    [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
    if strmatch('Root_TEC_Pi', importRxns{n}, 'exact')
        metToAdd = {'Nodule_Pi[C]'};
    else
        metToAdd = setdiff(metaboliteList, ...
            {'Root_Pi[C]', 'Root_ATP[C]', 'Root_ADP[C]', 'Root_PROTON[C]'});
        metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
    end
    nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
        rxnAbr, metToAdd, [1], 0, 0, 1000, 0);
end

% Add transfer from general nodule to zone I
for n = 1:length(importRxns)
    if strmatch('Root_TEC_Pi', importRxns{n}, 'exact')
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNI_');
        rxnFormula = 'Nodule_Pi[C] + 0.25 NoduleI_ATP[C] -> 0.25 NoduleI_ADP[C] + 0.25 NoduleI_PROTON[C] + 1.25 NoduleI_Pi[C]';
    else
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNI_');
        rxnFormula = printRxnFormula(nodulatedPlantDisconnected_import, importRxns{n});
        [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
        metToAdd = setdiff(metaboliteList, {'Root_Pi[C]', 'Root_ATP[C]', 'Root_ADP[C]', 'Root_PROTON[C]'});
        metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
        rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleI_');
        rxnFormula = strcat(cell2mat(metToAdd), ' +_____', cell2mat(rxnFormula));
        rxnFormula = strrep(rxnFormula, '_____', ' ');
    end
    nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
        rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
end

% Add transfer from general nodule to zone IId
for n = 1:length(importRxns)
    if strmatch('Root_TEC_Pi', importRxns{n}, 'exact')
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIId_');
        rxnFormula = 'Nodule_Pi[C] + 0.25 NoduleIId_ATP[C] -> 0.25 NoduleIId_ADP[C] + 0.25 NoduleIId_PROTON[C] + 1.25 NoduleIId_Pi[C]';
    else
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIId_');
        rxnFormula = printRxnFormula(nodulatedPlantDisconnected_import, importRxns{n});
        [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
        metToAdd = setdiff(metaboliteList, {'Root_Pi[C]', 'Root_ATP[C]', 'Root_ADP[C]', 'Root_PROTON[C]'});
        metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
        rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIId_');
        rxnFormula = strcat(cell2mat(metToAdd), ' +_____', cell2mat(rxnFormula));
        rxnFormula = strrep(rxnFormula, '_____', ' ');
    end
    nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
        rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
end

% Add transfer from general nodule to zone IIp
for n = 1:length(importRxns)
    if strmatch('Root_TEC_Pi', importRxns{n}, 'exact')
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIIp_');
        rxnFormula = 'Nodule_Pi[C] + 0.25 NoduleIIp_ATP[C] -> 0.25 NoduleIIp_ADP[C] + 0.25 NoduleIIp_PROTON[C] + 1.25 NoduleIIp_Pi[C]';
    else
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIIp_');
        rxnFormula = printRxnFormula(nodulatedPlantDisconnected_import, importRxns{n});
        [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
        metToAdd = setdiff(metaboliteList, {'Root_Pi[C]', 'Root_ATP[C]', 'Root_ADP[C]', 'Root_PROTON[C]'});
        metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
        rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIIp_');
        rxnFormula = strcat(cell2mat(metToAdd), ' +_____', cell2mat(rxnFormula));
        rxnFormula = strrep(rxnFormula, '_____', ' ');
    end
    nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
        rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
end

% Add transfer from general nodule to zone IZ
for n = 1:length(importRxns)
    if strmatch('Root_TEC_Pi', importRxns{n}, 'exact')
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIZ_');
        rxnFormula = 'Nodule_Pi[C] + 0.25 NoduleIZ_ATP[C] -> 0.25 NoduleIZ_ADP[C] + 0.25 NoduleIZ_PROTON[C] + 1.25 NoduleIZ_Pi[C]';
    else
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIZ_');
        rxnFormula = printRxnFormula(nodulatedPlantDisconnected_import, importRxns{n});
        [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
        metToAdd = setdiff(metaboliteList, {'Root_Pi[C]', 'Root_ATP[C]', 'Root_ADP[C]', 'Root_PROTON[C]'});
        metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
        rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIZ_');
        rxnFormula = strcat(cell2mat(metToAdd), ' +_____', cell2mat(rxnFormula));
        rxnFormula = strrep(rxnFormula, '_____', ' ');
    end
    nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
        rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
end

% Add transfer from general nodule to zone III
for n = 1:length(importRxns)
    if strmatch('Root_TEC_Pi', importRxns{n}, 'exact')
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIII_');
        rxnFormula = 'Nodule_Pi[C] + 0.25 NoduleIII_ATP[C] -> 0.25 NoduleIII_ADP[C] + 0.25 NoduleIII_PROTON[C] + 1.25 NoduleIII_Pi[C]';
    else
        rxnAbr = strrep(importRxns{n}, 'Root_TEC_', 'TNIII_');
        rxnFormula = printRxnFormula(nodulatedPlantDisconnected_import, importRxns{n});
        [metaboliteList] = parseRxnFormula(rxnFormula{1,1});
        metToAdd = setdiff(metaboliteList, {'Root_Pi[C]', 'Root_ATP[C]', 'Root_ADP[C]', 'Root_PROTON[C]'});
        metToAdd = strrep(metToAdd, 'Root_', 'Nodule_');
        rxnFormula = strrep(rxnFormula, 'Root_', 'NoduleIII_');
        rxnFormula = strcat(cell2mat(metToAdd), ' +_____', cell2mat(rxnFormula));
        rxnFormula = strrep(rxnFormula, '_____', ' ');
    end
    nodulatedPlantDisconnected_import = addReaction(nodulatedPlantDisconnected_import, ...
        rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
end

%% Allow compounds transferred between root and shoot to be transfered to the nodule

% Set new model name
nodulatedPlantDisconnected_transfer = nodulatedPlantDisconnected_import;

% Add transfer to general nodule from root
for n = 1:length(transferRxns)
    if strmatch('TSR_PROTON', transferRxns{n})
        rxnAbr = 'TRN_PROTON';
        rxnFormula = 'Root_PROTON[C] -> Nodule_PROTON[C]';
        if strmatch('TRN_PROTON', nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = 'TRN_2_PROTON';
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    elseif strmatch('TRS_PROTON', transferRxns{n})
        rxnAbr = 'TRN_PROTON';
        rxnFormula = 'Root_PROTON[C] -> Nodule_PROTON[C]';
        if strmatch('TRN_PROTON', nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = 'TRN_2_PROTON';
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    else
        unwantedMets = {'Leave_ATP[C]'; 'Leave_Pi[C]'; 'Leave_PROTON[C]'; 'Leave_ADP[C]'; ...
            'Root_ATP[C]'; 'Root_Pi[C]'; 'Root_PROTON[C]'; 'Root_ADP[C]'};
        transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
            logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
            nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
        transferredMet = transferredMet(strmatch('Root_', transferredMet));
        rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TRN_');
        rxnAbr = strrep(rxnAbr, 'TRS_', 'TRN_');
        receivedMet = strrep(transferredMet, 'Root_', 'Nodule_');
        rxnFormula = [cell2mat(transferredMet) ' -> ' cell2mat(receivedMet)];
        if strmatch(rxnAbr, nodulatedPlantDisconnected_import.rxns, 'exact')
            rxnAbr = strrep(rxnAbr, 'TRN_', 'TRN_2_');
        end
        nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
            rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
    end
end

% Add transfer from general nodule to zone I
for n = 1:length(transferRxns)
    rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TNI_');
    rxnAbr = strrep(rxnAbr, 'TRS_', 'TNI_');
    if strmatch(rxnAbr, nodulatedPlantDisconnected_transfer.rxns, 'exact')
    else
        if strmatch('TSR_PROTON', transferRxns{n})
            rxnFormula = '0.25 NoduleI_ATP[C] + Nodule_PROTON[C] -> 0.25 NoduleI_ADP[C] + 0.25 NoduleI_Pi[C] + 1.25 NoduleI_PROTON[C]';
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        else
            unwantedMets = {'Leave_ATP[C]'; 'Leave_Pi[C]'; 'Leave_PROTON[C]'; 'Leave_ADP[C]'; ...
                'Root_ATP[C]'; 'Root_Pi[C]'; 'Root_PROTON[C]'; 'Root_ADP[C]'};
            transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
                logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
                nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
            transferredMet = transferredMet(strmatch('Root_', transferredMet));
            donorMet = strrep(transferredMet, 'Root_', 'Nodule_');
            receivedMet = strrep(transferredMet, 'Root_', 'NoduleI_');
            rxnFormula = ['0.25 NoduleI_ATP[C] +_____' cell2mat(donorMet) ...
                ' -> 0.25 NoduleI_ADP[C] + 0.25 NoduleI_Pi[C] + 0.25 NoduleI_PROTON[C] +_____' ...
                cell2mat(receivedMet)];
            rxnFormula = strrep(rxnFormula, '+_____', '+ ');
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        end
    end
end

% Add transfer from general nodule to zone IId
for n = 1:length(transferRxns)
    rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TNIId_');
    rxnAbr = strrep(rxnAbr, 'TRS_', 'TNIId_');
    if strmatch(rxnAbr, nodulatedPlantDisconnected_transfer.rxns, 'exact')
    else
        if strmatch('TSR_PROTON', transferRxns{n})
            rxnFormula = '0.25 NoduleIId_ATP[C] + Nodule_PROTON[C] -> 0.25 NoduleIId_ADP[C] + 0.25 NoduleIId_Pi[C] + 1.25 NoduleIId_PROTON[C]';
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        else
            unwantedMets = {'Leave_ATP[C]'; 'Leave_Pi[C]'; 'Leave_PROTON[C]'; 'Leave_ADP[C]'; ...
                'Root_ATP[C]'; 'Root_Pi[C]'; 'Root_PROTON[C]'; 'Root_ADP[C]'};
            transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
                logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
                nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
            transferredMet = transferredMet(strmatch('Root_', transferredMet));
            donorMet = strrep(transferredMet, 'Root_', 'Nodule_');
            receivedMet = strrep(transferredMet, 'Root_', 'NoduleIId_');
            rxnFormula = ['0.25 NoduleIId_ATP[C] +_____' cell2mat(donorMet) ...
                ' -> 0.25 NoduleIId_ADP[C] + 0.25 NoduleIId_Pi[C] + 0.25 NoduleIId_PROTON[C] +_____' ...
                cell2mat(receivedMet)];
            rxnFormula = strrep(rxnFormula, '+_____', '+ ');
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        end
    end
end

% Add transfer from general nodule to zone IIp
for n = 1:length(transferRxns)
    rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TNIIp_');
    rxnAbr = strrep(rxnAbr, 'TRS_', 'TNIIp_');
    if strmatch(rxnAbr, nodulatedPlantDisconnected_transfer.rxns, 'exact')
    else
        if strmatch('TSR_PROTON', transferRxns{n})
            rxnFormula = '0.25 NoduleIIp_ATP[C] + Nodule_PROTON[C] -> 0.25 NoduleIIp_ADP[C] + 0.25 NoduleIIp_Pi[C] + 1.25 NoduleIIp_PROTON[C]';
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        else
            unwantedMets = {'Leave_ATP[C]'; 'Leave_Pi[C]'; 'Leave_PROTON[C]'; 'Leave_ADP[C]'; ...
                'Root_ATP[C]'; 'Root_Pi[C]'; 'Root_PROTON[C]'; 'Root_ADP[C]'};
            transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
                logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
                nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
            transferredMet = transferredMet(strmatch('Root_', transferredMet));
            donorMet = strrep(transferredMet, 'Root_', 'Nodule_');
            receivedMet = strrep(transferredMet, 'Root_', 'NoduleIIp_');
            rxnFormula = ['0.25 NoduleIIp_ATP[C] +_____' cell2mat(donorMet) ...
                ' -> 0.25 NoduleIIp_ADP[C] + 0.25 NoduleIIp_Pi[C] + 0.25 NoduleIIp_PROTON[C] +_____' ...
                cell2mat(receivedMet)];
            rxnFormula = strrep(rxnFormula, '+_____', '+ ');
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        end
    end
end

% Add transfer from general nodule to zone IZ
for n = 1:length(transferRxns)
    rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TNIZ_');
    rxnAbr = strrep(rxnAbr, 'TRS_', 'TNIZ_');
    if strmatch(rxnAbr, nodulatedPlantDisconnected_transfer.rxns, 'exact')
    else
        if strmatch('TSR_PROTON', transferRxns{n})
            rxnFormula = '0.25 NoduleIZ_ATP[C] + Nodule_PROTON[C] -> 0.25 NoduleIZ_ADP[C] + 0.25 NoduleIZ_Pi[C] + 1.25 NoduleIZ_PROTON[C]';
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        else
            unwantedMets = {'Leave_ATP[C]'; 'Leave_Pi[C]'; 'Leave_PROTON[C]'; 'Leave_ADP[C]'; ...
                'Root_ATP[C]'; 'Root_Pi[C]'; 'Root_PROTON[C]'; 'Root_ADP[C]'};
            transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
                logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
                nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
            transferredMet = transferredMet(strmatch('Root_', transferredMet));
            donorMet = strrep(transferredMet, 'Root_', 'Nodule_');
            receivedMet = strrep(transferredMet, 'Root_', 'NoduleIZ_');
            rxnFormula = ['0.25 NoduleIZ_ATP[C] +_____' cell2mat(donorMet) ...
                ' -> 0.25 NoduleIZ_ADP[C] + 0.25 NoduleIZ_Pi[C] + 0.25 NoduleIZ_PROTON[C] +_____' ...
                cell2mat(receivedMet)];
            rxnFormula = strrep(rxnFormula, '+_____', '+ ');
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        end
    end
end

% Add transfer from general nodule to zone III
for n = 1:length(transferRxns)
    rxnAbr = strrep(transferRxns{n}, 'TSR_', 'TNIII_');
    rxnAbr = strrep(rxnAbr, 'TRS_', 'TNIII_');
    if strmatch(rxnAbr, nodulatedPlantDisconnected_transfer.rxns, 'exact')
    else
        if strmatch('TSR_PROTON', transferRxns{n})
            rxnFormula = '0.25 NoduleIII_ATP[C] + Nodule_PROTON[C] -> 0.25 NoduleIII_ADP[C] + 0.25 NoduleIII_Pi[C] + 1.25 NoduleIII_PROTON[C]';
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        else
            unwantedMets = {'Leave_ATP[C]'; 'Leave_Pi[C]'; 'Leave_PROTON[C]'; 'Leave_ADP[C]'; ...
                'Root_ATP[C]'; 'Root_Pi[C]'; 'Root_PROTON[C]'; 'Root_ADP[C]'};
            transferredMet = setdiff(nodulatedPlantDisconnected_transfer.mets(...
                logical(nodulatedPlantDisconnected_transfer.S(:,findRxnIDs(...
                nodulatedPlantDisconnected_transfer, transferRxns{n})))), unwantedMets);
            transferredMet = transferredMet(strmatch('Root_', transferredMet));
            donorMet = strrep(transferredMet, 'Root_', 'Nodule_');
            receivedMet = strrep(transferredMet, 'Root_', 'NoduleIII_');
            rxnFormula = ['0.25 NoduleIII_ATP[C] +_____' cell2mat(donorMet) ...
                ' -> 0.25 NoduleIII_ADP[C] + 0.25 NoduleIII_Pi[C] + 0.25 NoduleIII_PROTON[C] +_____' ...
                cell2mat(receivedMet)];
            rxnFormula = strrep(rxnFormula, '+_____', '+ ');
            nodulatedPlantDisconnected_transfer = addReaction(nodulatedPlantDisconnected_transfer, ...
                rxnAbr, rxnFormula, [], 0, 0, 1000, 0);
        end
    end
end

% Rename the model
nodulatedPlantConnected = nodulatedPlantDisconnected_transfer;

%% Deal with carbon dioxide

% Add carbon dioxide entry to nodule
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TRN_CARBON-DIOXIDE', ...
    {'Nodule_CARBON-DIOXIDE[C]'}, [1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNI_CARBON-DIOXIDE', ...
    {'Nodule_CARBON-DIOXIDE[C]', 'NoduleI_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIId_CARBON-DIOXIDE', ...
    {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIId_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIIp_CARBON-DIOXIDE', ...
    {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIIp_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIZ_CARBON-DIOXIDE', ...
    {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIZ_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIII_CARBON-DIOXIDE', ...
    {'Nodule_CARBON-DIOXIDE[C]', 'NoduleIII_CARBON-DIOXIDE[C]'}, [-1 1], 0, 0, 1000, 0);

%% Transfer asparagine and glutamine from zone III to the root

% Add asparagine transfer
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIIIR_ASN', ...
    {'NoduleIII_ASN[C]', 'NoduleIII_ATP[C]', 'NoduleIII_ADP[C]', ...
    'NoduleIII_Pi[C]', 'NoduleIII_PROTON[C]', 'Root_ASN[C]'}, ...
    [-1 -0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);

% Add glutamine transfer
nodulatedPlantConnected = addReaction(nodulatedPlantConnected, 'TNIIIR_GLN', ...
    {'NoduleIII_GLN[C]', 'NoduleIII_ATP[C]', 'NoduleIII_ADP[C]', ...
    'NoduleIII_Pi[C]', 'NoduleIII_PROTON[C]', 'Root_GLN[C]'}, ...
    [-1 -0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);

%% Remove reaction adding ammonium to the nodule from the soil

nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'TRN_AMMONIUM');
nodulatedPlantConnected = tncore_remove_reactions(nodulatedPlantConnected, 'TRN_NITRATE');

