%% Rename reactions, metabolites, and genes in the nodule zone I

% Add 'Nodule_' to the start of all plant nodule reactions
for n = 1:length(medicagoModel.rxns)
    rxnID = findRxnIDs(zoneI, medicagoModel.rxns{n});
    zoneI.rxns{rxnID} = strcat('NoduleI_', zoneI.rxns{rxnID});
end

% Add 'Nodule_' to the start of all plant nodule metabolites
for n = 1:length(medicagoModel.mets)
    metID = findMetIDs(zoneI, medicagoModel.mets{n});
    zoneI.mets{metID} = strcat('NoduleI_', zoneI.mets{metID});
end

% Add 'Nodule_' to the start of all plant nodule genes
for n = 1:length(medicagoModel.genes)
    geneID = findGeneIDs(zoneI, medicagoModel.genes{n});
    zoneI.genes{geneID} = strcat('NoduleI_', zoneI.genes{geneID});
end

% Update the grRules based on new gene names
for n = 1:length(medicagoModel.grRules)
    for m = 1:length(medicagoModel.genes)
        zoneI.grRules{n} = strrep(zoneI.grRules{n}, medicagoModel.genes{m}, zoneI.genes{m});
    end
end

%% Rename reactions, metabolites, and genes in the nodule zone IId

% Add 'Nodule_' to the start of all plant nodule reactions
for n = 1:length(medicagoModel.rxns)
    rxnID = findRxnIDs(zoneIId, medicagoModel.rxns{n});
    zoneIId.rxns{rxnID} = strcat('NoduleIId_', zoneIId.rxns{rxnID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.rxns)
    rxnID = findRxnIDs(zoneIId, melilotiModel.rxns{n});
    zoneIId.rxns{rxnID} = strcat('BacteroidIId_', zoneIId.rxns{rxnID});
end

% Add 'Nodule_ to the start of the EXTC reactions
zoneIId.rxns = strrep(zoneIId.rxns, 'EXCT_', 'NoduleIId_EXCT_');

% Add 'Nodule_' to the start of all plant nodule metabolites
for n = 1:length(medicagoModel.mets)
    metID = findMetIDs(zoneIId, medicagoModel.mets{n});
    zoneIId.mets{metID} = strcat('NoduleIId_', zoneIId.mets{metID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule metabolites
for n = 1:length(melilotiModel.mets)
    metID = findMetIDs(zoneIId, melilotiModel.mets{n});
    zoneIId.mets{metID} = strcat('BacteroidIId_', zoneIId.mets{metID});
end

% Add 'Nodule_' to the start of all plant nodule genes
for n = 1:length(medicagoModel.genes)
    geneID = findGeneIDs(zoneIId, medicagoModel.genes{n});
    zoneIId.genes{geneID} = strcat('NoduleIId_', zoneIId.genes{geneID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule genes
for n = 1:length(melilotiModel.genes)
    geneID = findGeneIDs(zoneIId, melilotiModel.genes{n});
    zoneIId.genes{geneID} = strcat('BacteroidIId_', zoneIId.genes{geneID});
end

% Update the grRules based on new gene names
for n = 1:length(noduleModel.grRules)
    for m = 1:length(noduleModel.genes)
        zoneIId.grRules{n} = strrep(zoneIId.grRules{n}, noduleModel.genes{m}, zoneIId.genes{m});
    end
end

%% Rename reactions, metabolites, and genes in the nodule zone IIp

% Add 'Nodule_' to the start of all plant nodule reactions
for n = 1:length(medicagoModel.rxns)
    rxnID = findRxnIDs(zoneIIp, medicagoModel.rxns{n});
    zoneIIp.rxns{rxnID} = strcat('NoduleIIp_', zoneIIp.rxns{rxnID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.rxns)
    rxnID = findRxnIDs(zoneIIp, melilotiModel.rxns{n});
    zoneIIp.rxns{rxnID} = strcat('BacteroidIIp_', zoneIIp.rxns{rxnID});
end

% Add 'Nodule_ to the start of the EXTC reactions
zoneIIp.rxns = strrep(zoneIIp.rxns, 'EXCT_', 'NoduleIIp_EXCT_');

% Add 'Nodule_' to the start of all plant nodule metabolites
for n = 1:length(medicagoModel.mets)
    metID = findMetIDs(zoneIIp, medicagoModel.mets{n});
    zoneIIp.mets{metID} = strcat('NoduleIIp_', zoneIIp.mets{metID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.mets)
    metID = findMetIDs(zoneIIp, melilotiModel.mets{n});
    zoneIIp.mets{metID} = strcat('BacteroidIIp_', zoneIIp.mets{metID});
end

% Add 'Nodule_' to the start of all plant nodule genes
for n = 1:length(medicagoModel.genes)
    geneID = findGeneIDs(zoneIIp, medicagoModel.genes{n});
    zoneIIp.genes{geneID} = strcat('NoduleIIp_', zoneIIp.genes{geneID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule genes
for n = 1:length(melilotiModel.genes)
    geneID = findGeneIDs(zoneIIp, melilotiModel.genes{n});
    zoneIIp.genes{geneID} = strcat('BacteroidIIp_', zoneIIp.genes{geneID});
end

% Update the grRules based on new gene names
for n = 1:length(noduleModel.grRules)
    for m = 1:length(noduleModel.genes)
        zoneIIp.grRules{n} = strrep(zoneIIp.grRules{n}, noduleModel.genes{m}, zoneIIp.genes{m});
    end
end

%% Rename reactions, metabolites, and genes in the nodule zone IZ

% Add 'Nodule_' to the start of all plant nodule reactions
for n = 1:length(medicagoModel.rxns)
    rxnID = findRxnIDs(zoneIZ, medicagoModel.rxns{n});
    zoneIZ.rxns{rxnID} = strcat('NoduleIZ_', zoneIZ.rxns{rxnID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.rxns)
    rxnID = findRxnIDs(zoneIZ, melilotiModel.rxns{n});
    zoneIZ.rxns{rxnID} = strcat('BacteroidIZ_', zoneIZ.rxns{rxnID});
end

% Add 'Nodule_ to the start of the EXTC reactions
zoneIZ.rxns = strrep(zoneIZ.rxns, 'EXCT_', 'NoduleIZ_EXCT_');

% Add 'Nodule_' to the start of all plant nodule metabolites
for n = 1:length(medicagoModel.mets)
    metID = findMetIDs(zoneIZ, medicagoModel.mets{n});
    zoneIZ.mets{metID} = strcat('NoduleIZ_', zoneIZ.mets{metID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.mets)
    metID = findMetIDs(zoneIZ, melilotiModel.mets{n});
    zoneIZ.mets{metID} = strcat('BacteroidIZ_', zoneIZ.mets{metID});
end

% Add 'Nodule_' to the start of all plant nodule genes
for n = 1:length(medicagoModel.genes)
    geneID = findGeneIDs(zoneIZ, medicagoModel.genes{n});
    zoneIZ.genes{geneID} = strcat('NoduleIZ_', zoneIZ.genes{geneID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule genes
for n = 1:length(melilotiModel.genes)
    geneID = findGeneIDs(zoneIZ, melilotiModel.genes{n});
    zoneIZ.genes{geneID} = strcat('BacteroidIZ_', zoneIZ.genes{geneID});
end

% Update the grRules based on new gene names
for n = 1:length(noduleModel.grRules)
    for m = 1:length(noduleModel.genes)
        zoneIZ.grRules{n} = strrep(zoneIZ.grRules{n}, noduleModel.genes{m}, zoneIZ.genes{m});
    end
end

%% Rename reactions, metabolites, and genes in the nodule zone III

% Add 'Nodule_' to the start of all plant nodule reactions
for n = 1:length(medicagoModel.rxns)
    rxnID = findRxnIDs(zoneIII, medicagoModel.rxns{n});
    zoneIII.rxns{rxnID} = strcat('NoduleIII_', zoneIII.rxns{rxnID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.rxns)
    rxnID = findRxnIDs(zoneIII, melilotiModel.rxns{n});
    zoneIII.rxns{rxnID} = strcat('BacteroidIII_', zoneIII.rxns{rxnID});
end

% Add 'Nodule_ to the start of the EXTC reactions
zoneIII.rxns = strrep(zoneIII.rxns, 'EXCT_', 'NoduleIII_EXCT_');

% Add 'Nodule_' to the start of all plant nodule metabolites
for n = 1:length(medicagoModel.mets)
    metID = findMetIDs(zoneIII, medicagoModel.mets{n});
    zoneIII.mets{metID} = strcat('NoduleIII_', zoneIII.mets{metID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule reactions
for n = 1:length(melilotiModel.mets)
    metID = findMetIDs(zoneIII, melilotiModel.mets{n});
    zoneIII.mets{metID} = strcat('BacteroidIII_', zoneIII.mets{metID});
end

% Add 'Nodule_' to the start of all plant nodule genes
for n = 1:length(medicagoModel.genes)
    geneID = findGeneIDs(zoneIII, medicagoModel.genes{n});
    zoneIII.genes{geneID} = strcat('NoduleIII_', zoneIII.genes{geneID});
end

% Add 'Bacteroid_' to the start of all bacterial nodule genes
for n = 1:length(melilotiModel.genes)
    geneID = findGeneIDs(zoneIII, melilotiModel.genes{n});
    zoneIII.genes{geneID} = strcat('BacteroidIII_', zoneIII.genes{geneID});
end

% Update the grRules based on new gene names
for n = 1:length(noduleModel.grRules)
    for m = 1:length(noduleModel.genes)
        zoneIII.grRules{n} = strrep(zoneIII.grRules{n}, noduleModel.genes{m}, zoneIII.genes{m});
    end
end
