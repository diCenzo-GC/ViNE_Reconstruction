%% Import Medicago model

% Import model with the script of Thomas Pfau
model = importMedicago();

% Set all bounds to 1000000 and -1000000
for n = 1:length(model.lb)
    if model.lb(n) <= -1000;
        model.lb(n) = -1000000;       
    end
    if model.ub(n) >= 1000;
        model.ub(n) = 1000000;        
    end
end

%% Update the model to the new genome version

% Add ATP usage to all transport reactions
updateTransportRxns;

% Update the model to the new genome version
updateModel;
updateGeneIDs;

%% Produce Tissues

% Build the tissue specific model
medicagoModel = BuildTissueModel(updatedModel);

% Add energy requirement for inter-tissue transfer
updateTransferRxns;

% Missing transport reactions required for nitrogen fixation
medicagoModel = addReaction(medicagoModel, 'Root_TEC_CO+2', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_CO+2[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);
medicagoModel = addReaction(medicagoModel, 'Root_TEC_CPD-3', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_CPD-3[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);

% Missing transport reactions required for bacterial growth
medicagoModel = addReaction(medicagoModel, 'Root_TEC_MN+2', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_MN+2[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);
medicagoModel = addReaction(medicagoModel, 'Root_TEC_ZN+2', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_ZN+2[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);
medicagoModel = addReaction(medicagoModel, 'Root_TEC_CA+2', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_CA+2[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);
medicagoModel = addReaction(medicagoModel, 'Root_TEC_K+', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_K+[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);
medicagoModel = addReaction(medicagoModel, 'Root_TEC_NA+', ...
    {'Root_ATP[C]', 'Root_ADP[C]', 'Root_Pi[C]', 'Root_PROTON[C]', 'Root_NA+[C]'}, ...
    [-0.25 0.25 0.25 0.25 1], 0, 0, 1000, 0);

%% Adjust asparagine in the biomass

ModelAdjustedForAsn = medicagoModel;
RootBiomass = find(ismember(ModelAdjustedForAsn.rxns,'Root_BiomassRoot'));
LeaveBiomass = find(ismember(ModelAdjustedForAsn.rxns,'Leave_BiomassShoot'));
LeaveBiomassWithoutStarch = find(ismember(ModelAdjustedForAsn.rxns,'Leave_BiomassShootWithOutStarch'));
AsparagineLeave = find(ismember(ModelAdjustedForAsn.mets,{'Leave_ASN[C]'}));
AsparagineRoot = find(ismember(ModelAdjustedForAsn.mets,{'Root_ASN[C]'}));
ASNMolWeight = 132.12;

%First do it for the Asparagine without starch
ASNLeaveWOStarchAmounts = ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
AsnLeaveWOStarchWeights = abs(ASNMolWeight/1e6 * ASNLeaveWOStarchAmounts);
AsnLeaveWOStarchWeightChange = AsnLeaveWOStarchWeights - AsnLeaveWOStarchWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch) = 0.05 * ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
ModelAdjustedForAsn.S(:,LeaveBiomassWithoutStarch) = ModelAdjustedForAsn.S(:,LeaveBiomassWithoutStarch) / (1-AsnLeaveWOStarchWeightChange);

% And repeat it for that with starch
ASNLeaveAmounts = ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomassWithoutStarch);
AsnLeaveWeights = abs(ASNMolWeight/1e6 * ASNLeaveAmounts);
AsnLeaveWeightChange = AsnLeaveWeights - AsnLeaveWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomass) = 0.05 * ModelAdjustedForAsn.S(AsparagineLeave,LeaveBiomass);
ModelAdjustedForAsn.S(:,LeaveBiomass) = ModelAdjustedForAsn.S(:,LeaveBiomass) / (1-AsnLeaveWeightChange);

%And repeat the process for the root
ASNRootAmounts = ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass);
AsnRootWeights = abs(ASNMolWeight/1e6 * ASNRootAmounts);
AsnRootWeightChange = AsnRootWeights - AsnRootWeights * 0.05;
ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass) = 0.05 * ModelAdjustedForAsn.S(AsparagineRoot,RootBiomass);
ModelAdjustedForAsn.S(:,RootBiomass) = ModelAdjustedForAsn.S(:,RootBiomass) / (1-AsnRootWeightChange);

plantModel = ModelAdjustedForAsn;

%% Make tissue specific gene names

% Make models for each tissue
rootModel = plantModel;
shootModel = plantModel;

% Rename genes in each model
for n = 1:length(plantModel.genes)
    geneID = findGeneIDs(rootModel, plantModel.genes{n});
    rootModel.genes{geneID} = strcat('Root_', rootModel.genes{geneID});
    shootModel.genes{geneID} = strcat('Leave_', shootModel.genes{geneID});
end
for n = 1:length(plantModel.grRules)
    for m = 1:length(plantModel.genes)
        rootModel.grRules{n} = strrep(rootModel.grRules{n}, plantModel.genes{m}, rootModel.genes{m});
        shootModel.grRules{n} = strrep(shootModel.grRules{n}, plantModel.genes{m}, shootModel.genes{m});
    end
end

% Get list or root and shoot reactions
rootRxns = plantModel.rxns(strmatch('Root_', plantModel.rxns));
shootRxns = plantModel.rxns(strmatch('Leave_', plantModel.rxns));

% Get necessary info to add root to the shoot
rootRxnNames = plantModel.rxnNames(strmatch('Root_', plantModel.rxns));
plantFormulas = printRxnFormula(plantModel);
rootFormulas = plantFormulas(strmatch('Root_', plantModel.rxns));
rootRev = plantModel.rev(strmatch('Root_', plantModel.rxns));
rootLb = plantModel.lb(strmatch('Root_', plantModel.rxns));
rootUb = plantModel.ub(strmatch('Root_', plantModel.rxns));
rootC = plantModel.c(strmatch('Root_', plantModel.rxns));
rootSubsystem = plantModel.subSystems(strmatch('Root_', plantModel.rxns));
rootGrrules = rootModel.grRules(strmatch('Root_', plantModel.rxns));

% Delete root reactions from the shoot models
shootModel = tncore_remove_reactions(shootModel, rootRxns);

% Add root reactions back to the shoot model
rebuiltModel = shootModel;
for n = 1:length(rootRxns)
    rebuiltModel = addReaction(rebuiltModel, {rootRxns{n}, rootRxnNames{n}}, rootFormulas{n}, ...
        [], rootRev(n), rootLb(n), rootUb(n), rootC(n), rootSubsystem{n}, ...
        rootGrrules{n}, [], [], false);
end
rebuiltModel.metChEBIID = cell(1712,1);

% Remove any unused genes (if they exist)
plantModel = removeUnusedGenes(rebuiltModel);

% Set objective to biomass
plantModel = changeObjective(plantModel, 'Biomass');

% Fix bounds
for n = 1:length(plantModel.lb)
    if plantModel.lb(n) <= -1000;
        plantModel.lb(n) = -1000;       
    end
    if plantModel.ub(n) >= 1000;
        plantModel.ub(n) = 1000;        
    end
end

% Test growth
final_sol = optimizeCbModel(plantModel)

%% Save and clean workspace

save('allWorkspace.mat');
save('plantModel.mat', 'plantModel');
clear;

