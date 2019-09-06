%% Modify the model based on the importMedicago script of Thomas Pfau

% Fix naming conventions
origstrings = {'__45__','__46__','__40__','__41__','__43__'};
targetstrings = {'-','.','(',')','+'};
for i=1:numel(targetstrings)
    noduleModel.rxns = regexprep(noduleModel.rxns,origstrings{i},targetstrings{i});
    noduleModel.mets= regexprep(noduleModel.mets,origstrings{i},targetstrings{i});
end
noduleModel.rxns = regexprep(noduleModel.rxns,'^_','');
noduleModel.mets = regexprep(noduleModel.mets,'^_','');

% Replace the name of the ATPase, which is originally RXN-11135_C
atpasepos = find(ismember(noduleModel.rxns,'RXN-11135_C'));
noduleModel.rxns{atpasepos} = 'ATPase';

% Change the BiomassWithOutStarch to BiomassShootWithOutStarch
biomasswostarchshootpos = find(ismember(noduleModel.rxns,'BiomassWithOutStarch'));
noduleModel.rxns{biomasswostarchshootpos} = 'BiomassShootWithOutStarch';

%% Adjust the asparagine amount in biomass

ModelAdjustedForAsn = noduleModel;
RootBiomass = find(ismember(ModelAdjustedForAsn.rxns,'BiomassRoot'));
LeaveBiomass = find(ismember(ModelAdjustedForAsn.rxns,'BiomassShoot'));
LeaveBiomassWithoutStarch = find(ismember(ModelAdjustedForAsn.rxns,'BiomassShootWithOutStarch'));
AsparagineLeave = find(ismember(ModelAdjustedForAsn.mets,{'ASN[C]'}));
AsparagineRoot = find(ismember(ModelAdjustedForAsn.mets,{'ASN[C]'}));
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

noduleModel = ModelAdjustedForAsn;
