%% Prepare the Medicago model

% Taken from importMedicago script from Thomas Pfau
origstrings = {'__45__','__46__','__40__','__41__','__43__'};
targetstrings = {'-','.','(',')','+'};
for i=1:numel(targetstrings)
    medicagoModel.rxns = regexprep(medicagoModel.rxns,origstrings{i},targetstrings{i});
    medicagoModel.mets= regexprep(medicagoModel.mets,origstrings{i},targetstrings{i});
end
medicagoModel.rxns = regexprep(medicagoModel.rxns,'^_','');
medicagoModel.mets = regexprep(medicagoModel.mets,'^_','');
%Replace the name of the ATPase, which is originally RXN-11135_C
atpasepos = find(ismember(medicagoModel.rxns,'RXN-11135_C'));
medicagoModel.rxns{atpasepos} = 'ATPase';

%And change the BiomassWithOutStarch to BiomassShootWithOutStarch
biomasswostarchshootpos = find(ismember(medicagoModel.rxns,'BiomassWithOutStarch'));
medicagoModel.rxns{biomasswostarchshootpos} = 'BiomassShootWithOutStarch';

% Reduce asparagine in root biomass
RootBiomass = find(ismember(medicagoModel.rxns,'BiomassRoot'));
AsparagineRoot = find(ismember(medicagoModel.mets,{'ASN[C]'}));
ASNMolWeight = 132.12;
ASNRootAmounts = medicagoModel.S(AsparagineRoot,RootBiomass);
AsnRootWeights = abs(ASNMolWeight/1e6 * ASNRootAmounts);
AsnRootWeightChange = AsnRootWeights - AsnRootWeights * 0.05;
medicagoModel.S(AsparagineRoot,RootBiomass) = 0.05 * medicagoModel.S(AsparagineRoot,RootBiomass);
medicagoModel.S(:,RootBiomass) = medicagoModel.S(:,RootBiomass) / (1-AsnRootWeightChange);

%% Prepare the meliloti model

melilotiModel.mets = strrep(melilotiModel.mets, '_c0', '[c]');
melilotiModel.mets = strrep(melilotiModel.mets, '_e0', '[e]');
melilotiModel.mets = strrep(melilotiModel.mets, '_b', '[b]');
melilotiModel.mets = strrep(melilotiModel.mets, '_p0', '[p]');

%% Fix nodule model

pos = findMetIDs(noduleModel, 'cpd00067_p0[c]');
noduleModel.mets{pos} = strrep(noduleModel.mets{pos}, '_p0[c]', '[p]');

%% Remove unnecessary unknowns from the nodule model

% Find gene index number for Unknown
unknownID = findGeneIDs(noduleModel,'Unknown');

% Replace the unknowns that are not needed
A = [' | x(' num2str(unknownID) ')'];
B = ['x(' num2str(unknownID) ') | '];
noduleModel.grRules = strrep(noduleModel.grRules,' or Unknown','');
noduleModel.grRules = strrep(noduleModel.grRules,'Unknown or ','');
noduleModel.rules = strrep(noduleModel.rules,A,'');
noduleModel.rules = strrep(noduleModel.rules,B,'');

% Fix the rxnGeneMat
rxnGeneMatNewFull = cell(length(noduleModel.rxns),length(noduleModel.genes));
for n = 1:length(noduleModel.rxns)
    if ~isempty(noduleModel.rules{n})
        rulesTemp = strrep(noduleModel.rules{n}, '&', '|');
        rules = strsplit(rulesTemp, '|');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            rxnGeneMatNewFull{n,str2num(rules{m})} = 1;
        end
    end
end
isEmpty = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(isEmpty) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
noduleModel.rxnGeneMat = rxnGeneMatNew;

%% Remove unnecessary unknowns from the meliloti model

% Find gene index number for Unknown
unknownID = findGeneIDs(melilotiModel,'Unknown');

% Replace the unknowns that are not needed
A = [' | x(' num2str(unknownID) ')'];
B = ['x(' num2str(unknownID) ') | '];
melilotiModel.grRules = strrep(melilotiModel.grRules,' or Unknown','');
melilotiModel.grRules = strrep(melilotiModel.grRules,'Unknown or ','');
melilotiModel.rules = strrep(melilotiModel.rules,A,'');
melilotiModel.rules = strrep(melilotiModel.rules,B,'');

% Fix the rxnGeneMat
rxnGeneMatNewFull = cell(length(melilotiModel.rxns),length(melilotiModel.genes));
for n = 1:length(melilotiModel.rxns)
    if ~isempty(melilotiModel.rules{n})
        rulesTemp = strrep(melilotiModel.rules{n}, '&', '|');
        rules = strsplit(rulesTemp, '|');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            rxnGeneMatNewFull{n,str2num(rules{m})} = 1;
        end
    end
end
isEmpty = cellfun('isempty',rxnGeneMatNewFull);
rxnGeneMatNewFull(isEmpty) = {0};
rxnGeneMatNewFull = cell2mat(rxnGeneMatNewFull);
rxnGeneMatNew = sparse(double(rxnGeneMatNewFull));
melilotiModel.rxnGeneMat = rxnGeneMatNew;

%% Fix the nodule biomass

% Fix reaction in the medicago model
medicagoModel = tncore_remove_reactions(medicagoModel, 'BiomassRoot');
rootBiomass = printRxnFormula(plantModel, 'Root_BiomassRoot');
rootBiomass = strrep(rootBiomass, 'Root_', '');
medicagoModel = addReaction(medicagoModel, 'BiomassRoot', ...
    rootBiomass{1,1}, [], 0, 0, 1000, 0, [], [], [], [], [], 0);
medicagoModel = addReaction(medicagoModel, 'EX_BiomassRoot', ...
    {'BiomassRoot'}, [-1], 0, 0, 1000, 0, [], [], [], [], [], 0);
medicagoModel = changeObjective(medicagoModel, 'BiomassRoot');

% Fix the reaction in the nodule model
noduleModel = tncore_remove_reactions(noduleModel, 'BiomassRoot');
rootBiomass = printRxnFormula(plantModel, 'Root_BiomassRoot');
rootBiomass = strrep(rootBiomass, 'Root_', '');
noduleModel = addReaction(noduleModel, 'BiomassRoot', ...
    rootBiomass{1,1}, [], 0, 0, 1000, 0, [], [], [], [], [], 0);
noduleModel = addReaction(noduleModel, 'EX_BiomassRoot', ...
    {'BiomassRoot'}, [-1], 0, 0, 1000, 0, [], [], [], [], [], 0);
noduleModel = changeObjective(noduleModel, 'BiomassRoot');
