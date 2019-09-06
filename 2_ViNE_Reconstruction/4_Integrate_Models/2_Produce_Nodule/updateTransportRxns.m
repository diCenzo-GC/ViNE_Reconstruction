%% Find transport reactions

% Save original model
modelOrig = model;

% Identify transport reactions
reactions = model.rxns(strmatch('TEC_', model.rxns));
reactions = vertcat(reactions, model.rxns(strmatch('TEH_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCE_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TGE_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('THE_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCG_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCH_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCR_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCV_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCX_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TGC_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('THC_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TMC_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TCR_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TVC_', model.rxns)));
reactions = vertcat(reactions, model.rxns(strmatch('TXC_', model.rxns)));
reactions = unique(reactions);

% Exclude gases and some others
reactionsToExclude = {'TCE_CARBON-DIOXIDE'; 'TCE_OXYGEN-MOLECULE'; 'TCE_WATER'; ...
    'TCR_CO2'; 'TCR_OXYGEN'; 'TCR_WATER'; 'TEC_CARBON-DIOXIDE'; 'TEC_OXYGEN-MOLECULE'; ...
    'TEC_WATER'; 'TEH_Light'; 'TGC_OXYGEN-MOLECULE'; 'TGC_WATER'; 'THC_CO2'; ...
    'THC_O2'; 'THC_WATER'; 'TMC_CO2'; 'TMC_O2'; 'TMC_WATER'; 'TXC_CO2'; 'TXC_O2';...
    'TEC_NITROGEN_GAS'; 'TCE_HYDROGEN_GAS'};
reactions = setdiff(reactions, reactionsToExclude);

%% Break reactions into lists

% Identify reactions with one compound
model.S = full(model.S);
reactionsTemp = {};
for n = 1:length(reactions)
    pos = findRxnIDs(model, reactions{n});
    if sum(abs(model.S(:,pos))) < 3
        reactionsTemp = vertcat(reactionsTemp, reactions{n});
    end
end
reactions = reactionsTemp;

% Separate reactions with a proton
protonRxns = findRxnsFromMets(model, 'PROTON[C]');
protonRxns = intersect(protonRxns, reactions);
reactions = setdiff(reactions, protonRxns);

% Separate reactions with a phosphate molecule
phosphateRxns = findRxnsFromMets(model, 'Pi[C]');
phosphateRxns = intersect(phosphateRxns, reactions);
reactions = setdiff(reactions, phosphateRxns);

% Separate reactions with a proton
atpRxns = findRxnsFromMets(model, 'ATP[C]');
atpRxns = intersect(atpRxns, reactions);
reactions = setdiff(reactions, atpRxns);

% Separate reactions with a proton
adpRxns = findRxnsFromMets(model, 'ADP[C]');
adpRxns = intersect(adpRxns, reactions);
reactions = setdiff(reactions, adpRxns);

%% Add ATP to reactions with one compound exchanged

ID_atp = findMetIDs(model, 'ATP[C]');
ID_adp = findMetIDs(model, 'ADP[C]');
ID_proton = findMetIDs(model, 'PROTON[C]');
ID_phosphate = findMetIDs(model, 'Pi[C]');
for n = 1:length(reactions)
    formula = printRxnFormula(model, reactions{n}, false);
    if contains(formula, '<=>')
        pos = findRxnIDs(model, reactions{n});
        mets = model.mets(logical(abs(model.S(:,pos))));
        comp1 = reactions{n}(2);
        comp2 = reactions{n}(3);
        metName = strsplit(reactions{n}, '_');
        rxn1 = strcat('T', comp1, comp2, model.rxns{pos}(4:end));
        rxn2 = strcat('T', comp2, comp1, model.rxns{pos}(4:end));
        rxnName1 = strcat(model.rxnNames{pos}, '_1');
        rxnName2 = strcat(model.rxnNames{pos}, '_2');
        model = tncore_remove_reactions(model, reactions{n});
        model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        pos = findRxnIDs(model, rxn1);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 0.25;
        pos = findRxnIDs(model, rxn2);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 0.25;
    else
        pos = findRxnIDs(model, reactions{n});
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 0.25;
    end
end

%% Deal with proton containing reactions

for n = 1:length(protonRxns)
    if strmatch('TC', protonRxns{n})
        pos = findRxnIDs(model, protonRxns{n});
        mets = model.mets(logical(abs(model.S(:,pos))));
        comp1 = protonRxns{n}(2);
        comp2 = protonRxns{n}(3);
        metName = strsplit(protonRxns{n}, '_');
        rxn1 = strcat('T', comp1, comp2, model.rxns{pos}(4:end));
        rxn2 = strcat('T', comp2, comp1, model.rxns{pos}(4:end));
        rxnName1 = strcat(model.rxnNames{pos}, '_1');
        rxnName2 = strcat(model.rxnNames{pos}, '_2');
        if length(mets) == 1
            model = tncore_remove_reactions(model, protonRxns{n});
            model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}}, [-1], ...
                false, 0, [], [], [], [], [], [], false, 0);
            model = addReaction(model, {rxn2, rxnName2}, {mets{1,1}}, [1], ...
                false, 0, [], [], [], [], [], [], false, 0);
        else
            model = tncore_remove_reactions(model, protonRxns{n});
            model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
                false, 0, [], [], [], [], [], [], false, 0);
            model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
                false, 0, [], [], [], [], [], [], false, 0);
        end
        pos = findRxnIDs(model, rxn1);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = -0.75;
        pos = findRxnIDs(model, rxn2);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 1.25;
    elseif strmatch('TE', protonRxns{n})
        pos = findRxnIDs(model, protonRxns{n});
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 0.25;
        model.S(ID_proton,pos) = 1.25;
    end
end

%% Deal with phosphate containing reactions

for n = 1:length(phosphateRxns)
    if strmatch('TEC', phosphateRxns)
        pos = findRxnIDs(model, phosphateRxns{n});
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 1.25;
        model.S(ID_proton,pos) = 0.25;
    else
        pos = findRxnIDs(model, phosphateRxns{n});
        mets = model.mets(logical(abs(model.S(:,pos))));
        comp1 = phosphateRxns{n}(2);
        comp2 = phosphateRxns{n}(3);
        metName = strsplit(phosphateRxns{n}, '_');
        rxn1 = strcat('T', comp1, comp2, model.rxns{pos}(4:end));
        rxn2 = strcat('T', comp2, comp1, model.rxns{pos}(4:end));
        rxnName1 = strcat(model.rxnNames{pos}, '_1');
        rxnName2 = strcat(model.rxnNames{pos}, '_2');
        model = tncore_remove_reactions(model, reactions{n});
        model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
            false, 0, [], [], [], [], [], [], false, 0);
        pos = findRxnIDs(model, rxn1);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = -0.75;
        model.S(ID_proton,pos) = 0.25;
        pos = findRxnIDs(model, rxn2);
        model.S(ID_atp,pos) = -0.25;
        model.S(ID_adp,pos) = 0.25;
        model.S(ID_phosphate,pos) = 1.25;
        model.S(ID_proton,pos) = 0.25;
    end
end

%% Deal with ATP containing reactions

for n = 1:length(atpRxns)
    pos = findRxnIDs(model, atpRxns{n});
    mets = model.mets(logical(abs(model.S(:,pos))));
    comp1 = atpRxns{n}(2);
    comp2 = atpRxns{n}(3);
    metName = strsplit(atpRxns{n}, '_');
        rxn1 = strcat('T', comp1, comp2, model.rxns{pos}(4:end));
        rxn2 = strcat('T', comp2, comp1, model.rxns{pos}(4:end));
    rxnName1 = strcat(model.rxnNames{pos}, '_1');
    rxnName2 = strcat(model.rxnNames{pos}, '_2');
    model = tncore_remove_reactions(model, atpRxns{n});
    model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    pos = findRxnIDs(model, rxn1);
    model.S(ID_atp,pos) = -1.25;
    model.S(ID_adp,pos) = 0.25;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
    pos = findRxnIDs(model, rxn2);
    model.S(ID_atp,pos) = 0.75;
    model.S(ID_adp,pos) = 0.25;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
end

%% Deal with ADP containing reactions

for n = 1:length(adpRxns)
    pos = findRxnIDs(model, adpRxns{n});
    mets = model.mets(logical(abs(model.S(:,pos))));
    comp1 = adpRxns{n}(2);
    comp2 = adpRxns{n}(3);
    metName = strsplit(adpRxns{n}, '_');
        rxn1 = strcat('T', comp1, comp2, model.rxns{pos}(4:end));
        rxn2 = strcat('T', comp2, comp1, model.rxns{pos}(4:end));
    rxnName1 = strcat(model.rxnNames{pos}, '_1');
    rxnName2 = strcat(model.rxnNames{pos}, '_2');
    model = tncore_remove_reactions(model, adpRxns{n});
    model = addReaction(model, {rxn1, rxnName1}, {mets{1,1}, mets{2,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    model = addReaction(model, {rxn2, rxnName2}, {mets{2,1}, mets{1,1}}, [-1 1], ...
        false, 0, [], [], [], [], [], [], false, 0);
    pos = findRxnIDs(model, rxn1);
    model.S(ID_atp,pos) = -0.25;
    model.S(ID_adp,pos) = -0.75;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
    pos = findRxnIDs(model, rxn2);
    model.S(ID_atp,pos) = -0.25;
    model.S(ID_adp,pos) = 1.25;
    model.S(ID_phosphate,pos) = 0.25;
    model.S(ID_proton,pos) = 0.25;
end
model.S = sparse(model.S);