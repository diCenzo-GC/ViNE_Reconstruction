%% Prepare truncatula model

% load truncatula model
model = readCbModel('MedicagoTruncatula.xml');
model = changeObjective(model, 'BiomassShoot');

% Set all bounds to 1000000 and -1000000
for n = 1:length(model.lb)
    if model.lb(n) <= -1000;
        model.lb(n) = -1000000;       
    end
    if model.ub(n) >= 1000;
        model.ub(n) = 1000000;        
    end
end

% Add a N2 gas exchange
model = addReaction(model, 'TEC_NITROGEN_GAS', ...
    {'NITROGEN__45__MOLECULE[C]'}, [1], false, 0, 1000);

% Add a H2 gas exchange
model = addReaction(model, 'TCE_HYDROGEN_GAS', ...
    {'HYDROGEN__45__MOLECULE[C]'}, [-1], false, 0, 1000);

% Add missing transport reactions
model = addReaction(model, 'TEC_MN+2', {'MN+2[C]'}, [1], 0, 0, 1000, 0);
model = addReaction(model, 'TEC_ZN+2', {'ZN+2[C]'}, [1], 0, 0, 1000, 0);
model = addReaction(model, 'TEC_CA+2', {'CA+2[C]'}, [1], 0, 0, 1000, 0);
model = addReaction(model, 'TEC_K+', {'K+[C]'}, [1], 0, 0, 1000, 0);
model = addReaction(model, 'TEC_NA+', {'NA+[C]'}, [1], 0, 0, 1000, 0);

% Add biotin demand reaction
model = addReaction(model, 'Demand_BIOTIN', {'BIOTIN[C]'}, [1], 0, 0, 1000, 0);

% Update transport reactions to have ATP usage
updateTransportRxns;

%% Update the truncatula model

updateModel;
model = updatedModel;

% Add the homocitrate reaction
model = addReaction(model, 'HOMOCITRATE_SYNTHASE_RXN_C', ...
    {'WATER[C]', '2__45__KETOGLUTARATE[C]', 'ACETYL__45__COA[C]', ...
    'PROTON[C]', 'HOMO__45__CIT[C]', 'CO__45__A[C]'}, ...
    [-1 -1 -1 1 1 1], true, -1000, 1000, 0, [], 'MtrunA17Chr1g0213481');

%%
exchangers = ~cellfun(@isempty, regexp(model.rxns,'^T([A-Z]E|E[A-Z])|Cons'));
exchangersRed = ~cellfun(@isempty, regexp(model.rxns,'^TE[A-Z]|Cons'));

% model = changeObjective(model, model.rxns(findRxnIDs(model,'BiomassWithStarch')))
% sol = optimizeCbModel(model)



modelwithClosedBounds = changeRxnBounds(model,model.rxns(exchangers),0,'b');


%reset any existing objective and set the BiomassWithStarch as objective
modelwithClosedBounds.c(:) = 0;
modelwithClosedBounds.c(end-1) = 1;
sol = optimizeCbModel(modelwithClosedBounds);


% Get truncatula cellular metabolites (only those in the cytosol can serve
% as input for the sino model, right?)
trunca_boundary_mets = ~cellfun(@isempty, regexp(model.mets,'\[C\]'));
truncatula_cytosol_mets = model.mets(trunca_boundary_mets);

truncatula_cytosol_mets_ok = [];

for i=1:length(truncatula_cytosol_mets)
    
        
      
      NewMetNameTrunca = strrep(truncatula_cytosol_mets(i), '[', '_');
      NewMetNameTrunca = strrep(NewMetNameTrunca, ']', '');

     
      truncatula_cytosol_mets_ok{end+1} = cell2mat(NewMetNameTrunca);
end
    
fileID = fopen('Truncatula_Cytosol.mets','w');
fprintf(fileID,'%s\n' ,truncatula_cytosol_mets_ok{:});
fclose(fileID);

IntegratedTrunca = model;
