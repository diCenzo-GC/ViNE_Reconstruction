%% Find the transfer reactions

reactions = medicagoModel.rxns(strmatch('TRS_', medicagoModel.rxns));
reactions = vertcat(reactions, medicagoModel.rxns(strmatch('TSR_', medicagoModel.rxns)));
rxnsToExclude = {'TRS_WATER'; 'TRS_PROTON'; 'TSR_WATER'; 'TSR_PROTON'};
reactions = setdiff(reactions, rxnsToExclude);

%% Add the ATP requirement

% Get IDs of metabolites to add
ID_atp_root = findMetIDs(medicagoModel, 'Root_ATP[C]');
ID_adp_root = findMetIDs(medicagoModel, 'Root_ADP[C]');
ID_proton_root = findMetIDs(medicagoModel, 'Root_PROTON[C]');
ID_phosphate_root = findMetIDs(medicagoModel, 'Root_Pi[C]');
ID_atp_shoot = findMetIDs(medicagoModel, 'Leave_ATP[C]');
ID_adp_shoot = findMetIDs(medicagoModel, 'Leave_ADP[C]');
ID_proton_shoot = findMetIDs(medicagoModel, 'Leave_PROTON[C]');
ID_phosphate_shoot = findMetIDs(medicagoModel, 'Leave_Pi[C]');

% Add the metabolties
medicagoModel.S = full(medicagoModel.S);
for n = 1:length(reactions)
    if strmatch('TRS_AMMONIUM', reactions{n}, 'exact')
        pos = findRxnIDs(medicagoModel, reactions{n});
        medicagoModel.S(ID_atp_root,pos) = -0.25;
        medicagoModel.S(ID_adp_root,pos) = 0.25;
        medicagoModel.S(ID_phosphate_root,pos) = 0.25;
        medicagoModel.S(ID_proton_root,pos) = 0.25;
        medicagoModel.S(ID_atp_shoot,pos) = -0.25;
        medicagoModel.S(ID_adp_shoot,pos) = 0.25;
        medicagoModel.S(ID_phosphate_shoot,pos) = 0.25;
        medicagoModel.S(ID_proton_shoot,pos) = 0.25;
    elseif strmatch('TRS_Pi', reactions{n}, 'exact')
        pos = findRxnIDs(medicagoModel, reactions{n});
        medicagoModel.S(ID_atp_root,pos) = -0.25;
        medicagoModel.S(ID_adp_root,pos) = 0.25;
        medicagoModel.S(ID_phosphate_root,pos) = -0.75;
        medicagoModel.S(ID_proton_root,pos) = 0.25;
        medicagoModel.S(ID_atp_shoot,pos) = -0.25;
        medicagoModel.S(ID_adp_shoot,pos) = 0.25;
        medicagoModel.S(ID_phosphate_shoot,pos) = 1.25;
        medicagoModel.S(ID_proton_shoot,pos) = 0.25;
    else
        pos = findRxnIDs(medicagoModel, reactions{n});
        medicagoModel.S(ID_atp_root,pos) = -0.25;
        medicagoModel.S(ID_adp_root,pos) = 0.25;
        medicagoModel.S(ID_phosphate_root,pos) = 0.25;
        medicagoModel.S(ID_proton_root,pos) = 0.25;
        medicagoModel.S(ID_atp_shoot,pos) = -0.25;
        medicagoModel.S(ID_adp_shoot,pos) = 0.25;
        medicagoModel.S(ID_phosphate_shoot,pos) = 0.25;
        medicagoModel.S(ID_proton_shoot,pos) = 0.25;
    end
end
medicagoModel.S = sparse(medicagoModel.S);
