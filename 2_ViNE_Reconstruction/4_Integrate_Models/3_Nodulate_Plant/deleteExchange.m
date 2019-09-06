%% Delete all exchange reactions from the nodule zone I

% Get exchange reactions
exchangeRxns = medicagoModel.rxns(strmatch('TEC_', medicagoModel.rxns));
exchangeRxns = vertcat(exchangeRxns, medicagoModel.rxns(strmatch('TCE_', medicagoModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, medicagoModel.rxns(strmatch('TGE_', medicagoModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, medicagoModel.rxns(strmatch('TEH_', medicagoModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, medicagoModel.rxns(strmatch('THE_', medicagoModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, medicagoModel.rxns(strmatch('EX_', medicagoModel.rxns)));
exchangeRxns = setdiff(exchangeRxns, 'EX_BiomassRoot');
exchangeRxns = setdiff(exchangeRxns, 'EX_cpd11416_e0');

% Add prefixes
for n = 1:length(exchangeRxns)
    if strmatch('EX_', exchangeRxns{n})
        exchangeRxns{n} = strcat('BacteroidI_', exchangeRxns{n});
    else
        exchangeRxns{n} = strcat('NoduleI_', exchangeRxns{n});
    end
end

% Delete the exchange reactions
zoneI_noEx = tncore_remove_reactions(zoneI, exchangeRxns);

%% Delete all exchange reactions from the nodule zone IId

% Get exchange reactions
exchangeRxns = noduleModel.rxns(strmatch('TEC_', noduleModel.rxns));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TCE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TGE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TEH_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('THE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('EX_', noduleModel.rxns)));
exchangeRxns = setdiff(exchangeRxns, 'EX_BiomassRoot');
exchangeRxns = setdiff(exchangeRxns, 'EX_cpd11416_e0');

% Add prefixes
for n = 1:length(exchangeRxns)
    if strmatch('EX_', exchangeRxns{n})
        exchangeRxns{n} = strcat('BacteroidIId_', exchangeRxns{n});
    else
        exchangeRxns{n} = strcat('NoduleIId_', exchangeRxns{n});
    end
end

% Delete the exchange reactions
zoneIId_noEx = tncore_remove_reactions(zoneIId, exchangeRxns);

%% Delete all exchange reactions from the nodule zone IIp

% Get exchange reactions
exchangeRxns = noduleModel.rxns(strmatch('TEC_', noduleModel.rxns));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TCE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TGE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TEH_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('THE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('EX_', noduleModel.rxns)));
exchangeRxns = setdiff(exchangeRxns, 'EX_BiomassRoot');
exchangeRxns = setdiff(exchangeRxns, 'EX_cpd11416_e0');

% Add prefixes
for n = 1:length(exchangeRxns)
    if strmatch('EX_', exchangeRxns{n})
        exchangeRxns{n} = strcat('BacteroidIIp_', exchangeRxns{n});
    else
        exchangeRxns{n} = strcat('NoduleIIp_', exchangeRxns{n});
    end
end

% Delete the exchange reactions
zoneIIp_noEx = tncore_remove_reactions(zoneIIp, exchangeRxns);

%% Delete all exchange reactions from the nodule zone IZ

% Get exchange reactions
exchangeRxns = noduleModel.rxns(strmatch('TEC_', noduleModel.rxns));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TCE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TGE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TEH_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('THE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('EX_', noduleModel.rxns)));
exchangeRxns = setdiff(exchangeRxns, 'EX_BiomassRoot');
exchangeRxns = setdiff(exchangeRxns, 'EX_cpd11416_e0');

% Add prefixes
for n = 1:length(exchangeRxns)
    if strmatch('EX_', exchangeRxns{n})
        exchangeRxns{n} = strcat('BacteroidIZ_', exchangeRxns{n});
    else
        exchangeRxns{n} = strcat('NoduleIZ_', exchangeRxns{n});
    end
end

% Delete the exchange reactions
zoneIZ_noEx = tncore_remove_reactions(zoneIZ, exchangeRxns);

%% Delete all exchange reactions from the nodule zone III

% Get exchange reactions
exchangeRxns = noduleModel.rxns(strmatch('TEC_', noduleModel.rxns));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TCE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TGE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('TEH_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('THE_', noduleModel.rxns)));
exchangeRxns = vertcat(exchangeRxns, noduleModel.rxns(strmatch('EX_', noduleModel.rxns)));
exchangeRxns = setdiff(exchangeRxns, 'EX_BiomassRoot');
exchangeRxns = setdiff(exchangeRxns, 'EX_cpd11416_e0');
exchangeRxns = setdiff(exchangeRxns, 'TCE_HYDROGEN_GAS');
exchangeRxns = setdiff(exchangeRxns, 'TEC_NITROGEN_GAS');

% Add prefixes
for n = 1:length(exchangeRxns)
    if strmatch('EX_', exchangeRxns{n})
        exchangeRxns{n} = strcat('BacteroidIII_', exchangeRxns{n});
    else
        exchangeRxns{n} = strcat('NoduleIII_', exchangeRxns{n});
    end
end

% Delete the exchange reactions
zoneIII_noEx = tncore_remove_reactions(zoneIII, exchangeRxns);