%% Initialize the system

% MATLAB R2016b
% Gurobi version 7.0.1
% libSBML version 5.13.0
% SBMLToolbox version 4.1.0
% COBRA Toolbox downloaded May 12, 2017
% Tn-Core Toolbox version 2.2

% Set up the system
clear all
addpath(genpath('../../Software/cobratoolbox/'));
rmpath(genpath('../../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../Software/SBMLToolbox-4.1.0/'));
addpath(genpath('../../Software/Tn-Core-v2.2/'));
initCobraToolbox;
changeCobraSolver('gurobi', 'all');

%% Expand the model

% Load the models
model = readCbModel('manually_expanded_model.xml');
iGD1575 = readCbModel('iGD1575b.xml');

% Fix format of iGD1575 metabolites
for n = 1:length(iGD1575.mets)
    iGD1575.mets{n} = strrep(iGD1575.mets{n}, '_c0[c0]', '[c0]');
    iGD1575.mets{n} = strrep(iGD1575.mets{n}, '_e0[e0]', '[e0]');
end

% Fix format of new model metabolites (I think not necessary, but I have
% just in case
for n = 1:length(model.mets)
    model.mets{n} = strrep(model.mets{n}, '_c0[c0]', '[c0]');
    model.mets{n} = strrep(model.mets{n}, '_e0[e0]', '[e0]');
end

% Run tncore_expand
model_expanded = tncore_expand(model, iGD1575, 0);

% Remove the biomass reactions transferred from iGD1575
model_expanded = removeRxns(model_expanded, 'biomass_bulk_c0');
model_expanded = removeRxns(model_expanded, 'biomass_rhizo_c0');

% Remove new demand reactions

% Remove dead-end producing reactions
model_expanded = tncore_deadends(model_expanded);

%% Save and export the model

save('allWorkspace_expansion.mat');
save('model_expanded_noBalancing.mat', 'model_expanded');
tncore_export(model_expanded, model_expanded);

%% Manual steps

% Make only one kind of unknown, and remove non-essential GPR
% Add the new reactions to the existing Excel file.
% Add the chemical formulas and charges
% Then proceed to the massCharge balancing step
