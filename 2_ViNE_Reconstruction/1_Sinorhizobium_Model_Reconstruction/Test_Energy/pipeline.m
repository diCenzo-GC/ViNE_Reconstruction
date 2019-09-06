% Set up the system
clear all
addpath(genpath('../../Software/cobratoolbox/'));
rmpath(genpath('../../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../../Software/SBMLToolbox-4.1.0/'));
addpath(genpath('../../Software/Tn-Core-v2.2/'));
initCobraToolbox;
changeCobraSolver('gurobi', 'all');

% Load model
model = readCbModel('testModel.xml');

% Get exchange reactions
exRxns = model.rxns(strmatch('EX_', model.rxns));

% Minimal medium composition
carbonFreeMinimal = { 'EX_cpd00001_e0', 'EX_cpd00007_e0', 'EX_cpd00009_e0', ...
'EX_cpd00149_e0', 'EX_cpd00013_e0', 'EX_cpd00030_e0', 'EX_cpd00034_e0', ...
'EX_cpd00048_e0', 'EX_cpd00058_e0', 'EX_cpd00063_e0', 'EX_cpd00099_e0', ...
'EX_cpd00104_e0', 'EX_cpd00149_e0', 'EX_cpd00205_e0', 'EX_cpd00254_e0', ...
'EX_cpd00971_e0', 'EX_cpd10515_e0', 'EX_cpd10516_e0', 'EX_cpd00305_e0' };

% Set growth medium
model = changeRxnBounds(model, exRxns, 0, 'l')'
model = changeRxnBounds(model, carbonFreeMinimal, -100, 'l');
model = changeRxnBounds(model, 'EX_cpd00027_e0', -2.41, 'l');

% Set obj fun
model = changeObjective(model, 'rxnBIOMASS');

% Change NGAM
model1 = changeRxnBounds(model, 'rxnNGAM', 0, 'b');
model2 = changeRxnBounds(model, 'rxnNGAM', 8.39, 'b');
model3 = changeRxnBounds(model, 'rxnNGAM', 839, 'b');

% Test growth and record fluxes
sol1 = optimizeCbModel(model1);
sol2 = optimizeCbModel(model2);
%sol3 = optimizeCbModel(model3);
%fluxes = table(model1.rxns, sol1.x, sol2.x, sol3.x);
sol1.f
sol2.f
