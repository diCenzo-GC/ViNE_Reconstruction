%% Initialize the system

% MATLAB R2016b
% libSBML version 5.13.0
% SBMLToolbox version 4.1.0
% COBRA Toolbox downloaded May 12, 2017
% Tn-Core Toolbox version 2.1b
% TIGER Toolbox version 1.2-beta
% FASTCORE version 1.0
% iLOG CPLEX Studio version 12.7.1 

% Set up the system
clear all
addpath(genpath('../Software/cobratoolbox/'));
rmpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/SBMLToolbox-4.1.0/'));
addpath(genpath('../Software/Tn-Core-v2.1b/'));
addpath(genpath('../Software/tiger/'));
addpath(genpath('../Software/Utilities'));
addpath(genpath('../Software/FastSL'));
rmpath(genpath('/home/marcof/Data/Software/opencobra-cobratoolbox-bb71768/'));
initCobraToolbox;
addpath(genpath('../Software/FastCore/'));
rmpath(genpath('../Software/cobratoolbox/src/dataIntegration/transcriptomics/FASTCORE'));
rmpath(genpath('/home/marcof/Data/Software/opencobra-cobratoolbox-bb71768/'));
changeCobraSolver('ibm_cplex', 'all');
addpath(genpath('../Software/DFBAlab/'));

%% Perform some preliminary model validation

cd 1_Model_Validation
modelValidationPipeline;
cd ..

%% Analyze trade-offs between N2-fixation and plant growth

cd 2_N2-fix_Tradeoffs
fixationTradeoffsPipeline;
cd ..

%% Evaluate gene and reaction essentiality

cd 3_Essential_Metabolism
essentialMetabolismPipeline;
cd ..

%% Examine bacteroid carbon source effect

cd 4_Carbon_Source
carbonSourcePipeline;
cd ..

%%

exit

