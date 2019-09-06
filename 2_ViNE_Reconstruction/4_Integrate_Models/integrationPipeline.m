%% Initialize the system

% MATLAB R2016b
% libSBML version 5.13.0
% SBMLToolbox version 4.1.0
% COBRA Toolbox downloaded May 12, 2017
% Tn-Core Toolbox version 2.2
% TIGER Toolbox version 1.2-beta
% FASTCORE version 1.0
% iLOG CPLEX Studio version 12.7.1 

% Set up the system
clear all
addpath(genpath('../Software/cobratoolbox/'));
rmpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/SBMLToolbox-4.1.0/'));
addpath(genpath('../Software/Tn-Core-v2.2/'));
addpath(genpath('../Software/tiger/'));
addpath(genpath('../Software/Utilities'));
rmpath(genpath('/home/marcof/Data/Software/opencobra-cobratoolbox-bb71768/'));
initCobraToolbox;
addpath(genpath('../Software/FastCore/'));
rmpath(genpath('../Software/cobratoolbox/src/dataIntegration/transcriptomics/FASTCORE'));
rmpath(genpath('/home/marcof/Data/Software/opencobra-cobratoolbox-bb71768/'));
changeCobraSolver('ibm_cplex', 'all');

%% Produce the plant (root/shoot) model

cd 1_Produce_Plant
producePlant;
cd ..
copyfile('1_Produce_Plant/plantModel.mat', '3_Nodulate_Plant/plantModel.mat');

%% Produce the initial base nodule (medicago/sinorhizobium) model

cd 2_Produce_Nodule
produceNodule;
cd ..
copyfile('2_Produce_Nodule/noduleModel.mat', '3_Nodulate_Plant/noduleModel.mat');
copyfile('2_Produce_Nodule/medicagoModel.mat', '3_Nodulate_Plant/medicagoModel.mat');
copyfile('2_Produce_Nodule/melilotiModel.mat', '3_Nodulate_Plant/melilotiModel.mat');

%% Add the nodule to the plant

cd 3_Nodulate_Plant
nodulatePlant;
cd ..
copyfile('3_Nodulate_Plant/nodulatedPlant.mat', '4_Constrain_Nodule/nodulatedPlant.mat');

%% Constrain all nodule zones based on RNA-seq data of Roux et al (2014)

cd 4_Constrain_Nodule
constrainNodule;
cd ..
copyfile('4_Constrain_Nodule/constrainedNodule.mat', '5_Finalize_Nodule/constrainedNodule.mat');

%% Finalize the nodule format

cd 5_Finalize_Nodule
finalizeNodule;
cd ..
copyfile('5_Finalize_Nodule/finalNodulatedPlant.mat', 'finalNodulatedPlant.mat');

exit
