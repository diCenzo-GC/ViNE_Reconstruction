%% Constrain all nodule zones based on RNA-seq data of Roux et al (2014)

cd Constrain_Nodule
constrainNodule;
cd ..
copyfile('Constrain_Nodule/constrainedNodule_sucrose.mat', 'Finalize_Nodule/constrainedNodule_sucrose.mat');

%% Finalize the nodule format

cd Finalize_Nodule
finalizeNodule;
cd ..
copyfile('Finalize_Nodule/finalNodulatedPlant_sucrose.mat', 'Combine_Models/finalNodulatedPlant_sucrose.mat');

%% Combine models

cd Combine_Models
combineModels;
cd ..
copyfile('Combine_Models/combinedModel.mat', 'combinedModel.mat');
