%% Zone-specific single gene deletion analysis in the nodule

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
zoneSpecific_singleGeneDeletion;

% Clean the workspace
clearvars -except output*

%% Zone-specific single reaction deletion analysis in the nodule

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
zoneSpecific_singleRxnDeletion;

% Clean the workspace
clearvars -except output*

%% Tissue-specific plant gene deletion in the shoot and root

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
tissueSpecific_singleGeneDeletion;

% Clean the workspace
clearvars -except output_*

%% Tissue-specific plant reaction deletion in the shoot and root

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
tissueSpecific_singleRxnDeletion;

% Clean the workspace
clearvars -except output_*

%% Tissue-independent gene deletion throughout the entire model

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
tissueIndependent_singleGeneDeletion;

% Clean the workspace
clearvars -except output_*

%% Tissue-independent reaction deletion throughout the entire model

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
tissueIndependent_singleRxnDeletion;

% Clean the workspace
clearvars -except output_*

%% Zone-specific double gene deletion in the nodule

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
zoneSpecific_doubleGeneDeletion;

% Clean the workspace
clearvars -except output_*

%% Tissue-independent double gene deletion throughout the entire model

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
tissueIndependent_doubleGeneDeletion;

% Clean the workspace
clearvars -except output_*

%% Flux distribution patterns throughout the entire model

% Load the model
load('../finalNodulatedPlant.mat');
model = finalNodulatedPlant;

% Run the analysis
wholePlant_fluxDistribution;

% Clean the workspace
clearvars -except output_*

%% Save the output

save('essentialMetabolismOutput.mat');
clear
