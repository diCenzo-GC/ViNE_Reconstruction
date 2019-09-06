%% Initialize the system

% MATLAB R2016b
% libSBML version 5.13.0
% SBMLToolbox version 4.1.0
% COBRA Toolbox downloaded May 12, 2017
% Tn-Core Toolbox version 2.2
% FASTCORE version 1.0
% TIGER Toolbox version 1.2.0-beta
% iLOG CPLEX Studio version 12.7.1

% Set up the system
clear all
addpath(genpath('../Software/cobratoolbox/'));
rmpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/cobratoolbox/external/libSBML-5.13.0-matlab_old'));
addpath(genpath('../Software/SBMLToolbox-4.1.0/'));
addpath(genpath('../Software/Tn-Core-v2.2/'));
addpath(genpath('../Software/tiger/'));
rmpath(genpath('/home/marcof/Data/Software/opencobra-cobratoolbox-bb71768/'));
initCobraToolbox;
addpath(genpath('../Software/FastCore/'));
rmpath(genpath('../Software/cobratoolbox/src/dataIntegration/transcriptomics/FASTCORE'));
rmpath(genpath('/home/marcof/Data/Software/opencobra-cobratoolbox-bb71768/'));
changeCobraSolver('ibm_cplex', 'all');

%% Load the models

% Load iGD1348
load('iGD1348_working.mat');

% Load iGD1575b
iGD1575 = readCbModel('iGD1575b.xml');

% Load iGD726
iGD726 = readCbModel('iGD726.xml');
iGD726 = deleteModelGenes(iGD726, 'NonEssential');
iGD726 = tncore_delete(iGD726);

%% Prepare iGD1348

% Get all exchange reactions
EX_list = iGD1348.rxns(strmatch('EX_', iGD1348.rxns));

% Medium exchange reactions
carbonFreeMinimal = {'EX_cpd00001_e0','EX_cpd00007_e0',...
    'EX_cpd00009_e0','EX_cpd00149_e0','EX_cpd00013_e0','EX_cpd00030_e0',...
    'EX_cpd00034_e0','EX_cpd00048_e0','EX_cpd00058_e0','EX_cpd00063_e0',...
    'EX_cpd00099_e0','EX_cpd00104_e0','EX_cpd00011_e0',...
    'EX_cpd00205_e0','EX_cpd00254_e0','EX_cpd00971_e0','EX_cpd10515_e0',...
    'EX_cpd10516_e0','EX_cpd00305_e0'};

% Set the medium composition without a carbon source
iGD1348 = changeRxnBounds(iGD1348, EX_list, 0, 'l');
iGD1348 = changeRxnBounds(iGD1348, carbonFreeMinimal, -10, 'l');

% Make sure the objective is set
iGD1348 = changeObjective(iGD1348, 'rxnBIOMASS');
iGD1348 = changeRxnBounds(iGD1348, 'rxnNGAM', 8.39, 'b');

%% Prepare iGD1575

% Get all exchange reactions
EX_list = iGD1575.rxns(strmatch('EX_', iGD1575.rxns));

% Medium exchange reactions
carbonFreeMinimal = {'EX_cpd00001_e0','EX_cpd00007_e0',...
    'EX_cpd00009_e0','EX_cpd00149_e0','EX_cpd00013_e0','EX_cpd00030_e0',...
    'EX_cpd00034_e0','EX_cpd00048_e0','EX_cpd00058_e0','EX_cpd00063_e0',...
    'EX_cpd00099_e0','EX_cpd00104_e0',...
    'EX_cpd00205_e0','EX_cpd00254_e0','EX_cpd00971_e0','EX_cpd10515_e0',...
    'EX_cpd10516_e0','EX_cpd00305_e0'};

% Set the medium composition without a carbon source
iGD1575 = changeRxnBounds(iGD1575, EX_list, 0, 'l');
iGD1575 = changeRxnBounds(iGD1575, carbonFreeMinimal, -10, 'l');

% Make sure the objective is set
iGD1575 = changeObjective(iGD1575, 'biomass_bulk_c0');

%% Prepare iGD726

% Get all exchange reactions
EX_list = iGD726.rxns(strmatch('EX_', iGD726.rxns));

% Medium exchange reactions
carbonFreeMinimal = {'EX_cpd00001_e0','EX_cpd00007_e0',...
    'EX_cpd00009_e0','EX_cpd00149_e0','EX_cpd00013_e0','EX_cpd00030_e0',...
    'EX_cpd00034_e0','EX_cpd00048_e0','EX_cpd00058_e0','EX_cpd00063_e0',...
    'EX_cpd00099_e0','EX_cpd00104_e0','EX_cpd00011_e0',...
    'EX_cpd00205_e0','EX_cpd00254_e0','EX_cpd00971_e0','EX_cpd10515_e0',...
    'EX_cpd10516_e0','EX_cpd00305_e0'};

% Set the medium composition without a carbon source
iGD726 = changeRxnBounds(iGD726, EX_list, 0, 'l');
iGD726 = changeRxnBounds(iGD726, carbonFreeMinimal, -10, 'l');

% Make sure the objective is set
iGD726 = changeObjective(iGD726, 'rxnBIOMASS');

%% Check effect of NGAM

% Storage variable
output_effectNGAM = cell(3,4);
output_effectNGAM{1,1} = 'Carbon_Source';
output_effectNGAM{1,2} = 'NGAM_0';
output_effectNGAM{1,3} = 'NGAM_8.39';
output_effectNGAM{1,4} = 'Ratio';
output_effectNGAM{2,1} = 'Glucose';
output_effectNGAM{3,1} = 'Succinate';

% Effect when grown with glucose
model2 = changeRxnBounds(iGD1348, 'EX_cpd00027_e0', -2.41, 'l');
model2 = changeRxnBounds(model2, 'rxnNGAM', 0, 'b');
sol = optimizeCbModel(model2, 'max');
output_effectNGAM{2,2} = sol.f;
model2 = changeRxnBounds(model2, 'rxnNGAM', 8.39, 'b');
sol = optimizeCbModel(model2, 'max');
output_effectNGAM{2,3} = sol.f;

% Effect when grown with succinate
model2 = changeRxnBounds(iGD1348, 'EX_cpd00036_e0', -3.921, 'l');
model2 = changeRxnBounds(model2, 'rxnNGAM', 0, 'b');
sol = optimizeCbModel(model2, 'max');
output_effectNGAM{3,2} = sol.f;
model2 = changeRxnBounds(model2, 'rxnNGAM', 8.39, 'b');
sol = optimizeCbModel(model2, 'max');
output_effectNGAM{3,3} = sol.f;

% Compare growth rates
output_effectNGAM{2,4} = output_effectNGAM{2,3} / output_effectNGAM{2,2};
output_effectNGAM{3,4} = output_effectNGAM{3,3} / output_effectNGAM{3,2};

% Clear extra variables
clearvars -except output_* iGD*;

%% Run the biolog analysis

% Define list of the models
testModels = {iGD1348; iGD1575; iGD726};

% Model names
models = {'iGD1348'; 'iGD1575'; 'iGD726'};

% Define list of the biolog test compounds
biolog = {'EX_cpd00020_e0','EX_cpd00023_e0','EX_cpd00027_e0',...
    'EX_cpd00029_e0','EX_cpd00033_e0','EX_cpd00035_e0','EX_cpd00036_e0',...
    'EX_cpd00041_e0','EX_cpd00047_e0','EX_cpd00051_e0','EX_cpd00053_e0',...
    'EX_cpd00054_e0','EX_cpd00064_e0','EX_cpd00069_e0','EX_cpd00076_e0',...
    'EX_cpd00082_e0','EX_cpd00100_e0','EX_cpd00105_e0','EX_cpd00106_e0',...
    'EX_cpd00107_e0','EX_cpd00108_e0','EX_cpd00119_e0','EX_cpd00121_e0',...
    'EX_cpd00129_e0','EX_cpd00130_e0','EX_cpd00138_e0','EX_cpd00142_e0',...
    'EX_cpd00154_e0','EX_cpd00155_e0','EX_cpd00156_e0','EX_cpd00158_e0',...
    'EX_cpd00159_e0','EX_cpd00161_e0','EX_cpd00179_e0','EX_cpd00182_e0',...
    'EX_cpd00208_e0','EX_cpd00224_e0','EX_cpd00246_e0','EX_cpd00281_e0',...
    'EX_cpd00308_e0','EX_cpd00314_e0','EX_cpd00322_e0','EX_cpd00366_e0',...
    'EX_cpd00382_e0','EX_cpd00386_e0','EX_cpd00392_e0','EX_cpd00396_e0',...
    'EX_cpd00417_e0','EX_cpd00567_e0','EX_cpd00588_e0','EX_cpd00794_e0',...
    'EX_cpd00797_e0','EX_cpd00851_e0','EX_cpd01133_e0','EX_cpd01307_e0',...
    'EX_cpd02175_e0','EX_cpd03198_e0','EX_cpd05158_e0','EX_cpd11585_e0',...
    'EX_cpd11588_e0','EX_cpd11589_e0','EX_cpd11592_e0','EX_cpd00024_e0',...
    'EX_cpd00040_e0','EX_cpd00060_e0','EX_cpd00066_e0','EX_cpd00072_e0',...
    'EX_cpd00079_e0','EX_cpd00089_e0','EX_cpd00094_e0','EX_cpd00136_e0',...
    'EX_cpd00137_e0','EX_cpd00139_e0','EX_cpd00141_e0','EX_cpd00157_e0',...
    'EX_cpd00164_e0','EX_cpd00180_e0','EX_cpd00211_e0','EX_cpd00212_e0',...
    'EX_cpd00222_e0','EX_cpd00232_e0','EX_cpd00248_e0','EX_cpd00266_e0',...
    'EX_cpd00280_e0','EX_cpd00306_e0','EX_cpd00320_e0','EX_cpd00339_e0',...
    'EX_cpd00361_e0','EX_cpd00374_e0','EX_cpd00380_e0','EX_cpd00432_e0',...
    'EX_cpd00438_e0','EX_cpd00453_e0','EX_cpd00477_e0','EX_cpd00489_e0',...
    'EX_cpd00550_e0','EX_cpd00599_e0','EX_cpd00607_e0','EX_cpd00609_e0',...
    'EX_cpd00611_e0','EX_cpd00666_e0','EX_cpd00728_e0','EX_cpd00750_e0',...
    'EX_cpd01055_e0','EX_cpd01067_e0','EX_cpd01107_e0','EX_cpd01113_e0',...
    'EX_cpd01246_e0','EX_cpd01363_e0','EX_cpd01502_e0','EX_cpd01799_e0',...
    'EX_cpd02143_e0','EX_cpd02351_e0','EX_cpd02599_e0','EX_cpd03161_e0',...
    'EX_cpd03561_e0','EX_cpd03696_e0','EX_cpd03734_e0','EX_cpd03737_e0',...
    'EX_cpd05161_e0','EX_cpd05192_e0','EX_cpd05240_e0','EX_cpd09457_e0',...
    'EX_cpd10719_e0','EX_cpd11602_e0','EX_cpd11685_e0','EX_cpd11717_e0',...
    'EX_cpd11879_e0','EX_cpd13391_e0','EX_cpd13392_e0','EX_cpd16821_e0',...
    'EX_cpd00039_e0','EX_cpd00122_e0','EX_cpd00132_e0','EX_cpd00162_e0',...
    'EX_cpd00185_e0','EX_cpd00249_e0','EX_cpd00276_e0','EX_cpd00492_e0',...
    'EX_cpd00589_e0','EX_cpd00737_e0','EX_cpd00751_e0','EX_cpd00832_e0',...
    'EX_cpd01030_e0','EX_cpd01171_e0','EX_cpd01200_e0','EX_cpd01262_e0',...
    'EX_cpd02274_e0','EX_cpd03884_e0','EX_cpd04349_e0','EX_cpd11594_e0',...
    'EX_cpd11601_e0','EX_cpd11603_e0','EX_cpd11748_e0','EX_cpd00080_e0',...
    'EX_cpd00117_e0','EX_cpd00118_e0','EX_cpd00184_e0','EX_cpd00227_e0',...
    'EX_cpd01242_e0','EX_cpd01293_e0','EX_cpd00197_e0','EX_cpd01524_e0'};

% Prepare output matrix
results = cell(length(biolog)+1,length(testModels)+3);

% Add headers to output
results{1,1} = 'exchangeRxn';
results{1,2} = 'cpdNumber';
results{1,3} = 'metaboliteName';

% Get metabolite names
for n = 1:length(biolog)
    metabolite = strrep(biolog{n}, 'EX_', '');
    metabolite = strrep(metabolite, '_e0', '[e0]');
    pos = findMetIDs(iGD1575, metabolite);
    metNames{n,1} = iGD1575.metNames{pos};
    metNames{n,2} = iGD1575.mets{pos};
    metNames{n,2} = strrep(metNames{n,2}, '[e0]', '');
end

% Do phenotype microarray and save output
for n = 1:length(testModels)
    results{1,n+3} = models{n};
    for m = 1:length(biolog)
        tmpTest = changeRxnBounds(testModels{n}, biolog{m}, -5, 'l');
        growth = optimizeCbModel(tmpTest, 'max');
        results{m+1,n+3} = round(growth.f, 4);
        results{m+1,2} = metNames{m,2};
        results{m+1,3} = metNames{m,1};
        results{m+1,1} = biolog{m};
    end
end

% Store the results
output_biolog = results;

% Clear extra variables
clearvars -except output_* iGD*;

%% Perform the Tn-Core analysis

% Import data
tnseq = table2cell(readtable('TnSeqData.txt'));

% Set carbon source
iGD1348 = changeRxnBounds(iGD1348, 'EX_cpd00076_e0', -1.2, 'l');
iGD1575 = changeRxnBounds(iGD1575, 'EX_cpd00076_e0', -1.2, 'l');
iGD726 = changeRxnBounds(iGD726, 'EX_cpd00076_e0', -1.2, 'l');

% Run Tn-Core without RNA-seq
[output_coreModel_iGD1348] = tncore_core(iGD1348, tnseq, [], ...
    [], [], [], [], [], 1);
[output_coreModel_iGD1575] = tncore_core(iGD1575, tnseq, [], ...
    [], [], [], [], [], 1);

% Clear extra variables
clearvars -except output_* iGD*;

%% Analyze the Tn-Core output

% Set output variable
output_TnCore = cell(1,1);
output_TnCore{1,1} = 'Model';
output_TnCore{2,1} = 'iGD1348';
output_TnCore{3,1} = 'iGD1575';
output_TnCore{1,2} = 'Genes_in_iGD726';
output_TnCore{1,3} = 'Percent_core_in_iG726';
output_TnCore{1,4} = 'Ess_genes_in_iGD726';
output_TnCore{1,5} = 'Ess_genes_in_ess_iGD726';
output_TnCore{1,6} = 'Percent_ess_core_in_iG726';
output_TnCore{1,7} = 'Percent_ess_overlap';

% Compare gene sets
output_TnCore{2,2} = length(intersect(output_coreModel_iGD1348.genes, iGD726.genes));
output_TnCore{3,2} = length(intersect(output_coreModel_iGD1575.genes, iGD726.genes));
output_TnCore{2,3} = output_TnCore{2,2} / length(output_coreModel_iGD1348.genes);
output_TnCore{3,3} = output_TnCore{3,2} / length(output_coreModel_iGD1575.genes);

% Get essential iGD1348_core gene lists
grRatio = singleGeneDeletion(output_coreModel_iGD1348);
test = isnan(grRatio);
grRatio(test) = 0;
x = 0;
for n = 1:length(grRatio)
    if round(grRatio(n),2) >= 0.01
        grRatio(n) = 1;
    end
end
iGD1348_core_ess = output_coreModel_iGD1348.genes(~logical(grRatio));

% Get essential iGD1575_core gene lists
grRatio = singleGeneDeletion(output_coreModel_iGD1575);
test = isnan(grRatio);
grRatio(test) = 0;
x = 0;
for n = 1:length(grRatio)
    if round(grRatio(n),2) >= 0.01
        grRatio(n) = 1;
    end
end
iGD1575_core_ess = output_coreModel_iGD1575.genes(~logical(grRatio));

% Get essential iGD726 gene lists
grRatio = singleGeneDeletion(iGD726);
test = isnan(grRatio);
grRatio(test) = 0;
x = 0;
for n = 1:length(grRatio)
    if round(grRatio(n),2) >= 0.01
        grRatio(n) = 1;
    end
end
iGD726_core_ess = iGD726.genes(~logical(grRatio));

% Comapre more lists
output_TnCore{2,4} = length(intersect(iGD1348_core_ess, iGD726.genes));
output_TnCore{3,4} = length(intersect(iGD1575_core_ess, iGD726.genes));
output_TnCore{2,5} = length(intersect(iGD1348_core_ess, iGD726_core_ess));
output_TnCore{3,5} = length(intersect(iGD1575_core_ess, iGD726_core_ess));
output_TnCore{2,6} = output_TnCore{2,5} / length(output_coreModel_iGD1348.genes);
output_TnCore{3,6} = output_TnCore{3,5} / length(output_coreModel_iGD1575.genes);
output_TnCore{2,7} = output_TnCore{2,5} / length(iGD726_core_ess);
output_TnCore{3,7} = output_TnCore{3,5} / length(iGD726_core_ess);

% Tncore overlap variable
overlap_TnCore = cell(3,8);
overlap_TnCore{1,1} = 'Overlap_type';
overlap_TnCore{2,1} = 'Full';
overlap_TnCore{3,1} = 'Essential';
overlap_TnCore{1,2} = 'iGD726_unique';
overlap_TnCore{1,3} = 'iGD1575_unique';
overlap_TnCore{1,4} = 'iGD1348_unique';
overlap_TnCore{1,5} = 'iGD726_iGD1575';
overlap_TnCore{1,6} = 'iGD726_iGD1348';
overlap_TnCore{1,7} = 'iGD13575_iGD1348';
overlap_TnCore{1,8} = 'iGD726_iGD13575_iGD1348';

% Get overlap data for all genes
A = intersect(iGD726.genes, output_coreModel_iGD1575.genes);
overlap_TnCore{2,8} = length(intersect(A, output_coreModel_iGD1348.genes)) - 2;
overlap_TnCore{2,7} = length(intersect(output_coreModel_iGD1575.genes, output_coreModel_iGD1348.genes)) - overlap_TnCore{2,8} - 2;
overlap_TnCore{2,6} = length(intersect(iGD726.genes, output_coreModel_iGD1348.genes)) - overlap_TnCore{2,8} - 2;
overlap_TnCore{2,5} = length(intersect(output_coreModel_iGD1575.genes, iGD726.genes)) - overlap_TnCore{2,8} - 2;
overlap_TnCore{2,4} = length(output_coreModel_iGD1348.genes) - overlap_TnCore{2,8} - overlap_TnCore{2,7} - overlap_TnCore{2,6} - 2;
overlap_TnCore{2,3} = length(output_coreModel_iGD1575.genes) - overlap_TnCore{2,8} - overlap_TnCore{2,7} - overlap_TnCore{2,5} - 2;
overlap_TnCore{2,2} = length(iGD726.genes) - overlap_TnCore{2,8} - overlap_TnCore{2,6} - overlap_TnCore{2,5} - 2;

% Get overlap data for all genes
A = intersect(iGD726_core_ess, iGD1575_core_ess);
overlap_TnCore{3,8} = length(intersect(A, iGD1348_core_ess)) - 2;
overlap_TnCore{3,7} = length(intersect(iGD1575_core_ess, iGD1348_core_ess)) - overlap_TnCore{3,8} - 2;
overlap_TnCore{3,6} = length(intersect(iGD726_core_ess, iGD1348_core_ess)) - overlap_TnCore{3,8} - 2;
overlap_TnCore{3,5} = length(intersect(iGD1575_core_ess, iGD726_core_ess)) - overlap_TnCore{3,8} - 2;
overlap_TnCore{3,4} = length(iGD1348_core_ess) - overlap_TnCore{3,8} - overlap_TnCore{3,7} - overlap_TnCore{3,6} - 2;
overlap_TnCore{3,3} = length(iGD1575_core_ess) - overlap_TnCore{3,8} - overlap_TnCore{3,7} - overlap_TnCore{3,5} - 2;
overlap_TnCore{3,2} = length(iGD726_core_ess) - overlap_TnCore{3,8} - overlap_TnCore{3,6} - overlap_TnCore{3,5} - 2;

% Clear extra variables
clearvars -except output_* overlap_* iGD*;
clearvars *core_ess

% Save
save('modelComparisonOutput.mat');
clear
quit
