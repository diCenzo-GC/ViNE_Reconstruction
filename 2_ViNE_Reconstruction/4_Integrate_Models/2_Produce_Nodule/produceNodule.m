%%

fprintf('\n\n 1. Loading and updating Truncatula model\n\n');
LoadTruncatula
medicagoModel = model;
save('medicagoModel.mat', 'medicagoModel');
clear medicagoModel;

%%

fprintf('\n\n 2. Loading Sinorhizobium model\n\n');
LoadSino
melilotiModel = sino;
save('melilotiModel.mat', 'melilotiModel');
clear melilotiModel;

%%

fprintf('\n\n 3. Executing bash pipeline\n\n');
BashPipeline

%%

fprintf('\n\n 4. Preparing Truncatula model for integration\n\n');
PrepareTruncatula

%%

fprintf('\n\n 5. Preparing Sinorhizobium model for integration\n\n');
PrepareSino

%%

fprintf('\n\n 6. Integrating Truncatula and Sinorhizobium models\n\n');
CombineModels

fprintf('\n\n 7. Set bounds and test growth\n\n');
CheckGrowth

%%

fprintf('\n\n 8. Return metabolites to original names\n\n');
renameMets

%%

fprintf('\n\n 9. Finalize the nodule model\n\n');
finalizeNodule

%%

fprintf('\n\n 10. Save and clear workspace\n\n');

save('allWorkspace.mat');
save('noduleModel.mat','noduleModel');
clear;
