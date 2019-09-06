%% Prepare expression data and thresholds

% Remove Medicago genes with no expression
x = 0;
for n = 1:length(rna_medic)
    test = rna_medic{n,2} + rna_medic{n,3} + rna_medic{n,4} + rna_medic{n,5} + rna_medic{n,6};
    if test ~= 0
        x = x + 1;
        rna_medic_temp(x,:) = rna_medic(n,:);
    end
end
rna_medic = rna_medic_temp;

% Remove meliloti genes with no expression
x = 0;
for n = 1:length(rna_sino)
    test = rna_sino{n,2} + rna_sino{n,3} + rna_sino{n,4} + rna_sino{n,5} + rna_sino{n,6};
    if test ~= 0
        x = x + 1;
        rna_sino_temp(x,:) = rna_sino(n,:);
    end
end
rna_sino = rna_sino_temp;

% Set the expression threshold for Medicago
medic_thresh = mean(mean(cell2mat(rna_medic(:,2:end)))) * 1.1;

% Set the expression threshold for Sinorhizobium
sino_thresh = mean(mean(cell2mat(rna_sino(:,3:end)))) * 1.1;

%% Pre-process the RNA-seq data

% Pull out data for the model plant genes
plantGenes = nodulatedPlant.genes(strmatch('NoduleIId_', nodulatedPlant.genes));
for n= 1:length(plantGenes)
    plantGenes{n} = strrep(plantGenes{n}, 'NoduleIId_', '');
end
plant_rnaseq_model = {};
plant_rnaseq_model_avg = {};
for n = 1:length(plantGenes)
    pos = strmatch(plantGenes{n}, rna_medic_full(:,1), 'exact');
    if ~isempty(pos)
        if max(cell2mat(rna_medic_full(pos,2:end))) >= medic_thresh
            plant_rnaseq_model = vertcat(plant_rnaseq_model, rna_medic_full(pos,:));
            plant_rnaseq_model_avg = vertcat(plant_rnaseq_model_avg, rna_medic_original(pos,:));
        end
    end
end

% Pull out data for the model plant genes
bacterialGenes = nodulatedPlant.genes(strmatch('BacteroidIId_', nodulatedPlant.genes));
for n= 1:length(bacterialGenes)
    bacterialGenes{n} = strrep(bacterialGenes{n}, 'BacteroidIId_', '');
end
bacterium_rnaseq_model = {};
bacterium_rnaseq_model_avg = {};
for n = 1:length(bacterialGenes)
    pos = strmatch(bacterialGenes{n}, rna_sino_full(:,1), 'exact');
    if ~isempty(pos)
        if max(cell2mat(rna_sino_full(pos,2:end))) >= sino_thresh
            bacterium_rnaseq_model = vertcat(bacterium_rnaseq_model, rna_sino_full(pos,:));
            bacterium_rnaseq_model_avg = vertcat(bacterium_rnaseq_model_avg, rna_sino_original(pos,:));
        end
    end
end

% Write files to system
plant_rnaseq_model2 = cell2table(plant_rnaseq_model);
bacterium_rnaseq_model2 = cell2table(bacterium_rnaseq_model);
writetable(plant_rnaseq_model2, 'rnaSeqModel_plant.txt', 'Delimiter', '\t', ...
    'WriteVariableNames', false);
writetable(bacterium_rnaseq_model2, 'rnaSeqModel_bacterium.txt', 'Delimiter', '\t', ...
    'WriteVariableNames', false);

% Run the statistical analysis
system('./rnaSeqProcessing.r');

% Import the stats output
plant_stats = table2cell(readtable('rnaSeqModel_plant_stats.txt', 'ReadVariableNames', false));
bacterium_stats = table2cell(readtable('rnaSeqModel_bacterium_stats.txt', 'ReadVariableNames', false));

% Deal with the plant RNA-seq data
for n = 1:length(plant_rnaseq_model_avg)
    for o = 1:4
        for m = 2:6
            if plant_rnaseq_model_avg{n,m} >= medic_thresh
                groups = strsplit(plant_stats{n,m}, '');
                for i = 1:length(groups)
                    for j = 2:6
                        if strfind(plant_stats{n,j}, groups{i})
                            if plant_rnaseq_model_avg{n,j} < medic_thresh && plant_rnaseq_model_avg{n,j} >= 0.8 * medic_thresh
                                pos = strmatch(plant_rnaseq_model_avg{n,1}, rna_medic(:,1), 'exact');
                                rna_medic{pos,j} = plant_rnaseq_model_avg{n,m};
                            end
                        end
                    end
                end
            end
        end
    end
end

% Deal with the bacterial RNA-seq data
for n = 1:length(bacterium_rnaseq_model_avg)
    for o = 1:4
        for m = 2:6
            if bacterium_rnaseq_model_avg{n,m} >= sino_thresh
                groups = cellstr(bacterium_stats{n,m}');
                for i = 1:length(groups)
                    for j = 2:6
                        if strfind(bacterium_stats{n,j}, groups{i})
                            if bacterium_rnaseq_model_avg{n,j} < sino_thresh && bacterium_rnaseq_model_avg{n,j} >= 0.8 * sino_thresh
                                pos = strmatch(bacterium_rnaseq_model_avg{n,1}, rna_sino(:,1), 'exact');
                                rna_sino{pos,j} = bacterium_rnaseq_model_avg{n,m};
                                bacterium_rnaseq_model_avg{n,j} = bacterium_rnaseq_model_avg{n,m};
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Prepare the RNA-seq data

% Prepare variables for plant in each zone
medic_zI = horzcat(rna_medic(:,1), rna_medic(:,2));
medic_zIId = horzcat(rna_medic(:,1), rna_medic(:,3));
medic_zIIp = horzcat(rna_medic(:,1), rna_medic(:,4));
medic_zIZ = horzcat(rna_medic(:,1), rna_medic(:,5));
medic_zIII = horzcat(rna_medic(:,1), rna_medic(:,6));

% Prepare variables for plant in each zone
sino_zIId = horzcat(rna_sino(:,1), rna_sino(:,3));
sino_zIIp = horzcat(rna_sino(:,1), rna_sino(:,4));
sino_zIZ = horzcat(rna_sino(:,1), rna_sino(:,5));
sino_zIII = horzcat(rna_sino(:,1), rna_sino(:,6));

% Rename plant gene names in RNA seq variables
for n = 1:length(medic_zIId)
    medic_zI{n,1} = insertBefore(medic_zI{n,1}, 1, 'NoduleI_');
    medic_zIId{n,1} = insertBefore(medic_zIId{n,1}, 1, 'NoduleIId_');
    medic_zIIp{n,1} = insertBefore(medic_zIIp{n,1}, 1, 'NoduleIIp_');
    medic_zIZ{n,1} = insertBefore(medic_zIZ{n,1}, 1, 'NoduleIZ_');
    medic_zIII{n,1} = insertBefore(medic_zIII{n,1}, 1, 'NoduleIII_');
end

% Rename bacteria gene names in RNA seq variables
for n = 1:length(sino_zIId)
    sino_zIId{n,1} = insertBefore(sino_zIId{n,1}, 1, 'BacteroidIId_');
    sino_zIIp{n,1} = insertBefore(sino_zIIp{n,1}, 1, 'BacteroidIIp_');
    sino_zIZ{n,1} = insertBefore(sino_zIZ{n,1}, 1, 'BacteroidIZ_');
    sino_zIII{n,1} = insertBefore(sino_zIII{n,1}, 1, 'BacteroidIII_');
end

% Combine plant RNA-seq variables
medicago_rnaseq_data = vertcat(medic_zI, medic_zIId, medic_zIIp, medic_zIZ, medic_zIII);

% Combine bacterium RNA-seq variables
meliloti_rnaseq_data = vertcat(sino_zIId, sino_zIIp, sino_zIZ, sino_zIII);

%% Force the shoot and root genes to be on

% Identify the shoot and root genes
rootGenes = nodulatedPlant.genes(strmatch('Root_', nodulatedPlant.genes));
shootGenes = nodulatedPlant.genes(strmatch('Leave_', nodulatedPlant.genes));

% Set root and shoot expression above threshold
for n = 1:length(rootGenes)
    rootGenes{n,2} = 2 * medic_thresh;
end
for n = 1:length(shootGenes)
    shootGenes{n,2} = 2 * medic_thresh;
end

% Add root and shoot to the expression list
medicago_rnaseq_data = vertcat(medicago_rnaseq_data, rootGenes);
meliloti_rnaseq_data = vertcat(meliloti_rnaseq_data, shootGenes);
all_rnaseq_data = vertcat(medicago_rnaseq_data, meliloti_rnaseq_data);

%% Identify model RNA-seq data

% Pull out data for all model genes. If doesn't exist, set to 0
all_rnaseq_data_model = cell(length(nodulatedPlant.genes), 2);
for n = 1:length(nodulatedPlant.genes)
    pos = strmatch(nodulatedPlant.genes{n}, all_rnaseq_data(:,1), 'exact');
    if isempty(pos)
        all_rnaseq_data_model{n,1} = nodulatedPlant.genes{n,1};
        all_rnaseq_data_model{n,2} = 0;
    else
        all_rnaseq_data_model{n,1} = nodulatedPlant.genes{n,1};
        all_rnaseq_data_model{n,2} = all_rnaseq_data{pos,2};
    end
end

% Split into separate variables for the plant and bacteria
medicago_rnaseq_data_model = {};
meliloti_rnaseq_data_model = {};
for n = 1:length(all_rnaseq_data_model)
    if strmatch('Leave', all_rnaseq_data_model{n,1});
        medicago_rnaseq_data_model = vertcat(medicago_rnaseq_data_model, all_rnaseq_data_model(n,:));
    elseif strmatch('Root', all_rnaseq_data_model{n,1});
        medicago_rnaseq_data_model = vertcat(medicago_rnaseq_data_model, all_rnaseq_data_model(n,:));
    elseif strmatch('Nodule', all_rnaseq_data_model{n,1});
        medicago_rnaseq_data_model = vertcat(medicago_rnaseq_data_model, all_rnaseq_data_model(n,:));
    elseif strmatch('Bacteroid', all_rnaseq_data_model{n,1});
        meliloti_rnaseq_data_model = vertcat(meliloti_rnaseq_data_model, all_rnaseq_data_model(n,:));
    end
end
