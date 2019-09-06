%% Prepare sino model 1

load('iGD1348_working.mat');

%% Add reaction to convert cpdFixed

iGD1348 = addReaction(iGD1348, {'rxnConvFixed', 'Convert fixed ammonium'}, ...
    {'cpdFixed[c0]', 'cpdSymCoF[c0]', 'cpd00013[c0]'}, [-1 -0.001 1], 0, 0, 1000, 0);

%% Change fluxes to umol / gDW / hr

% Modify the biomass reaction
[A B] = parseRxnFormula(cell2mat(printRxnFormula(iGD1348, 'rxnBIOMASS')));
C = B * 1000;
for n = 1:length(A)
    if strmatch('cpd11416', A{n})
        C(n) = 1;
    end
end
iGD1348 = addReaction(iGD1348, 'rxnBIOMASS', A, C, 0, 0, 1000, 1, '', iGD1348.grRules{findRxnIDs(iGD1348, 'rxnBIOMASS')});

% Multiply bounds by 1000
iGD1348.lb = iGD1348.lb * 1000;
iGD1348.ub = iGD1348.ub * 1000;

%% Reformat models so that reaction names are the same format

% Remove the _c0 and _e0 from the end of the reaction names
for n = 1:length(iGD1348.rxns);
    iGD1348.rxns{n} = strrep(iGD1348.rxns{n}, '_c0', '');
    iGD1348.rxns{n} = strrep(iGD1348.rxns{n}, '_e0', '');
end

% Add _e0 for the exchange reactions
for n = 1:length(iGD1348.rxns);
    if strfind(iGD1348.rxns{n}, 'EX_')
        iGD1348.rxns{n} = strcat(iGD1348.rxns{n}, '_e0');
    end
end

% Change model name
sino = iGD1348;

%% Modify input

sino.mets = strrep(sino.mets,'[c0]','_c0');
sino.mets = strrep(sino.mets,'[e0]','_e0');
sino.mets = strrep(sino.mets,'[b]','_b');


%get sino boundary metabolites (to be used to feed ShareMetsPipeline.sh)

boundary_mets = ~cellfun(@isempty, regexp(sino.mets,'_e0'));
sino_boundary_mets = sino.mets(boundary_mets);


fileID = fopen('Sino_boundary.mets','w');
fprintf(fileID,'%s\n' ,sino_boundary_mets{:});
fclose(fileID);

% modify protons
for n = 1:length(sino.mets)
    if strmatch('cpd00067_e0', sino.mets{n}, 'exact')
        sino.mets{n} = strrep(sino.mets{n}, '_e0', '_p0');
        sino.metNames{n} = strrep(sino.mets{n}, '_e0', '_p0');
    end
end
sino = addMetabolite(sino, 'cpd00067_e0', ...
    'H_e0', 'H', {[]}, {[]}, {[]}, {[]}, 1);
sino = changeRxnMets(sino, {'cpd00067_p0'}, ...
    {'cpd00067_e0'}, 'EX_cpd00067_e0');


% list of mets in sinorhizobium whose names are to be changed to MNX ones


%% Modify the unknown and spontaneous gene names

% Set initial variables
newUnknownNumber = 0;
newSpontaneousNumber = 0;
oldUnknown = findGeneIDs(sino, 'Unknown');
oldUnknownA = ['x(' mat2str(oldUnknown) '\)'];
oldUnknownB = ['x(' mat2str(oldUnknown) ')'];
oldSpontaneous= findGeneIDs(sino, 'Spontaneous');
oldSpontaneousA = ['x(' mat2str(oldSpontaneous) '\)'];
oldSpontaneousB = ['x(' mat2str(oldSpontaneous) ')'];

% Update the rules
for n = 1:length(sino.rxns)
    if strfind(sino.rules{n}, oldUnknownB)
        pos = strfind(sino.rules{n}, oldUnknownB);
        for m = 1:length(pos)
            newUnknownNumber = newUnknownNumber + 1;
            newUnknown = ['Unknown_' mat2str(newUnknownNumber)];
            newUnknownPos = length(sino.genes) + 1;
            newUnknownPos = ['x(' mat2str(newUnknownPos) '\)'];
            sino.genes{end+1} = newUnknown;
            sino.rules{n} = regexprep(sino.rules{n}, oldUnknownA, newUnknownPos, 'once');
        end
    end
    if strfind(sino.rules{n}, oldSpontaneousB)
        pos = strfind(sino.rules{n}, oldSpontaneousB);
        for m = 1:length(pos)
            newSpontaneousNumber = newSpontaneousNumber + 1;
            newSpontaneous = ['Spontaneous_' mat2str(newSpontaneousNumber)];
            newSpontaneousPos = length(sino.genes) + 1;
            newSpontaneousPos = ['x(' mat2str(newSpontaneousPos) '\)'];
            sino.genes{end+1} = newSpontaneous;
            sino.rules{n} = regexprep(sino.rules{n}, oldSpontaneousA, newSpontaneousPos, 'once');
        end
    end
end

% Update the grRules
sino.grRules = sino.rules;
for n = 1:length(sino.genes)
    pos = ['x(' mat2str(n) ')'];
    for m = 1:length(sino.grRules)
        sino.grRules{m} = strrep(sino.grRules{m}, pos, sino.genes{n});
    end
end
for n = 1:length(sino.rxns)
    sino.grRules{n} = strrep(sino.grRules{n}, '&', 'and');
    sino.grRules{n} = strrep(sino.grRules{n}, '|', 'or');
end

% Update the rxnGeneMat
sino.rxnGeneMat = cell(length(sino.rxns),length(sino.genes));
for n = 1:length(sino.rxns)
    if ~isempty(sino.rules{n})
        rules = strsplit(strrep(sino.rules{n}, '&', '|'), ' | ');
        for m = 1:length(rules)
            rules{m} = strrep(rules{m}, '(', '');
            rules{m} = strrep(rules{m}, ')', '');
            rules{m} = strrep(rules{m}, 'x', '');
            rules{m} = strrep(rules{m}, ' ', '');
        end
        for m = 1:length(rules)
            sino.rxnGeneMat{n,str2num(rules{m})} = 1;
        end
    end
end
isEmpty = cellfun('isempty',sino.rxnGeneMat);
sino.rxnGeneMat(isEmpty) = {0};
sino.rxnGeneMat = cell2mat(sino.rxnGeneMat);
sino.rxnGeneMat = sparse(double(sino.rxnGeneMat));

% Delete unused genes
unknownPos = findGeneIDs(sino, 'Unknown');
sino.genes(unknownPos) = [];
sino.rxnGeneMat(:,unknownPos) = [];
spontaneousPos = findGeneIDs(sino, 'Spontaneous');
sino.genes(spontaneousPos) = [];
sino.rxnGeneMat(:,spontaneousPos) = [];

% Update rules
sino.rules = sino.grRules;
for n = 1:length(sino.genes)
    pos = ['x(' mat2str(length(sino.genes) - n + 1) ')'];
    for m = 1:length(sino.rules)
        sino.rules{m} = strrep(sino.rules{m}, sino.genes{length(sino.genes) - n + 1}, pos);
    end
end
for n = 1:length(sino.rxns)
    sino.rules{n} = strrep(sino.rules{n}, 'and', '&');
    sino.rules{n} = strrep(sino.rules{n}, 'or', '|');
end

%%

IntegratedSino = sino;
