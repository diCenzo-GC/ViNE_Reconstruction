%% Load data

% Load gene ID conversion table
load('Data/finalConversions.mat');

% Load files to be updated
load('Data/GeneIDForPresenceCall.mat');
RootPresence = importdata(['Data' filesep 'RootPresenceCalls.csv']);
RootPresence = RootPresence.data;
StemPresence = importdata(['Data' filesep 'StemPresenceCalls.csv']);
StemPresence = StemPresence.data;
LeavePresence = importdata(['Data' filesep 'LeavePresenceCalls.csv']);
LeavePresence = LeavePresence.data;

%% Prepare updated files

% V5 gene list
genesNew = unique(finalConversions(:,3));

% Prepare variables
GeneIDnew = cell(length(genesNew),1);
LeavePresenceNew = cell2mat(cell(length(genesNew),1));
RootPresenceNew = cell2mat(cell(length(genesNew),1));
StemPresenceNew = cell2mat(cell(length(genesNew),1));

% Make a call table for v5 genes
for n = 1:length(genesNew)
    pos = strmatch(genesNew{n}, finalConversions(:,3), 'exact');
    oldGenes = finalConversions(pos, 1);
    for m = 1:length(oldGenes)
        pos2 = strmatch(oldGenes{m}, GeneID, 'exact');
        if ~isempty(pos2)
            GeneIDnew{n} = genesNew{n};
            LeaveCalls = LeavePresence(pos2);
            RootCalls = RootPresence(pos2);
            StemCalls = StemPresence(pos2);
            if max(LeaveCalls) == 1
                LeavePresenceNew(n) = 1;
            elseif min(LeaveCalls) == 1
                LeavePresenceNew(n) = -1;
            else
                LeavePresenceNew(n) = 0;
            end
            if max(RootCalls) == 1
                RootPresenceNew(n) = 1;
            elseif min(RootCalls) == 1
                RootPresenceNew(n) = -1;
            else
                RootPresenceNew(n) = 0;
            end
            if max(StemCalls) == 1
                StemPresenceNew(n) = 1;
            elseif min(StemCalls) == 1
                StemPresenceNew(n) = -1;
            else
                StemPresenceNew(n) = 0;
            end
        else
            GeneIDnew{n} = genesNew{n};
            LeavePresenceNew(n) = 0;
            RootPresenceNew(n) = 0;
            StemPresenceNew(n) = 0;
        end
    end
end
LeavePresenceNew = transpose(LeavePresenceNew);
RootPresenceNew = transpose(RootPresenceNew);
StemPresenceNew = transpose(StemPresenceNew);
LeavePresence = LeavePresenceNew;
RootPresence = RootPresenceNew;
StemPresence = StemPresenceNew;
GeneID = GeneIDnew;

%% Save the new files

save('Data/RootPresenceCallsNew.mat', 'RootPresence');
save('Data/LeavePresenceCallsNew.mat', 'LeavePresence');
save('Data/StemPresenceCallsNew.mat', 'StemPresence');
save('Data/GeneIDForPresenceCallNew.mat', 'GeneID');


