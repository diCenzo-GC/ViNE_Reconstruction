%% Prepare sino model 2

sino = changeObjective(sino, 'rxn06874');
EXlist = sino.rxns(strmatch('EX_', sino.rxns));
sino = changeRxnBounds(sino,EXlist,0,'l');

sol = optimizeCbModel(sino);
PlantSinoEXlist = importdata('SinoNamesMapping_of_shared.ExListofShared');

sino_test = changeRxnBounds(sino, PlantSinoEXlist, -0.1, 'l');
sol = optimizeCbModel(sino_test)

%rename mets according to Metanetx codes (identified throught
%ShareMetsPipeline.sh)

    %import new and old names (files written by PrepareMetsSubstitution.pl)
    SinoNewNames = importdata('SinoNamesMapping_of_shared.NewMetNames');
    SinoOldNames = importdata('SinoNamesMapping_of_shared.OldMetNames');
    sharedMets = importdata('SinoNamesMapping_of_shared.ExListofShared_MNXcode');

    for i = 1:length(SinoOldNames)

        metindex=findMetIDs(IntegratedSino, SinoOldNames(i));
    
     % _b mets do not formally exist in the mat version of the model. the
     % following if statement avoids this issue
     
     if (metindex > 0)
          
                 IntegratedSino.mets(metindex) = SinoNewNames(i);
         
     end
     
    end
    
    %remove ex reactions (except demand reactions) and then add new EX
    %(EXPS) reactions for model cross-talk
    
    AllEx_IntegratedSino = IntegratedSino.rxns(findExcRxns(IntegratedSino));
    AllEx_IntegratedSino = AllEx_IntegratedSino(strmatch('EX', AllEx_IntegratedSino));
    sol = optimizeCbModel(IntegratedSino);
        

    %Close all EX reactions
    IntegratedSino_withEX = changeRxnBounds(IntegratedSino, AllEx_IntegratedSino, 0, 'l');
    IntegratedSino_Export = changeRxnBounds(IntegratedSino, AllEx_IntegratedSino, 0, 'l');
    %Find EX reactions involving a shared met
    
    AllReactionsWithMNXCompounds = findRxnsFromMets(IntegratedSino_withEX, sharedMets);
    AllEXWithMNXCompounds = AllReactionsWithMNXCompounds(~cellfun(@isempty, regexp(AllReactionsWithMNXCompounds,'EX_')));
    printRxnFormula(IntegratedSino_withEX, AllEXWithMNXCompounds, false);
    
    %Open them and test growth
    IntegratedSino_withEX = changeRxnBounds(IntegratedSino_withEX, AllEXWithMNXCompounds, -.1, 'l');
 
    IntegratedSino_withEX = changeObjective(IntegratedSino_withEX, 'rxn06874');
    sol = optimizeCbModel(IntegratedSino_withEX);

    IntegratedSino_Export = changeObjective(IntegratedSino_Export, 'rxn06874');
    sol = optimizeCbModel(IntegratedSino_Export);

    
    %change name to cross-talking EX reactions (from EX_ to EXCT_ and from cpd codes to MNX codes)
    noAtpRxns = {'EX_cpd00001_e0'; 'EX_cpd00011_e0'; 'EX_cpd00007_e0';...
        'EX_cpd00528_e0'; 'EX_cpd11640_e0'; 'EX_cpd00013_e0'};
    AllEXWithMNXCompoundsNewName = [];
    for j = 1:length(AllEXWithMNXCompounds)
        
       formula = printRxnFormula(IntegratedSino_withEX, AllEXWithMNXCompounds{j}, false);
       [mat,tok] = regexp(formula,'^(\w+)\_', 'match', 'tokens');
       mat{:};
       newname = strcat('EXCT_for_',mat{:},'e0');
       newnameB = strcat('EXCT_rev_',mat{:},'e0');
       met1 = strrep(mat{:}, '_', '_e0');
       met2 = strrep(mat{:}, '_', '[C]');
       met_atp = {'ATP[C]'};
       met_adp = {'ADP[C]'};
       met_pi = {'MNXM9[C]'};
       met_proton = {'MNXM01[C]'};
       rxnID = findRxnIDs( IntegratedSino_withEX,AllEXWithMNXCompounds{j});
       IntegratedSino_withEX.rxns{rxnID}  = newname{:};
       AllEXWithMNXCompoundsNewName{end+1} = cell2mat(newname);
       %This is the model that will be actually exported and integrated
       %into the truncatula one.
       if strmatch(AllEXWithMNXCompounds{j}, noAtpRxns, 'exact')
           if strmatch('EX_cpd00013_e0', AllEXWithMNXCompounds{j}, 'exact')
               met_proton2 = {'MNXM01_e0'};
               IntegratedSino_Export = addReaction(IntegratedSino_Export,newnameB{:},...
                   [met1, met_proton2, met2, met_proton],[-1 -1 1 1],false,0,...
                   [],[],[],[],[],[],[],0);
               IntegratedSino_Export = addReaction(IntegratedSino_Export,newname{:},...
                   [met2, met_atp, met1, met_adp, met_pi, met_proton],...
                   [-1 -0.25 1 0.25 0.25 0.25],false,0,[],[],[],[],[],[],[],0);
           else
               IntegratedSino_Export = addReaction(IntegratedSino_Export,newname{:},...
                   [met2, met1],[-1 1],false,0,[],[],[],[],[],[],[],0);
               IntegratedSino_Export = addReaction(IntegratedSino_Export,newnameB{:},...
                   [met1, met2],[-1 1],false,0,[],[],[],[],[],[],[],0);
           end
       else
           IntegratedSino_Export = addReaction(IntegratedSino_Export,newname{:},...
               [met2, met_atp, met1, met_adp, met_pi, met_proton],...
               [-1 -0.25 1 0.25 0.25 0.25],false,0,[],[],[],[],[],[],[],0);
           IntegratedSino_Export = addReaction(IntegratedSino_Export,newnameB{:},...
               [met1, met_atp, met2, met_adp, met_pi, met_proton],...
               [-1 -0.25 1 0.25 0.25 0.25],false,0,[],[],[],[],[],[],[],0);
       end

    end
    
    %% create array for addying exchange  reaction to Export model (for debugging purposes)
    AllNewMets = [];
    for j =  1:length(AllEXWithMNXCompoundsNewName);
        
     formula = printRxnFormula(IntegratedSino_withEX, AllEXWithMNXCompoundsNewName{j}, false);
       [mat,tok] = regexp(formula,'^(\w+)\_', 'match', 'tokens');

       met2 = strrep(mat{:}, '_', '[C]');
       AllNewMets{j} = met2{:};

end

    
    sol = optimizeCbModel(IntegratedSino_withEX);
    sol = optimizeCbModel(IntegratedSino_Export);


    CrossTalk_EX_reactions = IntegratedSino_Export.rxns(~cellfun(@isempty, regexp(IntegratedSino_Export.rxns,'EXCT_')));

    AllEX_sinoIntegrated =  IntegratedSino_withEX.rxns(findExcRxns(IntegratedSino_withEX));
    
    printRxnFormula(IntegratedSino_withEX, CrossTalk_EX_reactions, false);
    
    sol = optimizeCbModel(IntegratedSino_withEX);
    

    %% close all EX reactions except those that will be open during symbiosis and test growth
    
    TestSinoModel = IntegratedSino_withEX;
    AllEx_IntegratedSinoWithEX = TestSinoModel.rxns(findExcRxns(TestSinoModel));
    sol = optimizeCbModel(TestSinoModel);
    AllEX_reactions = AllEx_IntegratedSinoWithEX(~cellfun(@isempty, regexp(AllEx_IntegratedSinoWithEX,'EX')));

    TestSinoModel = changeRxnBounds(TestSinoModel, AllEX_reactions, 0, 'l');
    sol = optimizeCbModel(TestSinoModel);
    
    CrossTalk_in_CombinedModel = TestSinoModel.rxns(~cellfun(@isempty, regexp(TestSinoModel.rxns,'EXCT_')));
    TestSinoModel = changeRxnBounds(TestSinoModel , CrossTalk_in_CombinedModel, -.001, 'l');
    sol = optimizeCbModel(TestSinoModel);

    %% prepare model for models integration
    
    IntegratedSino_withEX.rxnsformula = printRxnFormula(IntegratedSino_withEX, TestSinoModel.rxns, false);
    IntegratedSino_Export.rxnsformula = printRxnFormula(IntegratedSino_Export, IntegratedSino_Export.rxns,false);
    sol = optimizeCbModel(IntegratedSino_Export);

    
    %% testing growth with all the nutrients that will be available inside the plant model
    AllLowers = repmat(-.001, 1, length(AllNewMets));
    AllUps = repmat(0, 1, length(AllNewMets));
    IntegratedSino_ExportDebug = addExchangeRxn(IntegratedSino_Export, AllNewMets, AllLowers, AllUps);
    sol = optimizeCbModel(IntegratedSino_ExportDebug);
