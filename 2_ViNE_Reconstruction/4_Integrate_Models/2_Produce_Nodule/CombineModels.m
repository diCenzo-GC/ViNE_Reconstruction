%  createModel Create a COBRA model from inputs or an empty model
%  structure if no inputs are provided.
%  
%   model = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList,...
%      lowerBoundList,upperBoundList,subSystemList,grRuleList,geneNameList,...
%      systNameList)
%  
%  INPUTS
%   rxnAbrList            List of names of the new reactions
%   rxnNameList           List of names of the new reactions
%   rxnList               List of reactions: format: {'A -> B + 2 C'}
%                         If the compartment of a metabolite is not
%                         specified, it is assumed to be cytoplasmic, i.e. [c]
%                         
%                         
%                         OPTIONAL INPUTS
%   revFlagList           List of reversibility flag (opt, default = 1)
%   lowerBoundList        List of lower bound (Default = 0 or -vMax)
%   upperBoundList        List of upper bound (Default = vMax)
%   subSystemList         List of subsystem (Default = '')
%   grRuleList            List of gene-reaction rule in boolean format (and/or allowed)
%                         (Default = '');
%   geneNameList          List of gene names (used only for translation
%                         from common gene names to systematic gene names)
%   systNameList          List of systematic names
%  
                        
rxnNameList = vertcat(ExportmodelwithOpenBounds.rxnNames, IntegratedSino_Export.rxnNames);                        
rxnAbrList = vertcat(ExportmodelwithOpenBounds.rxns, IntegratedSino_Export.rxns);
rxnList = vertcat(ExportmodelwithOpenBounds.rxnsformula, IntegratedSino_Export.rxnsformula);
revFlagList = vertcat(ExportmodelwithOpenBounds.rev, IntegratedSino_Export.rev);
lowerBoundList = vertcat(ExportmodelwithOpenBounds.lb, IntegratedSino_Export.lb);
upperBoundList= vertcat(ExportmodelwithOpenBounds.ub, IntegratedSino_Export.ub);
subSystemList= vertcat(ExportmodelwithOpenBounds.subSystems, IntegratedSino_Export.subSystems);
grRuleList = vertcat(ExportmodelwithOpenBounds.grRules, IntegratedSino_Export.grRules);
geneNameList = vertcat(ExportmodelwithOpenBounds.genes, IntegratedSino_Export.genes);

fprintf('\n\n 6. Combining models');
PSmodelOriginal = createModel(rxnAbrList,rxnNameList,rxnList,revFlagList, lowerBoundList, upperBoundList, subSystemList, grRuleList);
PSmodelOriginal.lb = transpose(PSmodelOriginal.lb);
PSmodelOriginal.ub = transpose(PSmodelOriginal.ub);
PSmodelOriginal.c = transpose(PSmodelOriginal.c);

PSmodelOriginal = changeObjective(PSmodelOriginal, 'BiomassRoot');
sol = optimizeCbModel(PSmodelOriginal, 'max');



% flux in ex reactions
sol.x(findRxnIDs(PSmodelOriginal, PSmodelOriginal.rxns(findExcRxns(PSmodelOriginal))));

PSmodelOriginal = changeObjective(PSmodelOriginal, 'rxn06874');
sol = optimizeCbModel(PSmodelOriginal);
sol.x(findRxnIDs(PSmodelOriginal, PSmodelOriginal.rxns(findExcRxns(PSmodelOriginal))));


%% Identify cross-talk reactions between sino and trunca (in original model)

CrossTalk_in_PSmodel = PSmodelOriginal.rxns(~cellfun(@isempty, regexp(PSmodelOriginal.rxns,'EXCT_')));
%printRxnFormula(PSmodelOriginal, CrossTalk_in_PSmodel);



%% change name to sino-derived mets

PSmodelOriginalMetsChange = PSmodelOriginal;

 for j=1:length(PSmodelOriginalMetsChange.mets)
    
        
      metIndex = findMetIDs(PSmodelOriginalMetsChange, PSmodelOriginalMetsChange.mets(j));
      
      NewMetName = strrep(PSmodelOriginalMetsChange.mets(metIndex), '_c0[c]', '[c]');
      NewMetName = strrep(NewMetName, '_e0[c]', '[e]');

           
      PSmodelOriginalMetsChange.mets(metIndex) = NewMetName;
 end
    


