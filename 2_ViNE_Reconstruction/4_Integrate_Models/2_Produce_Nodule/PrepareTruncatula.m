%% Prepare truncatula model 2


  TruncaNewNames = importdata('TruncaNamesMapping_of_shared.NewMetNames');
  TruncaOldNames = importdata('TruncaNamesMapping_of_shared.OldMetNames');

  
   for i = 1:length(TruncaOldNames)

       metindex=findMetIDs(IntegratedTrunca, TruncaOldNames(i));
    
     % _b mets do not formally exist in the mat version of the model. the following
     % if statement avoids this issue
     
     if (metindex > 0)
        IntegratedTrunca.mets(metindex) = TruncaNewNames(i);
     end
     
   end
    
    IntegratedTrunca.rxnsformula = printRxnFormula(IntegratedTrunca, IntegratedTrunca.rxns, false);

    sol = optimizeCbModel(IntegratedTrunca);

    exchangers = ~cellfun(@isempty, regexp(IntegratedTrunca.rxns,'^TE[A-Z]|Cons'));   
    ExportmodelwithOpenBounds = changeRxnBounds(IntegratedTrunca,IntegratedTrunca.rxns(exchangers),10,'u');
    
    ExportmodelwithOpenBounds = IntegratedTrunca;
    
    ExportmodelwithOpenBounds = changeObjective(ExportmodelwithOpenBounds, 'BiomassWithOutStarch');
    sol = optimizeCbModel(ExportmodelwithOpenBounds );

    ExportmodelwithOpenBounds.rxnsformula = printRxnFormula(ExportmodelwithOpenBounds , ExportmodelwithOpenBounds.rxns, false);
    sol = optimizeCbModel(ExportmodelwithOpenBounds);

