noduleModel = PSmodelOriginalMetsChange;

% Replace truncatula met names with original
   for i = 1:length(TruncaOldNames)

       metindex=findMetIDs(noduleModel, TruncaNewNames(i));
    
     % _b mets do not formally exist in the mat version of the model. the following
     % if statement avoids this issue
     
     if (metindex > 0)
        noduleModel.mets(metindex) = TruncaOldNames(i);
     end
     
   end
   
% Replace meliloti met names with original
    for i = 1:length(SinoOldNames)
        
        SinoOldNames{i} = strrep(SinoOldNames{i}, '_b', '[b]');
        SinoOldNames{i} = strrep(SinoOldNames{i}, '_c0', '[c]');
        SinoOldNames{i} = strrep(SinoOldNames{i}, '_e0', '[e]');
        SinoNewNames{i} = strrep(SinoNewNames{i}, '_C', '[b]');
        SinoNewNames{i} = strrep(SinoNewNames{i}, '_c0', '[c]');
        SinoNewNames{i} = strrep(SinoNewNames{i}, '_e0', '[e]');

        metindex=findMetIDs(noduleModel, SinoNewNames(i));
    
     % _b mets do not formally exist in the mat version of the model. the
     % following if statement avoids this issue
     
     if (metindex > 0)
          
                 noduleModel.mets(metindex) = SinoOldNames(i);
         
     end
     
    end
