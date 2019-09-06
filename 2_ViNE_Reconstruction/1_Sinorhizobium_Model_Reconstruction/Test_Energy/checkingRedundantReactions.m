%% Check identical reactions
for n = 1:length(model.rxns)
    for m = 1:length(model.rxns)
        total = sum(model.S(:,n) == model.S(:,m));
        if total == length(model.mets)
            if n ~= m
                n
                m
            end
        end
    end
end

%% Check products and reactants switched
for n = 1:length(model.rxns)
    for m = 1:length(model.rxns)
        total = sum(model.S(:,n) == -1 * model.S(:,m));
        if total == length(model.mets)
            if n ~= m
                n
                m
            end
        end
    end
end

