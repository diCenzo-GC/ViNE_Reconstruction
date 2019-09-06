function output = paretoFrontier(model, rxn1, rxn2, steps)

%% Set defaults

% Set steps default
if nargin < 4
    steps = 100;
elseif isempty(steps)
    steps = 100;
end

%% Find the reaction optimums

% Optimize for the first reaction
model = changeObjective(model, rxn1);
sol = optimizeCbModel(model);
rxn1_max = sol.f;

% Optimize for the second reaction
model = changeObjective(model, rxn2);
sol = optimizeCbModel(model);
rxn2_max = sol.f;

%% Calculate the pareto optimum

% Set output variable
output = cell(steps * 2 + 2, 2);

% Get reaction positions
pos1 = findRxnIDs(model, rxn1);
pos2 = findRxnIDs(model, rxn2);

% Vary the first reaction
model2 = changeObjective(model, rxn2);
for n = 0:steps
    if n == 0
        model2.lb(pos1) = 0;
        model2.ub(pos1) = 0;
    else
        model2.lb(pos1) = n * (rxn1_max / steps);
        model2.ub(pos1) = n * (rxn1_max / steps);
    end
    sol = optimizeCbModel(model2);
    if isempty(sol.x)
        output{n+1,1} = 0;
        output{n+1,2} = 0;
    else
        output{n+1,1} = sol.x(pos1);
        output{n+1,2} = sol.x(pos2);
    end
end

% Vary the second reaction
model2 = changeObjective(model, rxn1);
for n = 0:steps
    if n == 0
        model2.lb(pos2) = 0;
        model2.ub(pos2) = 0;
    else
        model2.lb(pos2) = n * (rxn2_max / steps);
        model2.ub(pos2) = n * (rxn2_max / steps);
    end
    sol = optimizeCbModel(model2);
    if isempty(sol.x)
        output{n+1,1} = 0;
        output{n+1,2} = 0;
    else
        output{n+2+steps,1} = sol.x(pos1);
        output{n+2+steps,2} = sol.x(pos2);
    end
end

% Add headers
output = vertcat({rxn1, rxn2}, output);

end