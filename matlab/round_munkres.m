function [Pis, costs] = round_munkres(Pis)

% re-arrange into solutions
% costMat = u*u';
costs = zeros(length(Pis),1);
for i=1:length(Pis)
    costMat = 1 - Pis{i};
    [assignment, cost] = munkres(costMat);
    Pis{i} = assignment;
    costs(i) = cost;
end

end