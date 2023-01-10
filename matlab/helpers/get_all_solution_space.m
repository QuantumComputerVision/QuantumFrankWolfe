function [allPossibleSolutions, allPossibleSolutionsPerm] = get_all_solution_space(N, nCams)

permIndicesP = uint16(perms(1:N));
numPerms = size(permIndicesP,1);
allPossibleSolutionsPerm =nmultichoosek(1:numPerms, nCams);

if (N<4)
    allPossibleSolutions = makebits_augmented(N*N*(nCams-1), N);
else
    allPossibleSolutions = [];
end

end