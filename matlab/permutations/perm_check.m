% check whether a matrix is a permutation
function [result] = perm_check(P)

[M, N] = size(P);
isSquare = M==N;

if (~isSquare)
    result = 0;
    return;
end

allPos = all(P(:)>=0);

if (~allPos)
    result = 0;
    return;
end

rowsSumToOne = all(P*ones(N,1)==1);

if (~rowsSumToOne)
    result = 0;
    return;
end

colsSumToOne = all(P'*ones(N,1)==1);

if (~colsSumToOne)
    result = 0;
    return;
end

result=1;

end