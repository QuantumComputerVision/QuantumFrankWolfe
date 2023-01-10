function [M, isConnected] = pMatch2permNC(pMatch, N, NC, firstIndex)

nFeature = zeros(NC,1);
filename = cell(NC,1);
for i = firstIndex:NC
    for j = i+1:NC
        if ~isempty(pMatch(i,j).nFeature)
            nFeature(i) = N;
            nFeature(j) = N;
        end
        if ~isempty(pMatch(i,j).filename)
            filename(i) = pMatch(i,j).filename(1);
            filename(j) = pMatch(i,j).filename(2);
        end
    end
end
cumIndex = cumsum([0; nFeature]);

nFeatConst = pMatch(1,2).X;
nFeatConst = sqrt(length(nFeatConst));
xInd = [];
for i=1:N
    xInd = [xInd (i-1)*10+1 : (i-1)*10 + N];
end

M = sparse(cumIndex(end),cumIndex(end));
isConnected=1;
for i = firstIndex:NC
    for j = i+1:NC
        if ~isempty(pMatch(i,j).matchInfo)
            X = pMatch(i,j).X(xInd);
            mInd = 1:length(X);
            if (sum(X(:))<N)
                isConnected = 0;
            end
            matchList = double(pMatch(i,j).matchInfo.match);
            M(cumIndex(i)+1:cumIndex(i+1),cumIndex(j)+1:cumIndex(j+1)) = ...
                sparse(matchList(1,mInd),matchList(2,mInd),X, nFeature(i), nFeature(j));
        end
    end
end
M = M + M';