
% Generates a synthetic camera graph with uniformly sampled relative
% permutations. (potentially noisy)
% N            : number of nodes or view-points or cameras
% completeness : fraction [0,1] of pairwise edges to keep
% swapRatio    : relative number of swaps to be made per relative
% (introduces wrong matches into the pairwise estimates)
% permutation (this adds noise)
% returns a sparse adjacency matrix G
% Pi: Absolute permutations
% Pij: Relative permutations (noisy)
% PijGnd: Ground truth relative permutations
% G: sparse graph adj. mat.
% I: edges of the graph [ [1 2]', [3 4]'...]
% Iupper: only edges of the upper triangular G
function [Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio)

%% generate the pose graph edges
Npair = nCams*(nCams-1)/2;
I = zeros(2, Npair);

% generate all pairwise edges
k=1;
for i=1:nCams-1
    for j=i+1:nCams
        I(:,k)=[i; j]; k=k+1;
    end
end

% now keep a portion of the edges
e = ceil(completeness*Npair);
ind = randperm(Npair, e);
I = I(:, ind);
vals = ones(1, e);
G = sparse(I(1,:), I(2,:), vals, nCams, nCams);
e = length(I);

%% synthesize absolute permutations and data
Pis = {};
for i =1:nCams
    
    % currently the absolute permutations are always identity.
    %Pi = perm_rand(N);
    Pi = eye(N);
    status = perm_check(Pi);
        
    % TODO: very rarely perm_rand can generate a partial perm due to 
    % munkres. we can fix it in the future maybe.
%     while (~status) % some
%         Pi = perm_rand(N);
%         status = perm_check(Pi);
%         if (status~=1)
%             disp('something was wrong');
%         end
%     end
    Pis = [Pis Pi];
end

%% compute relative permutations with noise
Pijs1 = {};
PijsGnd1 = {};
Pijs2 = {};
PijsGnd2 = {};

p = randperm(e*N);
maxSwap = uint32(round(swapRatio*double(length(p))));
isSwap = zeros(length(p),1);
isSwap(p(1:maxSwap)) = 1;

for k=1:e
    i = I(1, k); j = I(2, k);
    Pi = Pis{i};
    Pj = Pis{j};
    Pijgnd = perm_relative(Pi, Pj);
    Pij = Pijgnd;
    
    swapInd = isSwap((k-1)*N+1:k*N);
    for swap = 1:N
        if(swapInd(swap))
            i = randi(N);
            j = randi(N);
            Pij = swap_rows(Pij, i, j);
        end
    end
    
    % randomly generate wrong matches
%     numSwaps = uint32(swapRatio*N);
%     for swap = 1:numSwaps
%         % swap two random indices
%         i = randi(N);
%         j = randi(N);
%         Pij = swap_rows(Pij, i, j);
%     end
    
    Pijs1 = [Pijs1 Pij];
    Pijs2 = [Pijs2 Pij'];
    PijsGnd1 = [PijsGnd1 Pijgnd];
    PijsGnd2 = [PijsGnd2 Pijgnd'];
end

Pijs = [Pijs1 Pijs2];
PijsGnd = [PijsGnd1 PijsGnd2];

G = G+G'; % symmetrize
% [I1, I2] = find(G);
%I = [I'; [I(]'];
Iupper = I;
I = [I [I(2,:) ; I(1,:)]];

end
