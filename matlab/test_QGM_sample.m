
addpath permutations/
addpath algorithms/

data = load('./dataset/QGM/N3.mat'); % 9x9
% N4= load('./dataset/QGM/N4.mat'); % multiple 16x16 graph matching problems

numProblems = size(data.W1, 1);

% all problems stacked 
Qs = data.W1;
cs = data.c1;

i = 1; % I will only consider the first problem (loop if you wanna do all)
Q = squeeze(Qs(i, :, :)) + diag(squeeze(cs(i, :)));
Q = 0.5*(Q+Q'); % note that CoCoFW was using -Q.

% at this stage Q includes both W and c corresponding to the GM problem
% such that E = x'*Q*x (equivalently E = x'*Qs_i*x + cs_i' *x)
figure, imagesc(Q);

% now you want to solve this with respect to some constraints
% what are those? they are the permutation-ness constraints.
% the way these (Ax=b) type constraints are obtained is given in:
% Quantum Permutation Synchronization: https://arxiv.org/abs/2101.07755
% (Eq. 14)
N = int32(fix(sqrt(length(data.W1(1,:,:)))));
[A, b] = compute_constraint_system(N);

% let's try an exhaustive solution first:
qSol = solve_exhaustive_QGM(Q, N);
qSolExhaustive = qSol{1};

% now let's attempt to use QFW to solve this:
maxit = 20;
beta0 = 1;
verbose = 1;
isSave = 0;
isVisualize = 1;
[Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'FWAL', maxit, beta0, verbose, 0, isSave, isVisualize);

% let's retrieve the solution and compare with the ground truth:
disp(['The error is: ' num2str(norm(Pis{1}-qSolExhaustive))]);
disp(['The final cost is: ' num2str(costMunkres), ' and the constraint cost is: ' num2str(constraintResidual(end))]);