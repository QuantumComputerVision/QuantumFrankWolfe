addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');
addpath('./algorithms/');

close all;

rng(88); % so that we always get the same data

N = 2; % the number of the grid points. permutations are NxN
nCams = 3; % number of cameras (# nodes)
completeness = 1.0; % 1 is fully connected, 0 has no connectivity
swapRatio = 0; % like to add some noise? 0.01-0.15
lambda = 1.25; % regularizer (permutation constratins)

% set up the synth. problem
% Pijs are noisy relative permutations
% PijsGnd are non-noisy (true) relative permutations
% Pis are the solution whereas Pijs are the data.
[Pis, Pijs, PijsGnd, G, I, Iupper] = gen_synth_multi_graph_2d_no_grid(N, nCams, completeness, swapRatio);
% [Pis, Pijs, PijsGnd, G, I, Iupper, points] = gen_synth_multi_graph_2d(N, nCams, completeness, swapRatio);

% compute the unconstrained Q matrix from ground truth relative permutations
% without permutation constrints (unconstrained)
% QGnd induces the quadratic energy: -q'*QGnd*q where qs are vectorized.
QGnd = compute_Q(PijsGnd, I, nCams);

% verify that QGnd is maximized at q=vec(pis)
q = perms_cell_to_q(Pis);
-q'*QGnd*q

% noisy Q should lead to a higher cost
Q = compute_Q(Pijs, I, nCams);
-q'*Q*q

% now let's play with the constrained version
[Ai, ~] = compute_constraint_system(N);
A = kron(eye(nCams), Ai);
b = ones(size(A,1), 1);
% let's verify the constraints:
disp(['Constraint debugging: sum(A*q-b) = ', num2str(sum(A*q-b))]);

[Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
% [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
q = perms_cell_to_q(Pis);
dataErr = -q'*Qcons*q - s'*q;
consErr = sum(A*q-b);
disp(['Data error of regularized solver (lambda=' num2str(lambda) ') = ', num2str(dataErr)]);
disp(['Constraint error of regularized (lambda=' num2str(lambda) ') = ', num2str(consErr)]);

data = store_synth_data([], Q, Q, QGnd, Qcons, s, PijsGnd, Pijs, Pis, N, G, I, Iupper, nCams, lambda, swapRatio, completeness);

% % frank-wolfe (Alp's edits)
Aoperator = @(x) (A*x).^2;
ATranspose = @(y)  A' * bsxfun(@times, A, y);
bsquared = b.^2; % this shouldn't matter if b is always 0 or 1
maxit = 2e2;
beta0 = 10;

folder = '../datasets/synth';
[XFW] = HCGM_with_data(folder,data,Aoperator,ATranspose,bsquared,beta0,maxit);

[PisConsExh, qConsExh] = solve_exhaustive_constrained(Qcons, s, nCams, N);
-qConsExh'*Qcons*qConsExh - s'*qConsExh

figure, spy(QGnd);
figure, spy(Q);


%%%% IF YOU WANT TO SAVE %%%%


