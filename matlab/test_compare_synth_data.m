clear all;
rng(0,'twister');

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');
addpath('./algorithms/');

close all;

N = 3; % the number of the grid points. permutations are NxN
nCams = 3; % number of cameras (# nodes)
completeness = 1.0; % 1 is fully connected, 0 has no connectivity
swapRatio = 0.0; % like to add some noise? 0.01-0.15
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

% compute the exhaustive solution - this should align with the gnd-truth
[PisExhaustive, qResultExhaustive] = solve_exhaustive(Q, nCams, N);
-qResultExhaustive'*Q*qResultExhaustive
XResultExhasutive = qResultExhaustive*qResultExhaustive';

% now let's play with the constrained version
[Ai, ~] = compute_constraint_system(N);
A = kron(eye(nCams), Ai);
b = ones(size(A,1), 1);
% let's verify the constraints:
disp(['Cosntraint debugging: sum(A*q-b) = ', num2str(sum(A*q-b))]);

%bsquared = b.^2; % this shouldn't matter if b is always 0 or 1
maxit = 100;
beta0 = 1;
verbose = false;
dWaveMode = 2;
isSave = 1;
isVisualize = 0;
tic;
[PisQFWAL, XFW, costMunkres, constraintResidual] = CoCoFW(Q, N, A, b, 'dwaveFWAL', maxit, beta0, verbose, dWaveMode, isSave, isVisualize);
timeCoco = toc;

qsol = perms_cell_to_q(PisQFWAL);    
costCoco = -qsol'*Q*qsol;
costExh = -qResultExhaustive'*Q*qResultExhaustive;
cocoObj = XFW(:)'*Q(:);
disp(['cost CoCo: ' num2str(costCoco) ', CoCo-obj: ' num2str(cocoObj) ', cost exh: ' num2str(costExh)]);

%%
% saveDir = '../output/noisy_synchronization/Dwave_2_20220307T171857/';
% saveDir = '../output/noisy_synchronization/20220307T003746/';
% saveDir = 'results/Dwave_1_20220307T195635/';
saveDirs = {'../output/noiseless_synchronization/20220307T083349_noiseless/', 'results/Dwave_1_20220307T210300/'};

%saveDir = '../output/noiseless_synchronization/20220307T083349_noiseless/';

[fig, layout, hs] = open_multiplot(2, 10);

allAcc=[];
for k=1:2
    saveDir = saveDirs{k};
    j=1;
    for t = 1:maxit
        if (mod(t,10)==0)
            dataFWAL = load([saveDir,num2str(t)]);
            X = dataFWAL.X;
            Q = dataFWAL.grad(2:end,2:end);
            set(gcf, 'currentaxes', hs(k,j));
            imagesc(Q,[0,17]);
            if(j==10)
                colorbar;
            end
            if (j==1 && k==1)
                ylabel('FWAL');
            end
            if (j==1 && k==2)
                ylabel('Q-FWAL');
            end

            % [qx, energyMunkres,u,PsX] = project_relaxed(dataFWAL.grad, X);
            % PsX = perms_transform_first_to_eye(PsX);
            % [~, ~, ~, accQFWAL] = get_synth_solution_PR(PsX);
            % accQFWAL
            % allAcc=[allAcc, accQFWAL];
            j=j+1;
        end
    end
end
sgtitle('iterations (t)');
% mean(allAcc)
%%
% saveDir = '../output/noisy_synchronization/20220307T003746/';
% saveDir = '../output/noisy_synchronization/Dwave_2_20220307T171857/';
%saveDirs = 'results/Dwave_1_20220307T193912/';
%saveDir = 'results/Dwave_1_20220307T200953/';
%saveDir = 'results/Dwave_2_20220307T204616/';
%saveDir = 'results/Dwave_2_20220307T204822/';
%saveDir = '../output/noiseless_synchronization/20220307T083349_noiseless/';
saveDir = 'results/Dwave_2_20220307T205956/';
%saveDir = 'results/Dwave_1_20220307T210300/';

% saveDirs = {'../output/noiseless_synchronization/20220307T083349_noiseless/', 'results/Dwave_1_20220307T210300/'};
%saveDirs = {'../output/noiseless_synchronization/20220307T083349_noiseless/', 'results/Dwave_1_20220307T210300/', 'results/Dwave_2_20220307T212632/','../output/noisy_synchronization/20220307T003746/'};

saveDirs = {'../output/noiseless_synchronization/20220307T083349_noiseless/', 'results/Dwave_1_20220307T210300/', 'results/Dwave_2_20220307T212632/','../output/noisy_synchronization/20220307T003746/'};
Xgroundtruth = q*q';
maxit=100;
figure;
set(gcf, 'Color', [1 1 1]); % Sets figure background
set(gca, 'Color', [1 1 1]); % Sets axes background
clrs = lines(7);
markers = {'diamond', '+', 'v', '.'};
for i=1:length(saveDirs)
    saveDir = saveDirs{i};
    errors=zeros(maxit,1);
    for t = 1:maxit
        %if (mod(t,10)==0)
        dataFWAL = load([saveDir,num2str(t)]);
        Xt = dataFWAL.X(2:end,2:end);
        err = abs(Q(:)' * Xt(:) - Q(:)' * Xgroundtruth(:));
        errors(t)=err;
        %end
    end
    hold on, semilogy(1:maxit, errors,'LineWidth',3,'Marker',markers{i}, 'Color', clrs(i,:));
    if(i==1)
        ax = gca;
        ax.YMinorGrid = 'on';
        set(gca,'yscale','log')
        grid on;
    end
    drawnow;
end
%grid on;
legend({'FWAL','Q-FWAL','Q-FWAL-SA','FWAL (\sigma=0.2)'});
axis([6, 100, 0.01, 60]);
xlabel('iterations (time)','interpreter','latex');
ylabel('$\log(\varepsilon_\mathrm{CC})$','interpreter','latex');
return;

[Qcons, s] = compute_Q_constrained(Pijs, I, nCams, lambda);
% [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
q = perms_cell_to_q(Pis);
-q'*Qcons*q - s'*q

[PisConsExh, qConsExh] = solve_exhaustive_constrained(Qcons, s, nCams, N);
-qConsExh'*Qcons*qConsExh - s'*qConsExh

figure, spy(QGnd);
figure, spy(Q);
