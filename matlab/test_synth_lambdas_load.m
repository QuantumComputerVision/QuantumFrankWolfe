clear all;
clc;
close all;

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');

fileName = './dataset/synth2d_dataN3NC3';
fileNameQ = './dataset/QMatrixN3NC3';
fileNameAllData = './dataset/synth_data_outputN3NC3.mat';
data = load(fileNameAllData);
data = data.data;

QGnd = data.QGnd;
Q = data.Q;
Qcons = data.Qcons;
s = data.s;
PijsGnd = data.PijsGnd;
Pijs = data.Pijs;
Pis = data.Pis;
N = data.N; % # of points
G = data.G; % # of points
I = data.I; % # of points
Iupper = data.Iupper; % # of points
nCams = data.nCams;
lambda = data.lambda;
swapRatio = data.swapRatio;
completeness = data.completeness;

[Qcons0, s0] = compute_Q_constrained(PijsGnd, I, nCams, 0);

% 18 bits preceeded by 3x3 identity matrices
allPossibleSolutions = makebits_augmented(N*N*(nCams-1), N);

% debug lambda now
%lambdas = [0, 0.01, 0.015 0.025, 0.04 0.05, 0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.35 0.4 0.425 0.45 0.5 0.525 0.55 0.6 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 10 11]; 
lambdas = [0, 0.01, 0.015 0.025, 0.04 0.05, 0.075 0.1 0.125 0.15 0.2 0.25 0.3 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.4 0.425 0.45 0.5 0.525 0.55 0.6 0.75 1 1.5 2 2.5 3 3.5 4 4.5 5]; 

[PisCons, qCons] = solve_exhaustive_constrained(Qcons0, s0, nCams, N); % solves only permutations
currEnergy = -qCons'*Qcons0*qCons - s0'*qCons;

energies = [];
energies0 = [];
differences = [];
minDiffEnergy = 9999;
bestInd = 1;
lambdaSol = lambdas(bestInd);
bestq = zeros(length(qCons),1);
bestPis = PisCons;
for i=1:length(lambdas)
    lambda = lambdas(i);
    [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
    % [Qcons, s] = compute_Q_constrained(PijsGnd, I, nCams, lambda);
    [PisCons, qCons] = solve_exhaustive_constrained(Qcons, s, nCams, N);
    currEnergy = -qCons'*Qcons*qCons - s'*qCons;
    
    [PisExh, qConsExh] = solve_exhaustive_unconstrained(Qcons, s, nCams, N, allPossibleSolutions);

    currEnergyExh = -qConsExh'*Qcons*qConsExh - s'*qConsExh;
        
    diffEnergy = abs(currEnergyExh-currEnergy);
    disp(['lambda: ' num2str(lambda) ', diff: ', num2str(diffEnergy)]);
    
    if (diffEnergy<minDiffEnergy)
        lambdaSol = lambda;
        bestInd = i;
        bestq = qConsExh;
        bestPis = PisExh;
    end
    energies0 = [energies0, currEnergy];
    energies = [energies, currEnergyExh];
    differences = [differences, diffEnergy];
end

%gndEnergy = ;
figure, plot(lambdas, energies); 
hold on, plot(lambdas, energies0); beautify_plot; 
%hold on, yline(currEnergy);

figure, plot(lambdas, differences); beautify_plot; 

% qEye = eye(N);
% qSol = [qEye(:)' qEye(:)' qEye(:)' ]';