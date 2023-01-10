clear all;
clc;
close all;

addpath('./permutations/');
addpath('./dataset/');
addpath('./helpers/');
addpath('./sota/');

clrs = lines(7);
markers = {'o', '+', '*', 'x', 'square', 'diamond', 'pentagram', 'hexagram', 'v'};
allPossibleSolutions = [];

%% test VAR_C_SR 
%  N=4, NC=4, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3:
fileName = './output/dwave/04_VAR_C_SR.txt';

N = 4; NC = 4;
% [allPossibleSolutions, allPossibleSolutionsPerm] = get_all_solution_space(N, NC);

PR = [];
PRCons = [];
accuracies = [];
accuraciesCons = [];
accuraciesDWave = [];
energies = [];
energiesDWave = [];

fid = fopen(fileName, 'r');
qs = [];
prevCompl = -1;
i=0;
while (~feof(fid))
    SR = str2double(fgetl(fid));
    C = str2double(fgetl(fid));
    occ = str2double(fgetl(fid));
    qubits = str2double(fgetl(fid));
    chain = str2double(fgetl(fid));
    sol = (fgetl(fid));
    if (prevCompl==C)
        continue;
    else
        prevCompl = C;
    end
    i
    i=i+1;
    sol = erase(sol, "]");
    sol = erase(sol, "[");
    qDWave = str2num(sol);
    qDWave = qDWave(:);
    qs = [qs; qDWave];
    
    fileNameSynth = ['./dataset/synth_data/Q_N4_NC4_C' num2str(C) '_SR' num2str(SR) '_L2.5.mat'];
    data = load(fileNameSynth);
    data = data.data;
    
    Qcons = data.Qcons;
    s = data.s;
    currEnergyDWave = -qDWave'*Qcons*qDWave - s'*qDWave;
    
    [PisCons, qCons, currEnergyCons] = solve_exhaustive_constrained(Qcons, s, NC, N);
    PisDWave = perms_q_to_cell(qDWave, N);
    
    % run state of the art
    Pijs = data.Pijs;
    I = data.I;
    [PisMatchEIG, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchEIG'); % MatchEIG
    [PisMatchALS, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchALS'); % MatchALS
    [PisMatchLift, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchLift'); % MatchLift
    [PisBirdal, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchBirdal'); % Birdal'19
    
    if (~isempty(allPossibleSolutions))
        [PisExh, qConsExh] = solve_exhaustive_unconstrained(Qcons, s, NC, N, allPossibleSolutions);
        [precision, recall, ~, accExh] = get_synth_solution_PR(PisExh);
    end
    
    [precisionCons, recallCons, ~, accCons] = get_synth_solution_PR(PisCons);
    [precision, recall, ~, accDWave] = get_synth_solution_PR(PisDWave);
    [precision, recall, ~, accMatchEIG] = get_synth_solution_PR(PisMatchEIG);
    [precision, recall, ~, accMatchALS] = get_synth_solution_PR(PisMatchALS);
    [precision, recall, ~, accMatchLift] = get_synth_solution_PR(PisMatchLift);
    [precision, recall, ~, accBirdal] = get_synth_solution_PR(PisBirdal);
    
    energiesDWave = [energiesDWave, currEnergyDWave];
    energies = [energies, currEnergyCons];
    accuracies = [accuracies, [accCons, accDWave, accMatchEIG, accMatchLift, accMatchALS, accBirdal]'];
    PR = [PR [precision; recall]];
    PRCons = [PRCons [precisionCons; recallCons]];
         
end

fclose(fid);

figure;
t = tiledlayout(2,2);
t.Padding = 'none';
t.TileSpacing = 'compact';
h1=nexttile; beautify_plot;
h2=nexttile; beautify_plot;
h3=nexttile; beautify_plot;
h4=nexttile; beautify_plot;
set(gcf, 'currentaxes', h1); %subplot(2,2,1),
hold on; plotHandles(1) = plot(lambdas, energies, 'LineWidth', 2, 'Marker','o', 'Color', clrs(1,:));
hold on; plotHandles(2) = plot(lambdas, energiesDWave, 'LineWidth', 2, 'Color', clrs(2,:));

set(gcf, 'currentaxes', h2);  % subplot(2,2,2);
%hold on, plot(lambdas, abs(energies-energiesDWave), 'LineWidth', 2, 'Color', clrs(experiment,:), 'DisplayName', 'difference in energy levels');

hold on, plot(lambdas, accuraciesCons, 'LineWidth', 2, 'Color', clrs(1,:), 'DisplayName', 'exhaustive', 'Marker', markers{1});
% hold on, plot(PR(1,:), PR(2,:), 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)]);

% set(gcf, 'currentaxes', h4); % subplot(2,2,4);
hold on, plot(lambdas, accuraciesDWave, 'LineWidth', 2, 'Color', clrs(experiment,:), 'DisplayName', 'dwave', 'Marker', markers{2});
% hold on, plot(PRCons(1,:), PRCons(2,:), 'LineWidth', 1, 'Color', clrs(experiment,:), 'DisplayName', ['Exp - ', num2str(experiment)]);

set(gcf, 'currentaxes', h3); %subplot(2,2,3);

return;

%% test 01_VAR_N_NC 
% SR = 0, C = 1, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3
fileName = './output/dwave/01_VAR_N_NC.txt';

PR = [];
PRCons = [];
accuraciesCons = [];
accuraciesDWave = [];
energies = [];
energiesDWave = [];

fid = fopen(fileName, 'r');
qs = [];
while (~feof(fid))
    SR = str2double(fgetl(fid));
    C = str2double(fgetl(fid));
    occ = str2double(fgetl(fid));
    qubits = str2double(fgetl(fid));
    chain = str2double(fgetl(fid));
    sol = (fgetl(fid));
    sol = erase(sol, "]");
    sol = erase(sol, "[");
    qDWave = str2num(sol);
    qs = [qs; qDWave];
    
    fileNameSynth = ['./dataset/synth_data/Q_N4_NC4_C' num2str(C) '_SR' num2str(SR) '_L2.5.mat'];
    data = load(fileNameSynth);
    data = data.data;
    
    Qcons = data.Qcons;
    s = data.s;
    currEnergyDWave = -qDWave'*Qcons*qDWave - s'*qDWave;
    
    [PisCons, qCons, currEnergyCons] = solve_exhaustive_constrained(Qcons, s, NC, N);
    PisDWave = perms_q_to_cell(qDWave, N);
    if (~isempty(allPossibleSolutions))
        [PisExh, qConsExh] = solve_exhaustive_unconstrained(Qcons, s, NC, N, allPossibleSolutions);
        [precision, recall, ~, accExh] = get_synth_solution_PR(PisExh);
    end
    
    [precisionCons, recallCons, ~, accCons] = get_synth_solution_PR(PisCons);
    [precision, recall, ~, accDWave] = get_synth_solution_PR(PisDWave);
    
    energiesDWave = [energiesDWave, currEnergyDWave];
    energies = [energies, currEnergyCons];
    accuraciesDWave = [accuraciesDWave, accDWave];
    accuraciesCons = [accuraciesCons, accCons];
    PR = [PR [precision; recall]];
    PRCons = [PRCons [precisionCons; recallCons]];
         
end

fclose(fid);

return;


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