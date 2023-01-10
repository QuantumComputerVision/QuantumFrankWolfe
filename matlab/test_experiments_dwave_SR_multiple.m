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
dispNames = {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'matchQuantum (ours)'};

figure;
t = tiledlayout(2,2);
t.Padding = 'none';
t.TileSpacing = 'compact';
h1=nexttile; beautify_plot;
h2=nexttile; beautify_plot;
h3=nexttile; beautify_plot;
h4=nexttile; beautify_plot;


%% test VAR_C_SR
%  N=4, NC=4, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3:
folderDWave = './output/dwave/C1_SR_VAR_ALL';
fileName = './output/dwave/02_VAR_SR_1000.txt';
SRs = [0 0.05 0.1 0.15 0.2 0.25];

N = 4; NC = 4;
% [allPossibleSolutions, allPossibleSolutionsPerm] = get_all_solution_space(N, NC);

xValues = [];
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

folderNameSynth = ['./dataset/synth_data_multi_random_exp'];
folderNameDWaveOut = ['./output/dwave/C1_SR_VAR_ALL'];

for i=1:length(SRs)
    swapRatio = SRs(i);
    PRtmp = [];
    PRConstmp = [];
    accuraciestmp = [];
    accuraciesConstmp = [];
    accuraciesDWavetmp = [];
    energiestmp = [];
    energiesDWavetmp = [];
    
    for rndmExp = 1:10
        % exhaustive solution to synthetic problem
        fileNameSynth = [folderNameSynth '/Q_N4_NC4_C1_SR' num2str(swapRatio) '_L2.5_' num2str(rndmExp) '.mat'];
        data = load(fileNameSynth);
        data = data.data;
        Qcons = data.Qcons;
        s = data.s;
        Pijs = data.Pijs;
        I = data.I;
        
        [PisCons, qCons, currEnergyCons] = solve_exhaustive_constrained(Qcons, s, NC, N);
        
        folderNameDWave = folderNameDWaveOut;
        folderNameDWave = [folderNameDWave '/C1_SR' num2str(swapRatio)];
        folderNameDWave = erase(folderNameDWave,'.');
        fileNameDWave = [folderNameDWave '/Q_N4_NC4_C1_SR' num2str(swapRatio) '_L2.5_' num2str(rndmExp)];
        fileNameDWave = ['.' erase(fileNameDWave,'.')  '.txt'];
        [PisCells, energies, qubits] = load_dwave_text_output(fileNameDWave, N, NC);
        PisDWave = PisCells(:,1);
        qDWave = perms_cell_to_q(PisDWave);
        
        currEnergyDWave = -qDWave'*Qcons*qDWave - s'*qDWave;
        
        % run state of the art        
        [PisMatchEIG, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchEIG'); % MatchEIG
        [PisMatchALS, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchALS'); % MatchALS
        [PisMatchLift, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchLift'); % MatchLift
        [PisBirdal, ~] = perm_sync_sota(Pijs, I, N, NC, 'matchBirdal'); % Birdal'19
        
        if (~isempty(allPossibleSolutions))
            [PisExh, qConsExh] = solve_exhaustive_unconstrained(Qcons, s, NC, N, allPossibleSolutions);
            [~, recall, ~, accExh] = get_synth_solution_PR(PisExh);
        end
        
        [precisionCons, recallCons, ~, accCons] = get_synth_solution_PR(PisCons);
        [precision, recall, ~, accDWave] = get_synth_solution_PR(PisDWave);
        [precision, recall, ~, accMatchEIG] = get_synth_solution_PR(PisMatchEIG);
        [precision, recall, ~, accMatchALS] = get_synth_solution_PR(PisMatchALS);
        [precision, recall, ~, accMatchLift] = get_synth_solution_PR(PisMatchLift);
        [precision, recall, ~, accBirdal] = get_synth_solution_PR(PisBirdal);
        
        energiesDWavetmp = [energiesDWavetmp, currEnergyDWave];
        energiestmp = [energiestmp, currEnergyCons];
        % {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'd-wave (ours)'};
        accuraciestmp = [accuraciestmp, [accCons, accMatchEIG, accMatchALS, accMatchLift, accBirdal, accDWave]'];
        PRtmp = [PRtmp [precision; recall]];
        PRConstmp = [PRConstmp [precisionCons; recallCons]];
        
    end
    
    energiesDWave = [energiesDWave, mean(energiesDWavetmp)];
    energies = [energies, mean(energiestmp)];
    % {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'd-wave (ours)'};
    accuracies = [accuracies, mean(accuraciestmp,2)];
    PR = [PR mean(PRtmp,2)];
    PRCons = [PRCons mean(PRConstmp,2)];
    
    xValues = [xValues swapRatio];
end

fclose(fid);

set(gcf, 'currentaxes', h1); %subplot(2,2,1),
plot_multiple(xValues, accuracies, clrs, markers, dispNames, 1);
title('(a) Accuracy vs. noise', 'FontSize', 12);

