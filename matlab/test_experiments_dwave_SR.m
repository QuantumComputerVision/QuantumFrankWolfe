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
if(1)
%% test VAR_C_SR 
%  N=4, NC=4, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3:
fileName = './output/dwave/02_VAR_SR_1000.txt';

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
i=0;
while (~feof(fid))
    SR = str2double(fgetl(fid));
    occ = str2double(fgetl(fid));
    qubits = str2double(fgetl(fid));
    chain = str2double(fgetl(fid));
    sol = (fgetl(fid));
    if (prevCompl==SR)
        continue;
    else
        prevCompl = SR;
    end
    i
    i=i+1;
    sol = erase(sol, "]");
    sol = erase(sol, "[");
    qDWave = str2num(sol);
    qDWave = qDWave(:);
    qs = [qs; qDWave];
    
    fileNameSynth = ['./dataset/synth_data/Q_N4_NC4_C1_SR' num2str(SR) '_L2.5.mat'];
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
    % {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'd-wave (ours)'};
    accuracies = [accuracies, [accCons, accMatchEIG, accMatchLift, accMatchALS, accBirdal, accDWave]'];
    PR = [PR [precision; recall]];
    PRCons = [PRCons [precisionCons; recallCons]];
         
    xValues = [xValues SR];
end

fclose(fid);

set(gcf, 'currentaxes', h1); %subplot(2,2,1),
plot_multiple(xValues, accuracies, clrs, markers, dispNames, 1);
title('(a) Accuracy vs. noise', 'FontSize', 12);


%% N and NC
%  N=4, NC=4, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3:
fileName = './output/dwave/01_VAR_N_NC.txt';

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
prevN = -1; prevNC = -1;
i=0;
while (~feof(fid))
    N = str2double(fgetl(fid));
    NC = str2double(fgetl(fid));
    occ = str2double(fgetl(fid));
    qubits = str2double(fgetl(fid));
    chain = str2double(fgetl(fid));
    sol = (fgetl(fid));
    if (prevNC==NC && prevN==N)
        continue;
    else
        prevN=N;
        prevNC=NC;
    end
    disp([i N NC]);
    i=i+1;
    sol = erase(sol, "]");
    sol = erase(sol, "[");
    qDWave = str2num(sol);
    qDWave = qDWave(:);
    qs = [qs; qDWave];
    
    % fileNameSynth = ['./dataset/synth_data/Q_N4_NC4_C' num2str(C) '_SR0_L2.5.mat'];
    fileNameSynth = ['./dataset/synth_data/Q_N' num2str(N) '_NC' num2str(NC) '_C1_SR0_L2.5.mat'];
    data = load(fileNameSynth);
    data = data.data;
    
    Qcons = data.Qcons;
    s = data.s;
    currEnergyDWave = -qDWave'*Qcons*qDWave - s'*qDWave;
    
    %[PisCons, qCons, currEnergyCons] = solve_exhaustive_constrained(Qcons, s, NC, N);
    currEnergyCons = -999;
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
    
    %[precisionCons, recallCons, ~, accCons] = get_synth_solution_PR(PisCons);
    accCons = 1.0;
    [precision, ~, ~, accDWave] = get_synth_solution_PR(PisDWave);
    [precision, recall, ~, accMatchEIG] = get_synth_solution_PR(PisMatchEIG);
    [precision, recall, ~, accMatchALS] = get_synth_solution_PR(PisMatchALS);
    [precision, recall, ~, accMatchLift] = get_synth_solution_PR(PisMatchLift);
    [precision, recall, ~, accBirdal] = get_synth_solution_PR(PisBirdal);
    
    energiesDWave = [energiesDWave, currEnergyDWave];
    energies = [energies, currEnergyCons];
    % {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'd-wave (ours)'};
    accuracies = [accuracies, [accCons, accMatchEIG, accMatchALS, accMatchLift, accBirdal, accDWave]'];
    %PR = [PR [precision; recall]];
    %PRCons = [PRCons [precisionCons; recallCons]];
         
    xValues = [xValues N*N*(NC-1)];
end

fclose(fid);

[xSort, indSort] = sort(xValues);
accuracies = accuracies(:, indSort);

set(gcf, 'currentaxes', h2); %subplot(2,2,1),
plot_multiple(xSort, accuracies, clrs, markers, dispNames, 1);
title('(c) Accuracy vs. problem size', 'FontSize', 12);



%% completeness test
%  N=4, NC=4, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3:
fileName = './output/dwave/03_VAR_C_N4_NC4.txt';

N = 4 ; NC = 4;
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
i=0;
while (~feof(fid))
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
    
    fileNameSynth = ['./dataset/synth_data/Q_N4_NC4_C' num2str(C) '_SR0_L2.5.mat'];
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
    % {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'd-wave (ours)'};
    accuracies = [accuracies, [accCons, accMatchEIG, accMatchALS, accMatchLift, accBirdal, accDWave]'];
    PR = [PR [precision; recall]];
    PRCons = [PRCons [precisionCons; recallCons]];
         
    xValues = [xValues C];
end

fclose(fid);

set(gcf, 'currentaxes', h3); %subplot(2,2,1),
plot_multiple(xValues, accuracies, clrs, markers, dispNames, 1);
title('(b) Accuracy vs. completeness', 'FontSize', 12);

%return;

end
%% SR & C together
%  N=4, NC=4, N_samples = 200, lambda = 2.5, CHAIN STRENGTH = 3:
fileName = './output/dwave/04_VAR_C_SR.txt';

N=4; NC=4;

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
i=0;
while (~feof(fid))
    SR = str2double(fgetl(fid));
    C = str2double(fgetl(fid));
    occ = str2double(fgetl(fid));
    qubits = str2double(fgetl(fid));
    chain = str2double(fgetl(fid));
    sol = (fgetl(fid));
    
    if (SR~=0.15) % only for a particular noise case
        continue;
    end
    
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
    
    fileNameSynth = ['./dataset/synth_data/Q_N4_NC4_C' num2str(C) '_SR0_L2.5.mat'];
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
    % {'exhaustive', 'matchEIG', 'matchALS', 'matchLift', 'matchBirkhoff', 'd-wave (ours)'};
    accuracies = [accuracies, [accCons, accMatchEIG, accMatchALS, accMatchLift, accBirdal, accDWave]'];
    PR = [PR [precision; recall]];
    PRCons = [PRCons [precisionCons; recallCons]];
         
    xValues = [xValues C];
end

fclose(fid);

set(gcf, 'currentaxes', h4); %subplot(2,2,1),
plot_multiple(xValues, accuracies, clrs, markers, dispNames, 1);
title('(d) Accuracy vs. completeness (with noise)', 'FontSize', 12);
