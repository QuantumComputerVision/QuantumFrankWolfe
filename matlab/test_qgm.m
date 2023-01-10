addpath permutations/ 
addpath algorithms/

global pythonPath
global QPath

QPath = 'dummy.mat';
pythonPath = 'C:\Users\tolga\anaconda3\envs\jupenv\python.exe';

N3= load('./dataset/QGM/N3.mat'); % 9x9
N4= load('./dataset/QGM/N4.mat'); % 16x16

N = 3;
maxit = 200;
beta0 = 1;
verbose = 0;
isSave = 0;
isVisualize = 0;

%solversToRun = {'exh', 'dwaveSA', 'dwvave', 'FWAL', 'dwaveFWAL', 'HCGM', 'dwaveHCGM'};
solversToRun = {'dwaveFWAL'};

[A, b] = compute_constraint_system(N);

numProblems = size(N3.W1, 1);
% energiesGnd = zeros(numProblems, 1);
% energiesDWave = zeros(numProblems, 1);
% energiesDWaveSA = zeros(numProblems, 1);
% energiesCoCo = zeros(numProblems, 1);
% energiesCoCoObj = zeros(numProblems, 1);

if (N==3)
    Qs = N3.W1;
    cs = N3.c1;
%     Qsol3 = load('../output/QGM\test_qgm-results-N=3.mat');
%     energiesFWAL = Qsol3.energiesFWAL;
%     energiesFWALCoCoObj = Qsol3.energiesFWALCoCoObj;
%     energiesGnd = Qsol3.energiesGnd;
%     dwaveResults = load('QGGM_N3_I400_b1.mat');
%     energiesdwaveFWAL = dwaveResults.energiesdwaveFWAL;
%     energiesdwaveFWALCoCoObj = dwaveResults.energiesdwaveFWALCoCoObj;
elseif (N==4)
    Qs = N4.W1;
    cs = N4.c1;
%     Qsol4 = load('../output/QGM\test_qgm-results-N=4.mat');
%     energiesGnd = Qsol4.energiesGnd;
%     energiesFWAL = Qsol4.energiesFWAL;
%     energiesFWALCoCoObj = Qsol4.energiesFWALCoCoObj;
%     dwaveResults = load('QGGM_N4_I400_b1.mat');
%     energiesdwaveFWAL = dwaveResults.energiesdwaveFWAL;
%     energiesdwaveFWALCoCoObj = dwaveResults.energiesdwaveFWALCoCoObj;
end

for i=7:numProblems
    Q3 = squeeze(Qs(i, :, :)) + diag(squeeze(cs(i, :)));
    Q3 = 0.5*(Q3+Q3');

    for sj=1:length(solversToRun)
        solverName = solversToRun{sj};
        if (strcmp(solverName,'exh')==1)
            qSol = solve_exhaustive_QGM(Q3, N);
            PSol = qSol{1};
            qSol = PSol(:);
            costExh = qSol'*Q3*qSol;
            energiesGnd(i) = costExh;
            disp(['cost Gnd: ' num2str(costExh)]);
        elseif(strcmp(solverName,'dwave')==1)
            qSol = solve_dwave(Q3);
            costDWave = qSol'*Q3*qSol;
            energiesDWave(i) = costDWave;
            disp(['cost DWave: ' num2str(costDWave)]);
        elseif(strcmp(solverName,'dwaveSA')==1)
            qSol = solve_dwave(Q3, 1);
            Pis = perms_q_to_cell(qSol, N);
            [Pis, costs] = round_munkres(Pis);
            qSol = perms_cell_to_q(Pis);
            costDWaveSA = qSol'*Q3*qSol;
            energiesDWaveSA(i) = costDWaveSA;
            disp(['cost DWave: ' num2str(costDWaveSA)]);
        else
            if(strcmp(solverName,'HCGM')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q3, N, A, b, 'HCGM', maxit, beta0, verbose, 0, isSave, isVisualize);
            elseif(strcmp(solverName,'dwaveHCGM')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q3, N, A, b, 'dwaveHCGM', maxit, beta0, verbose, 1, isSave, isVisualize);                
            elseif(strcmp(solverName,'FWAL')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q3, N, A, b, 'FWAL', maxit, beta0, verbose, 0, isSave, isVisualize);
            elseif(strcmp(solverName,'dwaveFWAL')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q3, N, A, b, 'dwaveFWAL', maxit, beta0, verbose, 1, isSave, isVisualize);
            end
            Pis = Pis{1};
            q = Pis(:);
            cost = transpose(q)*Q3*q;
            % eval(['cost' solversToRun ' = transpose(q)*Q3*q']);
            eval(['energies' solverName '(i) = cost;']);
            % energiesCoCo(i) = costCoco;
            costCocoObj = Q3(:)'*XFW(:); %transpose(());
            eval(['energies' solverName 'CoCoObj(i) = costCocoObj;']);

            disp(['cost CoCo: ' num2str(cost) ', CoCo-obj: ' num2str(costCocoObj)]);
        end
        
    end

    %q3 = qubo_solvers(Q3);
    
end

%save('QGGM_N4_I400_b1.mat');
return;

% FWALBinResults = load('../output/QGM/test_qgm-results-N=3.mat');

base = -min(energiesGnd);
xpt = 1:10;
combined = [energiesGnd(:) energiesDWaveSA(:) energiesFWALCoCoObj(:) energiesFWAL(:) energiesdwaveFWALCoCoObj(:) energiesdwaveFWAL(:)];

combined = combined - repmat(energiesGnd(:),1,size(combined,2));
mean(combined,1);

figure, beautify_plot;
hb = bar(xpt, combined+0.01, 'grouped');



figure, beautify_plot;
bar(xpt, energiesGnd+base,0.125,'LineWidth',1);
hold on, bar(xpt, energiesDWaveSA(:)+base,0.25,'LineWidth',1);
hold on, bar(xpt, energiesFWALCoCoObj(:)+base,0.4,'LineWidth',1);
% hold on, plot(energiesFWAL,'LineWidth',3,'Marker','o');
hold on, bar(xpt, energiesFWAL(:)+base,0.55,'LineWidth',1);
hold on, bar(xpt, energiesdwaveFWALCoCoObj(:)+base,0.7,'LineWidth',1);
hold on, bar(xpt, energiesdwaveFWAL(:)+base,0.85,'LineWidth',1);
legend({'Gnd-Truth', 'SA', 'FWAL-Relaxed', 'FWAL', 'Q-FWAL-Relaxed', 'Q-FWAL'});
grid on;
xlabel('QGM experiment ID');
ylabel('Normalized Energy ($\mathbf{q}^\top\mathbf{Q}\mathbf{q}$)','interpreter','latex');

return;

figure, beautify_plot;
plot(energiesGnd+base,'LineWidth',3,'Marker','x');
hold on, plot(energiesFWALCoCoObj+base,'LineWidth',3,'Marker','o', 'LineStyle', '-.');
% hold on, plot(energiesFWAL,'LineWidth',3,'Marker','o');
hold on, plot(energiesFWAL+base,'LineWidth',3,'Marker','o', 'LineStyle', '-.');
hold on, plot(energiesdwaveFWALCoCoObj+base,'LineWidth',3,'Marker','+');
hold on, plot(energiesdwaveFWAL+base,'LineWidth',3,'Marker','+');
legend({'Gnd-Truth', 'FWAL-Relaxed', 'FWAL', 'Q-FWAL-Relaxed', 'Q-FWAL'});
grid on;
xlabel('QGM experiment ID');
ylabel('Normalized Energy ($\mathbf{q}^\top\mathbf{Q}\mathbf{q}$)','interpreter','latex');

axes('position',[.15 .175 .25 .25]);
box on; beautify_plot;
plot(energiesGnd,'LineWidth',3,'Marker','x');
hold on, plot(energiesFWALCoCoObj,'LineWidth',3,'Marker','o', 'LineStyle', '-.');
hold on, plot(energiesFWAL,'LineWidth',3,'Marker','o', 'LineStyle', '-.');
hold on, plot(energiesdwaveFWALCoCoObj,'LineWidth',3,'Marker','+');
hold on, plot(energiesdwaveFWAL,'LineWidth',3,'Marker','+');
axis tight;
grid on;

% try to compute the metric in table 1
energiesCoCoShift = energiesCoCo - energiesGnd;
energiesCoCoObjShift = energiesCoCoObj - energiesGnd;


% for i=1:numProblems
%     energiesCoCo(i) = energiesCoCo(i) - energiesGnd(i) ;
%     energiesCoCoObj(i) = energiesCoCoObj(i) - energiesGnd(i) ;
% end

% 
% % Q4 = N4.W1 + diag(N4.c1);
% numProblems = size(N4.W1, 1);
% for i=1:numProblems
%     Q4 = squeeze(N4.W1(i, :, :)) + diag(squeeze(N4.c1(i, :)));
%     q2 = solve_exhaustive(Q4, 1, 4);
%     q4 = q2{1}
% end

