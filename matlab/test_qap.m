%% Choose data
% NOTE: You need to download data from QAPLIB and locate them to under the
% "FilesQAP/data/qapdata/" folder (resp. TSPLIB, "FilesQAP/data/tspdata/").
% Links for QAPLIB: 
%   http://anjos.mgi.polymtl.ca/qaplib/inst.html
%   https://www.opt.math.tugraz.at/qaplib/inst.html
% Links for TSPLIB: 
%   http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/index.html
%   http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95

% dataname = 'qapdata/chr12a'; % you can choose other data files as well
% dataname = 'tspdata/ulysses16';

addpath permutations/ 
addpath algorithms/
addpath FilesQAP/

global pythonPath
global QPath

QPath = 'dummy.mat';
pythonPath = 'C:\Users\tolga\anaconda3\envs\jupenv\python.exe';

% datasets = {'chr12a','chr12b','chr12b','chr15a','chr15b','chr15c','esc16a','esc16b','esc16c','esc16d','esc16e','esc16g','esc16h','esc16i','esc16j','had12','had14','had16','nug12','nug14','nug15','nug16a','nug16b','nug17','rou12','rou15','scr12','scr15','tai12a','tai12b','tai15a','tai15b','tai17a'};
datasets = {'chr12a'};

maxit = 100;
beta0 = 1;
verbose = 0;
isSave = 0;
isVisualize = 0;

solversToRun = {'FWAL'}; % {'dwaveSA', 'dwvave', 'FWAL', 'dwaveFWAL', 'HCGM', 'dwaveHCGM'};
%solversToRun = {'dwaveFWAL'};

% numProblems = size(N3.W1, 1);
numProblems = length(datasets);

for i=1:numProblems
    
    dataname = datasets{i};
    fprintf(['*****',dataname,'*****\n']);

    [C,D,OPT,P_OPT] = qapread(['./dataset/QAP/',dataname]);
    energiesOPT(i) = OPT;

    N = size(C,1);
    Q = kron(D,C);
    Q = full(0.5*(Q + Q'));
    % Q = Q/norm(Q);

    [A, b] = compute_constraint_system(N);

    for sj=1:length(solversToRun)
        solverName = solversToRun{sj};
        if (strcmp(solverName,'exh')==1)
            qSol = solve_exhaustive_QGM(Q, N);
            PSol = qSol{1};
            qSol = PSol(:);
            costExh = qSol'*Q*qSol;
            energiesGnd(i) = costExh;
            disp(['cost Gnd: ' num2str(costExh)]);
        elseif(strcmp(solverName,'dwave')==1)
            qSol = solve_dwave(Q);
            costDWave = qSol'*Q*qSol;
            energiesDWave(i) = costDWave;
            disp(['cost DWave: ' num2str(costDWave)]);
        elseif(strcmp(solverName,'dwaveSA')==1)
            qSol = solve_dwave(Q, 1);
            Pis = perms_q_to_cell(qSol, N);
            [Pis, costs] = round_munkres(Pis);
            qSol = perms_cell_to_q(Pis);
            costDWaveSA = qSol'*Q*qSol;
            energiesDWaveSA(i) = costDWaveSA;
            disp(['cost DWave: ' num2str(costDWaveSA)]);
        else
            if(strcmp(solverName,'HCGM')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'HCGM', maxit, beta0, verbose, 0, isSave, isVisualize);
            elseif(strcmp(solverName,'dwaveHCGM')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'dwaveHCGM', maxit, beta0, verbose, 1, isSave, isVisualize);                
            elseif(strcmp(solverName,'FWAL')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'FWAL', maxit, beta0, verbose, 0, isSave, isVisualize);
            elseif(strcmp(solverName,'dwaveFWAL')==1)
                [Pis, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'dwaveFWAL', maxit, beta0, verbose, 1, isSave, isVisualize);
            end
            Pis = Pis{1};
            q = Pis(:);
            cost = transpose(q)*Q*q;
            % eval(['cost' solversToRun ' = transpose(q)*Q*q']);
            eval(['energies' solverName '(i) = cost;']);
            % energiesCoCo(i) = cost;
            costCocoObj = Q(:)'*XFW(:); %transpose(());
            eval(['energies' solverName 'CoCoObj(i) = costCocoObj;']);

            disp(['cost CoCo: ' num2str(cost) ', CoCo-obj: ' num2str(costCocoObj) ', OPT:' num2str(OPT)]);
        end
        
    end
    
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

