
addpath('./permutations/');
addpath('./dataset/willow');
addpath('./dataset/willow/utils');
addpath('./dataset/willow/utils/rrwm');
addpath algorithms/

%% Load data
%clear
%startup
imgsets = {'Car','Duck','Motorbike','Winebottle'};
plotTitles = {'(a) Car','(b) Duck','(c) Motorbike','(d) Winebottle'};
% if show match
showmatch = false;

baseFolder = './dataset/real_data/willow/';
%folderNameDWaveOut = ['./output/dwave/willow/'];
folderNameDWaveOut = ['./output/dwave/willow/2000Q/'];
saveData = 0;

maxit = 250;
beta0 = 1;
verbose = 0;
isVisualize = 0;
isSave = 1;

N = 4; NC = 4;
lambda = 2.5;
completeness = 1;
swapRatio = 0 ;
nCams = NC;

meanAccuracies = [];
stdAccuracies = [];
%%

solutions = cell(4, length(1:35));

gndSols = cell(4,1);
for i=1:length(gndSols)
    gndSols{i} = eye(4);
end
qGndTruth = perms_cell_to_q(gndSols);

solverName = 'dwaveSAFWAL';
for imgInd = 1:4 %1:length(imgsets)
    imgset = imgsets{imgInd};
    folderToSave = ['./dataset/real_data/willow/' imgset '/'];
    
    datapath = sprintf('./dataset/willow/WILLOW-ObjectClass/%s/',imgset);
    savepath = sprintf('./dataset/willow/WILLOW-ObjectClass/%s',imgset);
    savefile = sprintf('%s/match_kpts.mat',savepath);
    imgList = dir([datapath,'*.png']);
    imgList = {imgList.name};
    viewList = dir([datapath,'*hypercols_kpts.mat']);
    viewList = {viewList.name};
    
    folderNameDWave = folderNameDWaveOut;
    folderNameDWave = [folderNameDWave imgset];
    
    xValues = [];
    accuracies = [];
    
    %figure('Renderer', 'painters', 'Position', [272 102 1048 834]);
%     t = tiledlayout(5,1);
%     t.Padding = 'none';
%     t.TileSpacing = 'none';
    k=0;
    for startIndex=1:35
        startIndex
        k=k+1;
        
        pMatch = runGraphMatchBatch_N_NC(datapath,viewList,N,NC,startIndex,'all',[],'wEdge',0.0,'thscore',0);
        if (isempty(pMatch))
            continue;
        end
        % prepare the unconstrained Q, not the constrained
        jMatch = runJointMatch(pMatch, [],'Method','prepare_uncons','univsize',N,...
            'rank',3,'lambda',lambda, 'initonly', 0, 'excludegeom', 0);
        
        Q = full(jMatch.Q);
        A = jMatch.A ;
        b = jMatch.b ;
        
        if(strcmp(solverName,'HCGM')==1)
            [PisQFWAL, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'HCGM', maxit, beta0, verbose, 0, isSave);
        elseif(strcmp(solverName,'dwaveHCGM')==1)
            [PisQFWAL, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'dwaveHCGM', maxit, beta0, verbose, 1, isSave);
        elseif(strcmp(solverName,'FWAL')==1)
            [PisQFWAL, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'FWAL', maxit, beta0, verbose, 0, isSave);
        elseif(strcmp(solverName,'dwaveFWAL')==1)
            [PisQFWAL, XFW, costMunkres, constraintResidual] = CoCoFW(-Q, N, A, b, 'dwaveFWAL', maxit, beta0, verbose, 1, isSave, isVisualize);
        elseif(strcmp(solverName,'dwaveSAFWAL')==1)
            [PisQFWAL, XFW, costMunkres, constraintResidual] = CoCoFW(Q, N, A, b, 'dwaveFWAL', maxit, beta0, verbose, 2, isSave, isVisualize);
        end

        PisQFWAL = perms_transform_first_to_eye(PisQFWAL);
        solutions{imgInd, k}=PisQFWAL;
        [precision, recall, ~, accQFWAL] = get_synth_solution_PR(PisQFWAL);
        accQFWAL
        %accuracies = [accuracies, [accCons, accMatchEIG accMatchLift accDWave]'];
        accuracies = [accuracies, accQFWAL'];
        
        xValues = [xValues startIndex];
    end
    
    meanAccuracy = mean(accuracies, 2);
    stdAccuracy = std(accuracies, [], 2);
    meanAccuracies = [meanAccuracies meanAccuracy];
    stdAccuracies = [stdAccuracies stdAccuracy];
end
meanAccuracies
stdAccuracies
