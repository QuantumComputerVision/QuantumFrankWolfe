
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

qEye = eye(N);
q = [qEye(:)' ...
    1 0 0 0 1 0 0 0 1 1 0 0 0 1 0 0 0 1 ]';
qSol = [qEye(:)' qEye(:)' qEye(:)' ]';

disp(['ground-truth: ' num2str(qSol')])
disp(['solution:     ' num2str(q')])


disp (['ground-truth: ', num2str(-qSol'*Qcons*qSol-s'*qSol)]);
disp (['solution:     ', num2str(-q'*Qcons*q-s'*q)]);

disp (['ground-truth (gnd): ', num2str(-qSol'*QGnd*qSol)]);
disp (['solution (gnd):     ', num2str(-q'*QGnd*q)]);

% 
% -q'*QGnd*q
% -q'*Q*q
% -q'*Qcons*q
% -qSol'*Qcons*qSol
% 
