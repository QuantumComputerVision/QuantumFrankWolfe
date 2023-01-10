
isSA = 0; % use simulated annealing
Q = load('c:/Users/tolga/Downloads/Q.mat'); % some Q matrix
Q = Q.Q;
Q = Q(2:end,2:end); % ignore the first row/col

% dump the Q matrix to 'dummy'
data.Q = Q;
save('../datasets/synth/dummy.mat', 'data');

% initialize the dwave module as we always do
% to change the python path call as (QPath is not used at all) :
% [dwaveModule] = load_dwave_module(QPath, pythonPath);
dwave = load_dwave_module();

% this is the new function that is to be called in the beginning of FWAL to
% initialize the sparsity structure (this is expected to take long, but
% will be done only once thanks to the sparsity of Q remaining unchanged.)
dwave.prepare_sampler('dummy.mat'); 

% now sample. note that the last argument is 1 indicating that it should
% use the 'fixed' (a.k.a. precomputed) embedding
% note that the two lines below are exactly what we do if we call 
% solve_dwave_module
tic(); dwave.sample('dummy.mat',0,1); toc(); 

% sampled points are now written to the output csv. so load it and use x as
% desired :)
x = double(dlmread('../output/dwaveQFW/dummy.csv'))';

