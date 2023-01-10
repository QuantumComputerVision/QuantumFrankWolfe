function [x] = solve_dwave(Q,isSA,pythonPath)

if (~exist('pythonPath', 'var'))
    pythonPath = 'C:\Users\tolga\anaconda3\envs\jupenv\python.exe';
end

if (~exist('isSA', 'var')) isSA = 0; end

data.Q = Q;
save('../datasets/synth/dummy.mat', "data");
% run_dwave_python('dummy.mat',pythonPath,isSA);
dwaveModule = load_dwave_module();
x = solve_dwave_module(dwaveModule, Q,isSA);
% x = double(dlmread('../output/dwaveQFW/dummy.csv'))';
%x = double(dlmread('dummy.csv'));

end