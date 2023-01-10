function [x] = solve_dwave_module(dwaveModule, Q,isSA)

if (~exist('isSA', 'var')) isSA = 0; end

data.Q = Q;
save('../datasets/synth/dummy.mat', 'data');
%run_dwave_python('dummy.mat');
dwaveModule.sample('dummy.mat',isSA);
x = double(dlmread('../output/dwaveQFW/dummy.csv'))';
%x = double(dlmread('dummy.csv'));

end