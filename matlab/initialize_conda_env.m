function [] = run_dwave_python(QPath, pythonPath)
%% doesn't work, under development
% pyenv('Version', 'C:\Users\tolga\anaconda3\envs\jupenv\python.exe');

if (~exist('pythonPath', 'var'))
    pythonPath = 'C:\Users\tolga\anaconda3\envs\jupenv\python.exe';
end

if (~exist('QPath', 'var'))
    QPath = 'dummy.mat';
end
pathofDocument = fileparts(which('../python/run_dwave.py'));
if count(py.sys.path,pathofDocument) == 0
    insert(py.sys.path,int32(0),pathofDocument);
end

py_root_useFromMATLAB = fileparts(pythonPath);
ENV = getenv('PATH');
ENV = strsplit(ENV, ';');
items_to_add_to_path = {
    fullfile(py_root_useFromMATLAB, 'Library', 'mingw-w64', 'bin')
    fullfile(py_root_useFromMATLAB, 'Library', 'usr', 'bin')
    fullfile(py_root_useFromMATLAB, 'Library', 'bin')
    fullfile(py_root_useFromMATLAB, 'Scripts')
    };
ENV = [items_to_add_to_path(:); ENV(:)];
ENV = unique(ENV, 'stable');
ENV = strjoin(ENV, ';');
setenv('PATH', ENV);
% clear classes
%module_to_load = '../python/run_dwave.py FW_Q_t=2.mat';
module_to_load = 'runDWaveFunc';
%dwaveModule = py.importlib.import_module('runDWaveFunc')
% pyrunfile(module_to_load);
python_module_to_use = py.importlib.import_module(module_to_load);
%py.importlib.reload(python_module_to_use);
end