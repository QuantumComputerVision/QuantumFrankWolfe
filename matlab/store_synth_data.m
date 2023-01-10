function [data] = store_synth_data(fileName, Qorg, Q, QGnd, Qcons, s, PijsGnd, Pijs, Pis, N, G, I, Iupper, nCams, lambda, swapRatio, completeness)

%%%% IF YOU WANT TO SAVE %%%%
data.Qorg = Qorg;
data.Q = Q;
data.PijsGnd = PijsGnd;
data.Pijs = Pijs;
data.Pis = Pis;
data.N = N; % # of points
data.G = G; % # of points
data.I = I; % # of points
data.Iupper = Iupper; % # of points
data.nCams = nCams;
if (exist('QGnd','var') && ~isempty(QGnd))
    data.QGnd = QGnd;
end
if (exist('Qcons','var') && ~isempty(Qcons))
    data.Qcons = sparse(Qcons);
end
if (exist('s','var') && ~isempty(s))
    data.s = s;
end
if (exist('PijsGnd','var') && ~isempty(PijsGnd))
    data.PijsGnd = PijsGnd;
end
if (exist('Pijs','var') && ~isempty(Pijs))
    data.Pijs = Pijs;
end
if (exist('Pis','var') && ~isempty(Pis))
    data.Pis = Pis;
end

if (exist('lambda','var') && ~isempty(lambda))
    data.lambda = lambda;
end
data.lambda = lambda;
data.swapRatio = swapRatio;
data.completeness = completeness;

[data.QconsRows, data.QconsCols, data.QconsValues] = find(Qcons);
[data.QRows, data.QCols, data.QValues] = find(Q);
[data.QGndRows, data.QGndCols, data.QGndValues] = find(QGnd);

if (~isempty(fileName))
    save(fileName, 'data');
end

end