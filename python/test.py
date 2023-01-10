import numpy as np
from scipy.io import loadmat
from scipy.io import savemat

x = loadmat('../dataset/synth_data_output.mat')

x = x['data']

# The synchronization problem:
Pis = np.squeeze(x['Pis'])[()][0] # ground truth solution to the problem (absolute permutation matrices)
Pijs = np.squeeze(x['Pijs'])[()][0] # noisy observations (relative permutations)
PijsGnd = x['PijsGnd'] # ground truth solution to the problem in terms of relative permutations
N = np.squeeze(x['N'])[()] # number of points at each node (number of points in a grid)
G = x['G'] # the problem graph as an adjacency matrix (G_ij = 1 if i maps to j
I = np.squeeze(x['I'])[()] # edges as [rows' ; cols']
nCams = np.squeeze(x['nCams'])[()] # number of cameras

QGndSparse = x['QGnd'][0][0] # ground truth Q matrix to the problem as sparse matrix
QGnd = QGndSparse.todense() # ground truth solution to the problem as full matrix
QGndSparse

# access the sparse storage
rows = QGndSparse.tocoo().row
columns = QGndSparse.tocoo().col
values = QGndSparse.tocoo().data

Q = x['Q'][0][0] # unconstrained Q matrix (sparse), s is not required.
Qcons = x['Qcons'][0][0] # constrained Q matrix (sparse), s is required
s = x['s'][0][0] # s for the constrained case

lambdaReg = x['lambda'][0][0] # regularizer used in construction of Qcons

# to access the sparse storage of Q
rows = Q.tocoo().row
columns = Q.tocoo().col
values = Q.tocoo().data

rows = x['QconsRows'][0][0]
cols = x['QconsCols'][0][0]
values = x['QconsValues'][0][0]

# to access the sparse storage of Qcons
rows = Qcons.tocoo().row
columns = Qcons.tocoo().col
values = Qcons.tocoo().data
sizeQcons = Qcons.shape

rows = x['QGndRows'][0][0]
cols = x['QGndCols'][0][0]
values = x['QGndValues'][0][0]

# savemat example:
qresult = '1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,1,1,0,0,0,1,0,0,0,1'
numberOfBitsUsed = 47
chainStrength = 2
lambdaUsed = 2

outputDic = {'qresult': qresult, 'numberOfBitsUsed': numberOfBitsUsed, 'chainStrength': chainStrength, 'lambdaUsed': lambdaUsed}
savemat("../output/ouptut_N4_NC4_C1_SR0_L25.mat", outputDic)

X = np.fromfile('C:/Users/tolga/Box Sync/gitforks/QuantumSync/output/dwave/SAMPLES_N4_NC4_200RUNS.csv', sep="", dtype='%s')

np.save('C:/Users/tolga/Box Sync/gitforks/QuantumSync/output/dwave/SAMPLES_N4_NC4_200RUNS.npy', X)
tmp = {varname: mydata[varname] for varname in mydata.dtype.names}