import os
import sys
import math 
import numpy as np 
import json 
from ttictoc import tic,toc
from scipy.io import loadmat 
from scipy.io.matlab import savemat 
import dimod 
import dwavebinarycsp 
import dwave.inspector 
import minorminer 
from dwave.system import DWaveSampler, EmbeddingComposite 

# sample qubo on DWave 
# Q_weights: a real symmetric matrix of QUBO weights; its main diagonal contains qubit biases and non-diagonal elements are qubit couplings 
# CHAIN_STR: chain strength (if CHAIN_STR < 0, the maximum chain length criterion is used) 
# N_anneals: number of anneals 
def dwave_sample_qubo(solver_, sampler, Q_weights, CHAIN_STR, N_anneals): 
    
    # parameters that change rarely 
    epsilon         = 0.001                                                # threshold for an effective zero value 
    offset          = 0.0                                                  # bqm offset (should remaing zero) 
    vartype = dimod.BINARY                                                 # possible measurement outcomes are 0 and 1 
    anneal_params   = dict(anneal_schedule=[[0.0, 0.0], [20.0, 1.0]])      # default annealing schedule; no pause during the anneal 
    # 
    
    d = len(Q_weights) 
    Q = [] 
    C = [] 
    Q_RECT = Q_weights 
    C_RECT = np.diag(Q_RECT) 
    
    #print(d) 
    #print(Q_RECT) 
    #print(C_RECT) 
    
    # the QUBO weights (linear and quadratic) 
    linear = {} 
    quadratic = {} 
    
    # qubit bises 
    for i in range(0, d): 
        linear[i+1] = C_RECT[i] 
    
    # qubit couplings (only non-zero values); "quadratic" must be upper triangular matrix 
    for i in range(0, d): 
        for j in range(0, d): 
            if  ( (i < j) and (abs(Q_RECT[i, j]) > epsilon) ): 
                quadratic[(i+1, j+1)] =  2 * Q_RECT[i, j] 
    
    #print(linear) 
    #print(quadratic) 
    
    # initialise the problem and pass it to Adv4.1 
    bqm = dimod.BinaryQuadraticModel(linear, quadratic, offset, vartype) 
   
    # set chain strength according to the maximun chain length criterion 
    if CHAIN_STR < 0: 
        print('minorminer')
        tic()
        embedding_ = minorminer.find_embedding(list(bqm.quadratic.keys()), solver_.edgelist, return_overlap=1) 
        embedding_current = embedding_[0] 
        ch_lengths_ = {} 
        lenv = 0 
        for key in embedding_current: 
            ch_lengths_[lenv] =  len(embedding_current[key]) 
            lenv += 1
        CHAIN_STR = max(ch_lengths_.values()) + 0.5
        print(toc())
        
    # sample QUBO N_anneals times 
    print('sampling')
    tic()
    sampleset = sampler.sample(bqm, chain_strength = CHAIN_STR, num_reads = N_anneals, **anneal_params) 
    print(toc())
    
    print('lowest_energy_sample')
    tic()
    # select the lowest-energy sample 
    best_sample = sampleset.first.sample 
    print(toc())
    # return the lowest-energy solution 
    lowest_energy_sample = [0 for i in range(0, d)] 
    # lowest_energy_sample = np.asarray(np.zeros(d))
    
    for k in range(0, d): 
        lowest_energy_sample[k]  = best_sample[k+1] 
    
    
    return lowest_energy_sample; 
# end of dwave_sample_qubo 


class DWaveSolver:
    # default constructor
    def __init__(self):
        self.input_folder = "../datasets/synth/"
        self.output_folder = "../output/dwaveQFW/"
        self.CHAIN_STR_ = -1.0  # negative value: use the maximum chain length criterion; otherwise, [positive] CHAIN_STR_ is used
        self.Nsamples_ = 50  # number of samples/annealings

        # initialize sampler and solver
        self.solver_ = DWaveSampler(solver='Advantage_system4.1')
        self.sampler = EmbeddingComposite(DWaveSampler(self.solver_))
        
    def sample(self, input_file_name):
        input_file = self.input_folder + input_file_name
        ext = '.' + os.path.realpath(input_file_name).split('.')[-1:][0]
        output_file = self.output_folder + input_file_name.replace('.mat', '.csv')

        matfile = loadmat(input_file)
        x = matfile['data']
        Q = x['Q'][0][0]
        Q_weights_ = Q

        # sample qubo on DWave
        lowest_energy_sample = []
        lowest_energy_sample = dwave_sample_qubo(self.solver_, self.sampler, Q_weights_, self.CHAIN_STR_, self.Nsamples_)

        # save the result in a ".csv" file
        np.savetxt(output_file, lowest_energy_sample, fmt='%d', newline=" ")

def main():
    input_file_name     = sys.argv[1]
    solver = DWaveSolver()
    solver.sample(input_file_name)
    
#def main():
#
#    print('init')
#    
#    tic('tic0')
#    # parameters 
#    input_folder        = "../datasets/synth/" 
#    output_folder       = "../output/dwaveQFW/"
#    input_file_name     = sys.argv[1]
#    input_file          = input_folder +  input_file_name
#    CHAIN_STR_          = -1.0               # negative value: use the maximum chain length criterion; otherwise, [positive] CHAIN_STR_ is used 
#    Nsamples_           = 50                 # number of samples/annealings 
#
#    ext = '.'+ os.path.realpath(input_file_name).split('.')[-1:][0]
#    output_file = output_folder + input_file_name.replace('.mat','.csv')
#    
#    print('input:' + input_file)
#    print('output:' + output_file)
#
#    print('initialization')
#    tic()
#    # initialize sampler and solver
#    solver_ = DWaveSampler(solver='Advantage_system4.1')
#    sampler = EmbeddingComposite(DWaveSampler(solver_))
#    print(toc())
#    
#    # load the data 
#    matfile = loadmat(input_file) 
#    x = matfile['data'] 
#    Q = x['Q'][0][0] 
#    Q_weights_ = Q 
#
#    print('loading successful')
#    
#    # sample qubo on DWave 
#    lowest_energy_sample = [] 
#    lowest_energy_sample = dwave_sample_qubo(solver_, sampler, Q_weights_, CHAIN_STR_, Nsamples_) 
#    
#    # save the result in a ".csv" file
#    np.savetxt(output_file, lowest_energy_sample, fmt='%d', newline=" ") 
#    
#    print(output_file)
#    
#    print('ALL OF IT')
#    print(toc('tic0'))

if __name__ == "__main__":
    main()