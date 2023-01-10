import os
import sys
import math 
import numpy as np 
import json 
from ttictoc import tic,toc
from scipy.io import loadmat 
from scipy.io.matlab import savemat 
import dimod 
import neal
import dwavebinarycsp 
import dwave.inspector 
import minorminer 
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite 

input_folder = "../datasets/synth/"
output_folder = "../output/dwaveQFW/"
CHAIN_STR_ = 20  # negative value: use the maximum chain length criterion; otherwise, [positive] CHAIN_STR_ is used
Nsamples_ = 300  # number of samples/annealings
anneal_params = dict(anneal_schedule=[[0.0, 0.0], [10.0, 0.5], [110.0, 0.5], [120.0, 1.0]]) 
solver_ = DWaveSampler(solver='Advantage_system4.1')
sampler = EmbeddingComposite(DWaveSampler(solver_))
samplerNeal = neal.SimulatedAnnealingSampler()
fixedSampler = None
embeddingLoaded = 0

def load_Q(input_file_name):
    input_file = input_folder + input_file_name
    #print(input_file)
    matfile = loadmat(input_file)
    x = matfile['data']
    Q = x['Q'][0][0]
    Q_weights_ = Q
    
    return Q

def get_bqm(Q_weights):
    # parameters that change rarely 
    epsilon         = 0.001                                                # threshold for an effective zero value 
    offset          = 0.0                                                  # bqm offset (should remaing zero) 
    vartype = dimod.BINARY                                                 # possible measurement outcomes are 0 and 1 
    #anneal_params   = dict(anneal_schedule=[[0.0, 0.0], [20.0, 1.0]])      # default annealing schedule; no pause during the anneal 
    
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
    
    return bqm
    
def find_chain_strength_from_embedding(embedding_):
    embedding_current = embedding_[0] 
    ch_lengths_ = {} 
    lenv = 0 
    for key in embedding_current: 
        ch_lengths_[lenv] =  len(embedding_current[key]) 
        lenv += 1
    CHAIN_STR = max(ch_lengths_.values()) + 3
    return CHAIN_STR

def find_chain_strength(bqm,CHAIN_STR,isSA=False):
    # set chain strength according to the maximun chain length criterion 
    if CHAIN_STR < 0 and (not isSA):
        # print('minorminer')
        # tic()
        embedding_ = minorminer.find_embedding(list(bqm.quadratic.keys()), solver_.edgelist, return_overlap=1) 
        CHAIN_STR = find_chain_strength_from_embedding(embedding_)

    return CHAIN_STR

def load_sampler(embedding_file_name='emdebbing.emb'):
    embedding_retrieved = {} 
    with open(embedding_file_name) as json_file: 
        embedding_retrieved = json.load(json_file) 

    embedding_retrieved = {int(k):[int(i) for i in v] for k,v in embedding_retrieved.items()} 
    
    fixedSampler = FixedEmbeddingComposite(DWaveSampler(solver='Advantage_system4.1'), embedding_retrieved) 
    
    embeddingLoaded = 1
    
    return fixedSampler, embedding_retrieved
    
def prepare_sampler(input_file_name):
    
    Q = load_Q(input_file_name)
    bqm = get_bqm(Q)
    embedding_ = minorminer.find_embedding(list(bqm.quadratic.keys()), solver_.edgelist, return_overlap=1) 
    embedding_file_name = 'emdebbing.emb'
    with open(embedding_file_name, 'w') as fp: 
        json.dump(embedding_[0], fp) 
    
    fixedSampler = load_sampler()
    
    CHAIN_STR = find_chain_strength_from_embedding(embedding_)
    
    return fixedSampler
    
# sample qubo on DWave 
# Q_weights: a real symmetric matrix of QUBO weights; its main diagonal contains qubit biases and non-diagonal elements are qubit couplings 
# CHAIN_STR: chain strength (if CHAIN_STR < 0, the maximum chain length criterion is used) 
# N_anneals: number of anneals 
def dwave_sample_qubo(solver_, sampler, Q_weights, CHAIN_STR, N_anneals, isSA, useEmbedding=False): 
   
    d = len(Q_weights) 
    bqm = get_bqm(Q_weights)
    
    # print(toc())
        
    # sample QUBO N_anneals times 
    # print('sampling')
    # tic()
    if (isSA):
        sampleset = samplerNeal.sample(bqm)
    elif (useEmbedding):
        fixedSampler, embedding_ = load_sampler()        
        # CHAIN_STR = find_chain_strength_from_embedding(embedding_)
        # print('embedding loaded! running fixedSampler')
        sampleset = fixedSampler.sample(bqm, chain_strength = CHAIN_STR, num_reads = N_anneals, **anneal_params) 
    else:
        # sets the CHAIN_STR variable
        CHAIN_STR = find_chain_strength(bqm,CHAIN_STR,isSA)
        sampleset = sampler.sample(bqm, chain_strength = CHAIN_STR, num_reads = N_anneals, **anneal_params) 
   #  print(toc())
    #sampleset.wait()
    #print(sampleset.info["timing"])
    # print('lowest_energy_sample')
    # tic()
    # select the lowest-energy sample 
    best_sample = sampleset.first.sample 
    # print(toc())
    # return the lowest-energy solution 
    lowest_energy_sample = [0 for i in range(0, d)] 
    # lowest_energy_sample = np.asarray(np.zeros(d))
    
    for k in range(0, d): 
        lowest_energy_sample[k]  = best_sample[k+1] 
    
    
    return lowest_energy_sample; 
# end of dwave_sample_qubo 


def sample(input_file_name,isSA,useEmbedding):
    
    Q_weights_ = load_Q(input_file_name)

    # sample qubo on DWave
    lowest_energy_sample = []
    lowest_energy_sample = dwave_sample_qubo(solver_, sampler, Q_weights_, CHAIN_STR_, Nsamples_, isSA, useEmbedding)

    # ext = '.'+ os.path.realpath(input_file_name).split('.')[-1:][0]
    output_file = output_folder + input_file_name.replace('.mat','.csv')
    # save the result in a ".csv" file
    np.savetxt(output_file, lowest_energy_sample, fmt='%d', newline=" ")
    
def main():
    print('main')
    
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