import numpy as np
import umap
import umap.umap_ as umap

import cdlib
from cdlib import algorithms
import networkx as nx
from  scipy import sparse

import warnings
warnings.filterwarnings('ignore')

import time
import seaborn as sns;
import matplotlib.pyplot as plt

import os
import ray
from scipy.sparse import csr_matrix
import csv
import os
from networkx import Graph
import re
import pickle

# divides a sequence type into chunks

def chunker(seq, size):
    """
    input: sequence (eg. list)
    ouptput: sequence divided in chunks
    parameters: size: determines the length of the chunk
    """
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

# helper decorator function for benchmarking

def my_timer(func):
    """
    input: function to be benchmarked
    ouptput: execution time of benchmarked function
    """
    def wrapper(*args,**kwargs):
        t_start = time.time()
        result = func(*args,**kwargs)
        t_end = time.time() - t_start
        print('{} took {}s'.format(func.__name__, t_end))

        return result

    return wrapper

def build_graph(XYZ,k):
    """
    Input: an array of (x,y,z) coordinates
    Output: the weighted adjacency matrix of the UMAP graph representation
    Parameters: k is the most important parameter in the umap fuzzy_simplicial_set function.
    It will determine how sparse the final graph will be.
    """
#    umap.umap_.fuzzy_simplicial_set
    adj = umap.fuzzy_simplicial_set(
        XYZ,
        n_neighbors=k, # this parameter has to be fine-tuned
        random_state=np.random.RandomState(seed=42),
        metric='l2',
        metric_kwds={},
        knn_indices=None,
        knn_dists=None,
        angular=False,
        set_op_mix_ratio=1.0,
        local_connectivity=2.0,
        verbose=False,

        )

    return adj


def build_communities(adj):
    """
    Input: the weighted graph adjacency matrix
    Output: a list of communities, each one a represented as a list object
    leiden algorithm as implemented in the cdlib library.
    """
    # generate a graph networkx obj
    g = nx.from_scipy_sparse_matrix(adj)
    # get list of edges from graph
    eset = [(u, v) for (u, v, d) in g.edges(data=True)]
    # get list of weights from edges
    weights = [d['weight'] for (u, v, d) in g.edges(data=True)]
    # find communities using Leiden alg
    leiden_coms = algorithms.leiden(g,weights=weights) # check if the algo is stochastic, in that case set rnd generator
    # a list of lists of nodes
    return leiden_coms.communities


#parrarel graph aggregator

@ray.remote
def aggregate_graphs(graph1,graph2):
    return graph1 + graph2

@ray.remote

def read_and_prepare_graph_and_communities(folder,num,k,number_of_beads_per_structure):

    '''
    input: csv file with coordinates for single structure
    ouput: graph, list of communities and id of structure (from csv name)
    parameters: k - is the most important parameter in the umap fuzzy_simplicial_set function.
    It will determine how sparse the final graph will be.

    '''



    #obtain id of processed_structur


    filename = 'cf_' + str(num).zfill(6) + '.coords.csv'
    file = os.path.join(folder,filename)


    # read csv into np.array
    coordinates = np.genfromtxt(file, delimiter= ',')
    # get columns for x,y,z coordinates
    coordinates_xyz = coordinates[:,:3]
    # build a graph from x,y,z
    graph = build_graph(coordinates_xyz,k)
    # detect communities
    communities =  build_communities(graph[0])
    # communities are list of lists of lists : community / beads
    # obtained communities are used for as input to build a complete graph for given structure
    for community_index in range(len(communities)):
        # for the first community build graph
        if community_index == 0:
            community_graph = nx.complete_graph(communities[community_index])
        else:
         # for the following communities update graph
            community_graph.update(nx.complete_graph(communities[community_index]))
         #once done - return graph , list of communities and id of processed structure
    return nx.to_scipy_sparse_matrix(community_graph,nodelist=range(number_of_beads_per_structure)) , communities , num


@my_timer

# parallelizing the whole process

def process_csvs(folder,cores,k,number_of_beads_per_structure,number_of_structures):

    '''
    input: folder with csvs
    ouput: in silico HiC matrix for the population of structures in input folder, dictionary with communities

    parameters:

    k - is the most important parameter in the umap fuzzy_simplicial_set function.
    It will determine how sparse the final graph will be.
    cores - number of cores available
    number_of_beads_per_structure

    '''

    # list to accumulate matrices
    results = []

    #counter to follow progres
    counter = 0

    #dictionary for str:list of communities
    communities_ditc = {}

    #process multi
    for chunk in chunker(range(1,number_of_structures+1),cores):

        #initiate separate processes for each file in chunk
        ids = [read_and_prepare_graph_and_communities.remote(folder,num,k,number_of_beads_per_structure) for num in chunk]
        # list to accumulate matrices from each file in chunk
        partial_results = []
        # get results from processes
        partial_results_triple = ray.get(ids)
        # wait till all processes are done
        ready, not_ready = ray.wait(ids,num_returns= len(chunk))
        # for each matrix,communities,structure_id
        for triple in partial_results_triple:
            # add entry structure_id : communities to communities dictionary
            communities_ditc[triple[2]] = triple[1]
            # add matrix to partial results
            partial_results.append(triple[0])



        # aggregate matrices to one sigle matrix (to save memory)
        while len(partial_results) > 1:
            partial_results = partial_results[2:] + [aggregate_graphs.remote(partial_results[0], partial_results[1])]
        folder_chunk_results = ray.get(partial_results[0])

        #append aggregated matrix to final results
        results.append(folder_chunk_results)

        #track progress
        counter += cores
        if counter%1000 == 0:
            print(counter)

        # aggregate final matrices (to save memory)
        while len(results) > 1:
            results = results[2:] + [aggregate_graphs.remote(results[0], results[1])]


    return results,communities_ditc
