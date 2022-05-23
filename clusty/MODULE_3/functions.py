# load libraries

import numpy as np
import umap
import umap.umap_ as umap

import cdlib
from cdlib import algorithms
import networkx as nx


import warnings
warnings.filterwarnings('ignore')

import time
import seaborn as sns;
import matplotlib.pyplot as plt


import ray
from scipy.sparse import csr_matrix
import csv

import pickle

import scipy.spatial
from scipy import sparse
import itertools
import os

def read_from_pickle(path):
    file_to_read = open(path, "rb")
    loaded_object = pickle.load(file_to_read)
    file_to_read.close()
    return loaded_object


def get_indicies_for_cluster(cluster,product,primes):
    """
    gets structures in which given cluster occurs
    """
    # turn cluster beads to indicies
    cluster_indexes = [i - 1 for i in cluster]
    # turn indicies to primes product
    cluster_primes_product = np.prod(primes[cluster_indexes])
    # check occurence
    indicies = (np.where(product % cluster_primes_product == 0)[0])

    return indicies

def get_distributions_for_population(clusters_list,product_full,primes):
    clusters_in_structures = np.zeros((number_of_structures,number_of_beads_per_structure))
    for pair in clusters_list:
        clusters = pair[1]
        if len(clusters) > 0:
            bead = pair[0]
            bead_index = bead - 1
            for cl in clusters:
                str_indicies = get_indicies_for_cluster(cl,product_full[:,bead_index],primes)
                for str_index in str_indicies:
                    for b in cl:
                        clusters_in_structures[str_index,b-1] = 1

    return clusters_in_structures

# helper funcitons

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

@ray.remote
def build_communities(num,clusters_array,cutoff,dataset_folder):
    file = os.path.join(dataset_folder,'cf_' + str(num).zfill(6) + '.coords.csv')
    csv = np.genfromtxt(file,delimiter=',')
    coordinates = csv[:,:3]
    clustered = np.where(clusters_array == 1)[0]
    if len(clustered) == 0:
        return num, []
    coords_for_communities = coordinates[clustered]
    distances = scipy.spatial.distance.pdist(coords_for_communities)
    distances_matrix = scipy.spatial.distance.squareform(distances)
    adj_matrix = np.zeros((len(clustered),len(clustered)))
    for i in range(len(clustered)):
        for j in range(len(clustered)):
            if distances_matrix[i,j] <= cutoff:
                adj_matrix[i,j] = 1
    g = nx.from_numpy_matrix(adj_matrix)
    eset = [(u, v) for (u, v, d) in g.edges(data=True)] # get list of edges from graph
    weights = [d['weight'] for (u, v, d) in g.edges(data=True)] # get list of weights from edges
    # find communities
    # in this example we use the Leiden algorithm
    leiden_coms = algorithms.leiden(g,weights=weights) # check if the algo is stochastic, in that case set rnd generator
    leiden_coms.communities # a list of lists of nodes
    communities_in_str = []
    for community in leiden_coms.communities:
        communities_in_str.append(clustered[community])
    return num,communities_in_str
