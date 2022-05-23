import numpy as np
import scipy.spatial
from scipy import sparse
import ray
import os
import time
import re
import pickle
import itertools
import seaborn as sns
import matplotlib.pyplot as plt
import csv
import sys

def chunker(seq, size):
    """
    input: sequence (eg. list)
    ouptput: sequence divided in chunks
    parameters: size: determines the length of the chunk

    """
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def structure_file_to_distances_matrix(structure_file):

    """
    input: csv file with coordinates
    output: distance matrix
    """

    # load csv
    csv = np.genfromtxt(structure_file,delimiter=',')

    # limit to first theree columns
    coordinates = csv[:,:3]

    # calculate all distances
    distances = scipy.spatial.distance.pdist(coordinates)

    # convert distance matrix to square form
    distances_matrix = scipy.spatial.distance.squareform(distances)

    return distances_matrix

@ray.remote
def structure_to_neigbours(structure_file,neighbours,number):

    """
    takes a csv file of a single structure and returns
    matrix of k-nearrest beads for each bead as sparse matrix, together with csv id number

    """
    
    # read csv and return a matrix with pairwise distances

    dist_mat = structure_file_to_distances_matrix(structure_file)

    # sort according to distances

    distances_matrix_sorted = np.sort(dist_mat, axis = 1)

    # argsort according to distances, returns indexes, + 1 to obtain beads ids

    distances_matrix_sorted_args = np.argsort(dist_mat, axis = 1) + 1

    # first dimension should be equal to number of beads in structure

    first_dimension = dist_mat.shape[0]

    # second dimension is the k parameter

    second_dimension = neighbours

    # trim sorted beads ids -> [all beads, first neighbour:k-th neighbours]

    neighbours = distances_matrix_sorted_args[:,1:second_dimension+1]

    # convert to ushort - save mem

    neighbours = neighbours.astype(np.ushort)

    # covert to sparse form

    sparse_neighbours = sparse.csr_matrix(neighbours)

    return sparse_neighbours,number

def prepare_file_name(number):
    name = 'cf_' + number.zfill(6) + '.coords.csv'
    return name
