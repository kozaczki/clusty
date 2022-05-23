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
    """helper funciton for dividing
    sequence into chunks of equal size

    seq - list to divide
    size - size of chunk

    """
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def structure_file_to_distances_matrix(structure_file):

    """
    loads csv file with coordinates(x,y,z) in first three columnd
    returns matrix with distances

    structure_file - string with path to csv file

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

def get_average_distance_single_structure(borders,distances):

    """
    calculates the average distance between two neighbour beads in single structure - this value can differ dramatically across
    single structure as to different chromatin condensation

    borders - 1D array or list with first and last indicies of chromosomes coming in order 1,2,3...22,sex chromosome
    distances - distance matrix (square form)

    """

    # holders for accumulating distance and number of beads
    counter = 0
    cumulative_distance = 0

    # iterate over chromosome

    for i in range(23):
        first = borders[i*2] - 1
        last = borders[(i*2)+1] - 1
        # iterate over beads in chromosome (form 1 to last-1)
        for j in range(first,last):

            distance = distances[j,j+1]
            cumulative_distance += distance
            counter += 1

    # calculate mean

    mean_chromosomal_distance = cumulative_distance/counter

    return mean_chromosomal_distance

def get_average_distance(dataset_folder,n,borders_array):

    """
    calculates the average distance between two neighbour beads in subset of structes

    dataset_folder - location of structure files
    borders - 1D array or list with first and last indicies of chromosomes coming in order 1,2,3...22,sex chromosome
    n - number of structres to be used for calculating the average

    """

    # raed dataset folder content
    files = os.listdir(dataset_folder)

    #variable to accumulate sum
    average_distances_sum = 0

    #process n files from dataset folder

    for file in files[:n]:

        #build path for file
        structure_file = os.path.join(dataset_folder,file)

        #prepare distance matrix for the file(structure)
        distances = structure_file_to_distances_matrix(structure_file)

        #calculate avergae distance for the file (one structure)
        average_distance = get_average_distance_single_structure(borders_array,distances)

        #add obtain average structue to accumulating result
        average_distances_sum = average_distances_sum + average_distance

    #calculate and return average distance in n structures
    average_distance = average_distances_sum / n
    return average_distance

def get_max(sparse_matrices_list):
    maximum = 0
    for matrix in sparse_matrices_list:
        if matrix[0].shape[1] > maximum:
            maximum = matrix[0].shape[1]
    return maximum

@ray.remote
def structure_to_neigbours(structure_file,radius_c,num):

    """
    takes a csv file of a single structure and return a matrix with neighbours within radius
    as sparse matrix, together with csv id number

    """



    # distance matrix
    dist_mat = structure_file_to_distances_matrix(structure_file)
    distances_matrix_sorted = np.sort(dist_mat, axis = 1)

    # distance matrix args + 1 -> fbead numbers
    distances_matrix_sorted_args = np.argsort(dist_mat, axis = 1) + 1 # sort gives indicies of beads so to get ids -> + 1

    # get number of beads within treshold per bead
    numbers_of_beads_within_treshold = (distances_matrix_sorted <= radius_c).sum(axis = 1)


    # prepare matrix for neigbours
    first_dimension = dist_mat.shape[0] # number of beads in structure - could be replaced here
    second_dimension = numbers_of_beads_within_treshold.max() - 1 #because first arg is the origin bead
    neighbours_in_structure = np.zeros((first_dimension,second_dimension),dtype = np.ushort)

    for i in range(first_dimension):
            beads_within_treshold = numbers_of_beads_within_treshold[i]
            bead_neighbours = distances_matrix_sorted_args[i,1:beads_within_treshold] #start at 1 as 0 is the origin bead (self)
            neighbours_in_structure[i] = np.pad(bead_neighbours,(0,second_dimension - beads_within_treshold + 1),constant_values = 0) #adjust shape to second dimension
            sparse_neighbours = sparse.csr_matrix(neighbours_in_structure)
    return sparse_neighbours,num

def prepare_file_name(number):
    name = 'cf_' + number.zfill(6) + '.coords.csv'
    return name

def get_max(sparse_matrices_list):
    maximum = 0
    for matrix in sparse_matrices_list:
        if matrix[0].shape[1] > maximum:
            maximum = matrix[0].shape[1]
    return maximum
