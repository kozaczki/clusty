# load libraries

import numpy as np
import scipy.spatial
import ray
import os
import time
import re
import pickle
from scipy import sparse
import time
import matplotlib.pyplot as plt
import sys
import csv
import itertools
from itertools import combinations

# helper funcitons

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def get_beads_of_interest(bead_number,neigbours_array,borders_array,structures_fraction,number_of_beads_per_structure,primes_array):

    """
    takes a bead and returns indicies of beads that are in the neighbourhood of this bead in the desired fraction of structures

    """

    bead_index = bead_number - 1


    # read chromosomal borders - first and last bead for given chromosome
    chrom_limits = borders_array[bead_index]
    # bead = 0 -> [1,249]
    # bead NUMBERS

    # change from ids to indexes -> therefore -1
    first_bead_of_chromosome_index = chrom_limits[0] - 1 # beacuse it  reads x-1
    last_bead_of_chromosome = chrom_limits[1] - 1

    # pick row
    bead_of_interest = neigbours_array[bead_index]

    # get indicies of contacts above treshold
    bead_of_interest_indicies = np.where(bead_of_interest >= structures_fraction)[0] #structes fraction because it's frequencies array

    # filter for valus under these indicies


    right_slice =  np.where(bead_of_interest_indicies < first_bead_of_chromosome_index)[0]
    left_slice = np.where(bead_of_interest_indicies >last_bead_of_chromosome )[0]
    left_and_right = np.append(right_slice,left_slice)

    indicies = bead_of_interest_indicies[left_and_right]

    # RETURNS INDICES OF BEADS


    return indicies

@ray.remote
def pick_fours(bead_number,neighbours_array,borders,product_folder,structures_fraction,number_of_beads_per_structure,primes_array):

    # bead numbers         [1,2,3,4,5]
    # primes_array indexes [0,1,2,3,4]
    # primes               [2,3,5,7,11]

    # load a product for given bead (product was saved under bead id)

    prod = np.load(product_folder + '/prod_' + str(bead_number) +'.npy',allow_pickle=True)

    # add triplets here

    bead_triplets = []

    bead_index = bead_number - 1


    # get all indicies of bead in trans chromosomes
    bead_x_indicies = get_beads_of_interest(bead_number,neighbours_array,borders,structures_fraction,number_of_beads_per_structure,primes_array) # indices
    #get all possible doublets (in reality they are triplets, because of bead of origin)
    bead_0_indicies_all = [i for i in range(number_of_beads_per_structure) if (neighbours_array[bead_index,i] >= structures_fraction and i != bead_index)]

    # combine cis and trans
    bead_0_indicies_all_combined = list(combinations(bead_0_indicies_all,2))

    # filter combined for possible triplets
    bead_0_indicies_all_filtered = [(i[0],i[1]) for i in bead_0_indicies_all_combined if neighbours_array[i[0],i[1]] >= structures_fraction ]

    bead_0_indicies_all_filtered_subset = bead_0_indicies_all_filtered

    # bead numbers         [1,2,3,4,5]
    # primes_array indexes [0,1,2,3,4]
    # primes               [2,3,5,7,11]



    # neighbours for bead
    bead_0_str = prod
    primes_inter = primes_array[bead_x_indicies] # trans beads
    for i in bead_0_indicies_all_filtered_subset: # all doublets that passed >= fraction
        product = primes_array[i[0]]*primes_array[i[1]] # doublet as primes
        mods = np.mod(bead_0_str,product) # array with remainders
        ins = np.where(mods == 0)[0] # indicies of structures in which they are present
        ins_sum = ins.shape
        if ins_sum[0] >= structures_fraction: # if for given doublet there is more than hundred structures

        # subset of structures where doublet is present

            bead_0_str_subset = bead_0_str[ins]    # subset of strucutres to check

            potential_triplets_product = primes_inter * product  # THAT'S WHERE TRANS ARE ADDED

            potential_triplets_product_full = np.full((bead_0_str_subset.shape[0],potential_triplets_product.shape[0]),potential_triplets_product)
            mods_2 = np.mod(bead_0_str_subset,potential_triplets_product_full.T) # remainders for triplets

       #mods_2_reshape = np.reshape(mods_2(products_2.shape[0],))

            ins_2_0 = np.where(mods_2 == 0)[0] # by 3rd element [indexes of ]
            ins_2_1 = np.where(mods_2 == 0)[1] # by structure number

            val,counts = np.unique(ins_2_0,return_counts=True)


            counts_more_than_hundres = np.where(counts >= structures_fraction)

        #print("counts" + str(counts_more_than_hundres[0]))
            indexes_of_primes_sybset_passing_the_treshold = val[counts_more_than_hundres[0]]

        #print(indexes_of_primes_sybset_passing_the_treshold)
            pr = primes_inter[indexes_of_primes_sybset_passing_the_treshold]
        #print(pr)
        #print(i)
            for p in pr:
                t = (int(np.where(primes_array==p)[0]))
                if (i[0] != i[1] and i[0] != t and i[1] != t): # retruns index of prime, to get bead: add 1
                    #bead_triplet = [i[0],i[1],t]
                    bead_triplet = [i[0]+1,i[1]+1,t+1] # returns beads IDS

                    bead_triplet_sorted = tuple(sorted(bead_triplet))
                    bead_triplets.append(bead_triplet_sorted)


    if len(bead_triplets) > 0:
        print(bead_number,len(bead_triplets))


    return (bead_number, bead_triplets)

def filter_clusters(cluster):
    filtered_clusteres = []
    for triplet in cluster[1]:
        triplet_list = list(triplet)
        triplet_list_sorted = sorted(triplet_list)
        if triplet_list_sorted not in filtered_clusteres:
            filtered_clusteres.append(triplet_list_sorted)
    return cluster[0],filtered_clusteres
