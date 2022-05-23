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
import sys

import functions


path_to_parameters = sys.argv[1]

# load paramaters from csv file

# parse csv file with parameters
paramaters = []
with open(path_to_parameters, 'rt') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    paramaters.append(list(reader))
    csvfile.close()

#list with setup parameters
params = paramaters[0]

#assign setup variebles from params

home = params[0][1]
number_of_structures = int(params[1][1])
number_of_beads_per_structure = int(params[2][1])
fraction = float(params[3][1])
structures_fraction = number_of_structures * fraction
cores = int(params[4][1])
dataset_name =  params[6][1]
dataset_folder =  params[7][1]
a_type = params[8][1]
chromosomal_borders_file = params[9][1]
primes_file = params[10][1]
clustering_r = float(params[11][1])

# # compose analysis name
#
# analysis_name = dataset_name + '_neighbours_' + str(k)

# handle dataset_name depending on analysis_type

if a_type == 'fixed':
    r_factor  = float(params[5][1])
    analysis_name = dataset_name + '_fixed_radius_' + str(r_factor)

if a_type == 'neighbours':
    k  = int(params[5][1])
    analysis_name = dataset_name + '_neighbours_' + str(k)

# print setup variables for manual inspection
print("")
print("Running cluster detection in structures")
print("")

print("dataset name: " + dataset_name)

print("loaded setup variables")
print("")
print("home folder: " + home)
print("dataset folder: " + dataset_folder)
print("dataset name: " + dataset_name)
print("number of structures: " + str(number_of_structures))
print("number of beads per structure: " + str(number_of_beads_per_structure))
print("fraction: " + str(fraction))

if a_type == "fixed":
    print("radius factor: " + str(r_factor))
if a_type == 'neighbours':
    print("k: " + str(k))

print("cores: " + str(cores))
print("")
print("clustering_r: " + str(clustering_r))

# PATHS

helper_folder = os.path.join(home,'helper_data')

primes_array = np.load(os.path.join(helper_folder,primes_file),allow_pickle=True)

# LOAD clusters

clusters_path = os.path.join(home,'runs',analysis_name,'results' , analysis_name + '_clusters_simple_filtered')
clusters_list = functions.read_from_pickle(clusters_path)

# load product

product_full_path = os.path.join(home,'runs',analysis_name,'intermediate', 'products','product_full.npy')
product_full = np.load(product_full_path,allow_pickle=True)

# identify cluster

start_distributing = time.time()

print("distributing clusters in structures\n")

distro = functions.get_distributions_for_population(clusters_list,product_full,primes_array)

print(str(time.time() - start_distributing))

fig,ax = plt.subplots(figsize = (5,5))
ax.hist(distro.sum(axis=1))
fig.savefig(os.path.join(home,'runs',analysis_name,'figures',analysis_name + '_distro_hist.png'))

distro_full_path = os.path.join(home,'runs',analysis_name,'intermediate',analysis_name + '_distro_matrix')
np.save(distro_full_path,distro)

# PART II - identify clusters in structures:

# load list
list_path = os.path.join(home,'runs',analysis_name,'intermediate',analysis_name + '_sparse_neighbour_matrices_list')
list_of_str = functions.read_from_pickle(list_path)

# change list_of_str to lighter
list_str = [i[1] for i in list_of_str]


# loop through matrix and structures:

# begin multiprocessing

ray.init()

start = time.time()

cutoff = clustering_r

results = []

path_to_store = os.path.join(home,'runs',analysis_name,'intermediate',analysis_name + '_clusters_in_structures')

# process bead by bead

#for chunk in chunker(range(number_of_beads_per_structure),cores):

print("detecting clusters in individual structures")

for chunk in functions.chunker(range(number_of_structures),cores):
    communities = [functions.build_communities.remote(list_str[i],distro[i],cutoff,dataset_folder) for i in chunk]
    partial_results = ray.get(communities)
    results.append(partial_results)

    file_to_store = open(path_to_store, "wb")
    pickle.dump(results, file_to_store)

    file_to_store.close()


print(time.time() - start)
