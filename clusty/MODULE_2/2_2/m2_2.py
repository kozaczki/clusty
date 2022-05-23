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

import functions

print("")
print("********************* population-wide cluster detection *********************")
print("")



#### ENTER PATH HERE ####

path_to_parameters = sys.argv[1]

print("-------------> loading parameters from: " + path_to_parameters)
print("")


# load paramaters from csv file

paramaters = []
with open(path_to_parameters, 'rt') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    paramaters.append(list(reader))
    csvfile.close()

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

# handle dataset_name depending on analysis_type

if a_type == 'fixed':
    r_factor  = float(params[5][1])
    analysis_name = dataset_name + '_fixed_radius_' + str(r_factor)




if a_type == 'neighbours':
    k  = int(params[5][1])
    analysis_name = dataset_name + '_neighbours_' + str(k)

# print setup variables for manual inspection
print("")

print("analysis name: " + analysis_name)
print("dataset name: " + dataset_name + '\n')


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

# defin regex patter for csv files
pattern = '^(.*)cf_(.*).coords.csv$'

run_folder = os.path.join(home,'runs',analysis_name)
intermediate_files_folder = os.path.join(run_folder,'intermediate')
product_folder = os.path.join(intermediate_files_folder,"products")
results_folder = os.path.join(run_folder,'results')
figures_folder = os.path.join(run_folder,'figures')

# load stored input data

#  SAVE FILES
binary = np.load(os.path.join(intermediate_files_folder,analysis_name + '_binary.npy'))
frequencies = np.load(os.path.join(intermediate_files_folder,analysis_name + '_frequencies.npy'))

# HELPER DATA FILES

helper_folder = os.path.join(home,'helper_data')



# load and process helper data

primes_array = np.load(os.path.join(helper_folder,primes_file),allow_pickle=True)

# CHROMOSOMAL BORDERS ARRAY

chromosomal_borders = np.load(os.path.join(helper_folder,chromosomal_borders_file))

# get chromosomal borders

borders_array = np.unique(chromosomal_borders)

# begin multiprocessing

print("Clusters identification\n")

ray.init()

start = time.time()

results = []

path_to_store = intermediate_files_folder + analysis_name +  "_final_clusters_running"

# process bead by bead

#for chunk in chunker(range(number_of_beads_per_structure),cores):
for chunk in functions.chunker(range(number_of_beads_per_structure),cores):

    ids = [functions.pick_fours.remote((bi+1),frequencies,chromosomal_borders,product_folder,structures_fraction,number_of_beads_per_structure,primes_array) for bi in chunk]
    partial_results = ray.get(ids)
    results.append(partial_results)

    file_to_store = open(path_to_store, "wb")
    pickle.dump(results, file_to_store)

    file_to_store.close()

print("cluster identification took " + str(time.time() - start) + " s\n")

# put all the lists together

import itertools
clusters_list = results[0]
for i in range(1,len(results)):
    clusters_list = list(itertools.chain(clusters_list,results[i]))

filtered_clusters = [functions.filter_clusters(i) for i in clusters_list]

filtered_clusters_list_to_save_path = os.path.join(results_folder,analysis_name + '_clusters_simple_filtered')
file_to_store = open(filtered_clusters_list_to_save_path, "wb")
pickle.dump(filtered_clusters, file_to_store)
file_to_store.close()

print("clusters saves as : " + filtered_clusters_list_to_save_path + '\n')



clusters_list_to_save = list(itertools.chain(*results[0]))
clusters_list_to_save_path = os.path.join(results_folder,analysis_name + '_clusters_simple')
file_to_store = open(clusters_list_to_save_path, "wb")
pickle.dump(clusters_list_to_save, file_to_store)
file_to_store.close()

# reformat data for plot with clusters distribution
beads_numbers = [i[0] for i in clusters_list]
cluster_numbers = [len(i[1]) for i in clusters_list]

# reformat data for plot with FILTERED clusters distribution
beads_numbers_filtered = [i[0] for i in filtered_clusters]
cluster_numbers_filtered = [len(i[1]) for i in filtered_clusters]

# prepare plot data for plot with clusters distribution
fig, ax = plt.subplots(2,1,figsize = (20,5))
ax[0].scatter(beads_numbers,cluster_numbers)
ax[1].scatter(beads_numbers_filtered,cluster_numbers_filtered)


fig.savefig(os.path.join(figures_folder,analysis_name + '_clusters_distribution.png'))

print("Clusters distribution figure saved as : " + str(os.path.join(figures_folder,analysis_name + '_clusters_distribution.png')))


# save final clusters

path_to_store = intermediate_files_folder + analysis_name +  "final_clusters"


final_clusters = (results,beads_numbers,cluster_numbers)

file_to_store = open(path_to_store, "wb")
pickle.dump(final_clusters, file_to_store)

file_to_store.close()

print("done")
