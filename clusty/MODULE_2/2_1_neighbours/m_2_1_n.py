import functions

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

print("")
print("********************* preprocessing structures for the nearrest neighbours neighbourhood *********************")
print("")

path_to_parameters = sys.argv[1]

# start time to trace execution time

print("-------------> loading parameters from: " + path_to_parameters)
print("")

print("paramaters\n")

start_all = time.time()

# load paramaters from csv file

# parse csv file with parameters
paramaters = []
with open(path_to_parameters, 'rt') as csvfile:
    reader = csv.reader(csvfile, skipinitialspace=True)
    paramaters.append(list(reader))
    csvfile.close()

#list with setup parameters
params = paramaters[0]

# check if correct mode

mode = params[8][1]

if mode == 'fixed':
    print("incorrect mode - use either neighbours mode or 3_1_fixed")

#assign setup variebles from params

home = params[0][1]
number_of_structures = int(params[1][1])
number_of_beads_per_structure = int(params[2][1])
fraction = float(params[3][1])
structures_fraction = number_of_structures * fraction
cores = int(params[4][1])
k  = int(params[5][1])
dataset_name =  params[6][1]
dataset_folder =  params[7][1]
borders_file = params[9][1]
primes_file = params[10][1]


# compose analysis name

analysis_name = dataset_name + '_neighbours_' + str(k)

# print setup variables for manual inspection

print("loaded setup variables")
print("")
print("dataset name: " + dataset_name + '\n')
print("analysis name: " + analysis_name + '\n')

print("home folder: " + home)
print("dataset folder: " + dataset_folder)
print("number of structures: " + str(number_of_structures))
print("number of beads per structure: " + str(number_of_beads_per_structure))
print("fraction: " + str(fraction))
print("k: " + str(k))
print("cores: " + str(cores))
print("primes file: " + str(primes_file))
print("borders file: " + str(borders_file))

print("")

# build paths for helper data files

helper_folder = os.path.join(home,'helper_data')

# load and process helper data

    # primes array

primes_array = np.load(os.path.join(helper_folder,primes_file),allow_pickle=True)

    # chromosomal borders array

chromosomal_borders = np.load(os.path.join(helper_folder,borders_file),allow_pickle=True)

    # transform borders matrix into 1D

borders_array = np.unique(chromosomal_borders)

# build folders structure

run_folder = os.path.join(home,'runs',analysis_name)
intermediate_files_folder = os.path.join(run_folder,'intermediate')
product_folder = os.path.join(intermediate_files_folder,"products")
results_folder = os.path.join(run_folder,'results')
figures_folder = os.path.join(run_folder,'figures')

print("building folders structure for the run\n")


folders = [run_folder,intermediate_files_folder,product_folder,results_folder,figures_folder]

for folder in folders:
    try:
        os.mkdir(folder)
        print("Directory " , folder ,  " Created ")
    except FileExistsError:
        print("Directory " , folder ,  " already exists")

print("")



print("extracting neighbours\n")

# start multiprocessing

ray.init()

print("")

# process structures in batch:
# csv -> distances -> k neighbours
# and prepare list of tuples (neighbours matrix,id of structure)
# should be ordered by ids

files = os.listdir(dataset_folder)

sparse_matrices = []

for numbers_chunk in functions.chunker(list(range(1,(number_of_structures+1))),cores):
    ids = [functions.structure_to_neigbours.remote(os.path.join(dataset_folder,functions.prepare_file_name(str(number))),k,number) for number in numbers_chunk]
    partial_results = ray.get(ids)
    sparse_matrices.append(partial_results)

print("extracting neighbours finished\n")

# put all matrices in one list

final_list = sparse_matrices[0]
for i in range(1,len(sparse_matrices)):
    final_list = list(itertools.chain(final_list,sparse_matrices[i]))

# print length of the matrices list

print("Final list of neighbourhoods contains " + str(len(final_list)) + " elements\n")

print(time.time() - start_all)


# save matrices -> intermediate/dataset_name/_sparse_neighbour_matrices_list

print("saving neighbours matrices list \n")

file_to_store = open(os.path.join(intermediate_files_folder,analysis_name + '_sparse_neighbour_matrices_list') , "wb")
pickle.dump(final_list, file_to_store)
file_to_store.close()

print("saving neighbours matrices list finished\n")

# build matrix with neighbourhoods - all beads across all structures -> to be used for products
# values in  neigbhours correspond to beads ids
# keeps ids

print("saving neighbours matrix \n")

neigbhours = np.zeros((number_of_structures,number_of_beads_per_structure,k),dtype = int)

# iterate over list of neigbours in structures
for x in range(number_of_structures):

    # pick neighbours from list
    structure_x = final_list[x][0]
    # convert to dense form
    structure_x_dense = sparse.csc_matrix.todense(structure_x)
    # put in the neighbours matrix
    neigbhours[x] = structure_x_dense

print("saving neighbours matrix finished \n")

# here goes to indexes

print("saving frequency matrix \n")

frequencies = np.zeros((number_of_beads_per_structure,number_of_beads_per_structure),dtype = int)

# prepare coocurence matrix based on neigbhours matrix

#iterate over beads dimension in neighbourhood matrix
for x in range(number_of_beads_per_structure):

    # get neighbourhood for bead x across all structures
    bead_x_all_structures = neigbhours[:,x,:]
    # get occurences of beads in bead x neighbourhood across all structures
    val, counts = np.unique(bead_x_all_structures,return_counts=True)

    # assign values to frequencies matrix
    for i in range(len(val)):
        if val[i] != 0:
            frequencies[x,val[i]-1] = counts[i] # i - 1 becuese in neigbhours matrix values are correspondingto bead ids, whereas in frequencies it goes by indexes


print("saving frequency matrix finished\n")

# binary matrix <- used in module 3_2 to identify clusters, generate from frequencies matrix, if values is above treshold binary matirx get 1, otherwise 0

print("saving binary matrix \n")

binary = np.zeros((number_of_beads_per_structure,number_of_beads_per_structure),dtype = int)
for x in range(number_of_beads_per_structure):
    for y in range(number_of_beads_per_structure):
        if frequencies[x,y] >= structures_fraction:
            binary[x,y]=1

print("saving binary finished \n")

print("converting neighbourhoods to products of primes\n")

# get product

primes = primes_array

# bead numbers         [1,2,3,4,5]
# primes_array indexes [0,1,2,3,4]
# primes               [2,3,5,7,11]


# matrix to hold products

products = np.zeros((number_of_structures,number_of_beads_per_structure),dtype = object)


# iterate over structures

for x in range(len(final_list)): # structures

    # load structure x
    structure_x = final_list[x][0]
    # convert from sparse to dense
    structure_x_dense = sparse.csc_matrix.todense(structure_x)
    # get shape for conversion to primes
    shape = np.shape(structure_x_dense)
    second_dimension = shape[1]
    # matirx for primes
    beads_as_primes = np.zeros((number_of_beads_per_structure,second_dimension),dtype = object)

    # for each bead
    for i in range(number_of_beads_per_structure):
        # for each of the neighbours
        for j in range(second_dimension):
            bead = structure_x_dense[i,j]
            # if 0 than no neigbour
            if bead == 0:
                # 1 multiplying by 1 will not change product
                beads_as_primes[i,j] = 1
            else:
                beads_as_primes[i,j] = primes[bead - 1] # as below
                # bead numbers         [1,2,3,4,5]
                # primes_array indexes [0,1,2,3,4]
                # primes               [2,3,5,7,11]

    # obtain neighbourhoods as products
    product_str_x = np.prod(beads_as_primes,axis = 1)
    # assign to products matrix
    products[x] = product_str_x

# save full product

print("converting neighbourhoods to products of primes finished\n")

print("saving products of primes\n")

np.save(product_folder + '/product_full', products)

print("Full product saved as: " + product_folder + '/product_full\n')

# save products by beads (for module 3_2)
for b in range(number_of_beads_per_structure):
    prod_bead = products[:,b]
    np.save(product_folder + '/prod_' + str(b+1),prod_bead)

print("Products per bead saved in : " + product_folder + '\n')

print('saving figures and files\n')

# prepare figure for frequencies

sns.set_theme()
fig, ax = plt.subplots(figsize=(20,15))
ax = sns.heatmap(frequencies)
fig.savefig(os.path.join(figures_folder,analysis_name + '_frequencies.png'))

print("Frequencies matrix visualisation saved as : " + str(os.path.join(figures_folder,analysis_name + '_frequencies.png')) + '\n')

# prepare figure for binaryy

sns.set_theme()
fig, ax = plt.subplots(figsize=(20,15))
ax = sns.heatmap(binary)
fig.savefig(os.path.join(figures_folder,analysis_name + '_binary.png'))

print("Binary matrix visualisation saved as : " + str(os.path.join(figures_folder,analysis_name + '_binary.png')) + '\n')

#  SAVE FILES frequencies and binaries

np.save(os.path.join(intermediate_files_folder,analysis_name + '_binary'),binary)
np.save(os.path.join(intermediate_files_folder,analysis_name + '_frequencies'),frequencies)

print("Frequencies matrix values saved as : " + str(os.path.join(intermediate_files_folder,analysis_name + '_frequencies')) + '\n')
print("Binary matrix values saved as : " + str(os.path.join(intermediate_files_folder,analysis_name + '_binary')) + '\n')


print("done\n")
