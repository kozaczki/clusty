# import libraries

import sys
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

import functions



print("")
print("******************************* in silico HiC *******************************")
print("")

# load path to csv with paramaters from command line
path_to_parameters = sys.argv[1]

print("-------------> loading parameters from: " + path_to_parameters)
print("")

print("paramaters\n")

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
cores = int(params[3][1])
k = int(params[4][1])
dataset_name =  params[5][1]
dataset_folder =  params[6][1]

# print loaded parameters:

print('Home directory : ' + str(home))
print('Dataset name :' + dataset_name)
print('Dataset directory: ' + dataset_folder)
print('number of structures: ' + str(number_of_structures))
print('number of beads per structure: ' + str(number_of_beads_per_structure))
print('cores: ' + str(cores))
print('k for graph: ' + str(k))

# compose analysis name

analysis_name = dataset_name + '_inSilico_' + str(k)

# build folders structure

run_folder = os.path.join(home,'runs',analysis_name)
results_folder = os.path.join(run_folder,'results')
figures_folder = os.path.join(run_folder,'figures')

print("")
print("-------------> building folders structure for the run\n")

folders = [run_folder,results_folder,figures_folder]

for folder in folders:
    try:
        os.mkdir(folder)
        print("Directory " , folder ,  " Created ")
    except FileExistsError:
        print("Directory " , folder ,  " already exists")

print("")
print("-------------> initializing ray\n")
# initialize ray

ray.init()

print("")
print("-------------> calculating in silico HiC matrix\n")

# run algorith on chosen dataset

# process all structures

inSilicoHiC,communities_dict = functions.process_csvs(dataset_folder,cores,k,number_of_beads_per_structure,number_of_structures)

# get results

out = ray.get(inSilicoHiC[0])

# visualize matrix
print("")
print("-------------> Saving results\n")

# transform hic matrix to dense form
PA_dense = sparse.csr_matrix.todense(out)
# set up visualization theme
sns.set_theme()
# figsize in inches
fig, ax = plt.subplots(figsize=(20,15))
ax = sns.heatmap(PA_dense)
# set up figure title
title = "in silico HiC for " + dataset_name + " k=" + str(k)
fig.suptitle(title)
# save figure
fig.savefig(os.path.join(figures_folder,'hic_' + analysis_name))

print("in silico HiC matrix visualization saved as: " + str(os.path.join(figures_folder,'hic_' + analysis_name)) +"\n")

# save matrix
# set up path to matrix
file_to_store__HIC = os.path.join(results_folder,'HIC_matrix_' + analysis_name)
# save file
sparse.save_npz(file_to_store__HIC,out)
print("in silico HiC matrix values saved as: " + str(os.path.join(results_folder,'HIC_matrix_' + analysis_name)) +"\n")


# save communities
# set up path to communities
file_to_store_communities = open(os.path.join(results_folder,'communities_' + analysis_name),'wb')
# save and close the file
pickle.dump(communities_dict,file_to_store_communities)
file_to_store_communities.close()

print("communities saved as: " + str(os.path.join(results_folder,'communities_' + analysis_name)))
