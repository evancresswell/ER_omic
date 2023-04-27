# STEP 1: Import code
import os.path, sys
from pathlib import Path
# Define path to utility code and insert it into jupyter path for import
path_to_utilities = '/data/cresswellclayec/ER_omic/utilities/'
sys.path.insert(0, path_to_utilities)

# Import computational science suite
import numpy as np
np.random.seed(1) # set random seed for numpy for consistency

import pandas as pd
from scipy import linalg
from scipy.sparse import csr_matrix
from sklearn.preprocessing import OneHotEncoder
from scipy.spatial import distance_matrix

# Import Bipython
import warnings
import Bio.PDB
from Bio import SeqIO, pairwise2,BiopythonWarning
pdb_list = Bio.PDB.PDBList()
pdb_parser = Bio.PDB.PDBParser()
warnings.simplefilter('ignore', BiopythonWarning)

# Import parallel and time-count
from joblib import Parallel, delayed
import timeit
# %matplotlib inline

import matplotlib.pyplot as plt

# # --- Import our Code ---# #
#import emachine as EM
from direct_info import direct_info
# import data processing and general DCA_ER tools
from data_processing import data_processing_msa2pdb
import tools
# # -----------------------# #"
# Which family is being mutated.
family_id = int(sys.argv[1])
n_cpu = int(sys.argv[2])

# Define data directories
line1_dir = '/data/cresswellclayec/ER_omic/energy_landscape/line1'
data_dir = "%s/line1_data" % line1_dir


all_seqs = []
with open('%s/all_seqs.fa' % data_dir, 'rU') as f:
    seq_iter = SeqIO.parse(f,'fasta')
    for seq in seq_iter:
        all_seqs.append(seq)

lp_msas = []
lp_fa_prefix = '%s/' % data_dir
lp_names = ['LPa1', 'LPa2', 'LPa3', 'LPa4', 'LPa5', 'LPa6', 'LPa7' ]
total_len = 0
lp_files = ['cls1.1.fasta','cls1.2.fasta','cls123.3.fasta','cls1ab234.4.fasta','cls123.5.fasta','cls1abc234.6.fasta','cls123.7.fasta']
lp_ids = []
for i,filename in enumerate(lp_files):
    lp_msas.append([])
    lp_ids.append([])
    print('Loading MSA for ',lp_names[i])
    with open(lp_fa_prefix+filename, 'rU') as f:
        seq_iter = SeqIO.parse(f,'fasta')
        for seq in seq_iter:
            lp_msas[-1].append(seq.seq) 
            lp_ids[-1].append(seq.id)
    f.close()
    print(len(lp_msas[-1]))
    total_len += len(lp_msas[-1])
print('number of all individual LPa MSA sequences: ',total_len)
print('number of all aligned sequences: ',len(all_seqs))



non_evo_seqs = []
non_evo_ids = []
filename = '%s/151_cc_all_peps_fnl_ed.1.fa' % data_dir
print('Loading MSA for ',filename)
with open(filename, 'rU') as f:
    seq_iter = SeqIO.parse(f,'fasta')
    for seq in seq_iter:
#             print(seq)

        non_evo_seqs.append(seq.seq) 
        non_evo_ids.append(seq.id)
f.close()

non_evo_seqs = np.array([np.array(list(str(record))) for record in non_evo_seqs])


# match up sequences between individual family fasta files and full fasta files
# We do this instead of directly loading the familiy fasta files 
# to ensure that the family sequences exist in the full sequence set. (sanity check)
#    - print statements show that the files have all corresponding sequences.
family_ref = []
family_indx = []
for ii,ids in enumerate(lp_names):
    family_indx.append([])
for i,seq in enumerate(all_seqs):
    found = False
    for j,ids in enumerate(lp_ids):
        if seq.id in ids:
            family_ref.append(j)
            family_indx[j].append(i)
            found=True
            break
    if not found:
        print('could not categorize sequence!!')
        
family_ref = np.array(family_ref)
print('all familys contained: ',np.unique(family_ref))
# print('\\nnumber of ids attributed to each family:')
for i in range(7):
    family_seqs = family_ref==i
#     print(family_seqs.sum(axis=0))
#     print(len(family_indx[i]))


s0 = np.array([seq.seq for seq in all_seqs])
print(s0.shape)

full_s0 = np.concatenate((s0,non_evo_seqs), axis=0)

onehot_encoder = OneHotEncoder(sparse=False,categories='auto')
onehot_encoder.fit(full_s0)

s = onehot_encoder.transform(s0)
non_evo_s = onehot_encoder.transform(non_evo_seqs)
full_s = onehot_encoder.transform(full_s0)


# number of positions
n_var = full_s0.shape[1]
n_seq = full_s0.shape[0]

print("Number of residue positions:",n_var)
print("Number of sequences:",n_seq)

# number of aminoacids at each position
mx = np.array([len(np.unique(full_s0[:,i])) for i in range(n_var)])
#mx = np.array([m for i in range(n_var)])
print("Number of different amino acids at each position",mx)

mx_cumsum = np.insert(mx.cumsum(),0,0)
i1i2 = np.stack([mx_cumsum[:-1],mx_cumsum[1:]]).T
# print(\"(Sanity Check) Column indices of first and (\",i1i2[0],\") and last (\",i1i2[-1],\") positions\")
# print(\"(Sanity Check) Column indices of second and (\",i1i2[1],\") and second to last (\",i1i2[-2],\") positions\")


# number of variables
mx_sum = mx.sum()
print("Total number of variables",mx_sum)

# number of bias term
n_linear = mx_sum - n_var

from er_energy import w_seq_dist           # function which will mutate sequence
                                          # random seed is set in er_energy.
import random    
nwalk = 100
niter = 1000


# define files to save results to
w_file = "%s/%s_w.npy" % ( data_dir, lp_names[family_id])   
wsym_file = "%s/%s_w_sym.npy" % (data_dir, lp_names[family_id])   
b_file = "%s/%s_b.npy" % (data_dir, lp_names[family_id])        
 
# want to use asymmetric w for mutating sequence
w_ER = np.load(w_file)      
w_ER_sym = np.load(w_file)                                                                               
b = np.load(b_file)                        

# mean of family sequences
# this sequence also contains information on which aa are expressed in family
# the aa-position information for the family is used/assumed in w_seq_dist                                                                                        
gp_mean = np.mean(s[family_indx[family_id]], axis=0)
# get the sequience closest to the mean.
for i,seq in enumerate(s[family_indx[family_id]]):
    dist = np.linalg.norm(gp_mean-seq)
    if i == 0:
        min_indx = i
        min_dist = dist
    elif dist < min_dist:
        min_dist = dist
        min_indx = i
print('sequence %d closets to group mean with distance of %f' % (min_indx, min_dist))
mean_seq = s[family_indx[family_id]][min_indx] 

print(data_dir)
print(lp_names[family_id])
print(nwalk)
print(niter)
# depending on how many cores you have you may need to change n_jobs and ncpu...

#w_seq_dist(mean_seq, beta, i1i2, w_ER, b, s[family_indx[family_id]], n_iter=niter,seed=42, ncpu = 2)

betas = [10e-4, 10e-3, .1, 1., 10., 100., 1000., 10000.]
for i, beta in enumerate(betas):
    print('running sequence walk with beta = ', beta)

    start_time = timeit.default_timer()                                                       
    res_walk = Parallel(n_jobs = int(n_cpu/2))(delayed(w_seq_dist) 
            (mean_seq, beta, i1i2, w_ER, b, s[family_indx[family_id]], n_iter=niter,seed=i0, ncpu = 2)
            for i0 in range(nwalk)) 


    run_time = timeit.default_timer() - start_time   
    print('run time:',run_time) 

    sequence_walk_file = "%s/%s_dist_%dwalk_%di_beta%d.npy" % (data_dir, lp_names[family_id], nwalk, niter,i)
    print('saving to sequence walk to %s' % sequence_walk_file)

   
    SW_dist = []
    for walk in res_walk:
        SW_dist.append(walk)
    
    np.save(sequence_walk_file, SW_dist)
    

