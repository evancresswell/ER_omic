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
# Define data directories
line1_dir = '/data/cresswellclayec/ER_omic/energy_landscape/line1'
data_dir = "%s/line1_data" % line1_dir

n_cpu = int(sys.argv[1])
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


family_ref = []

for i,seq in enumerate(all_seqs):
    found = False
    for j,ids in enumerate(lp_ids):
        if seq.id in ids:
            family_ref.append(j)
            found=True
            break
    if not found:
        print('could not categorize sequence!!')
        
family_ref = np.array(family_ref)
print('all familys contained?: ',np.unique(family_ref))
print('\nnumber of ids attributed to each family:')

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

from joblib import Parallel, delayed                                                                     
import expectation_reflection as ER                                                                      
s_centered = s - s.mean(axis=0)                                                                                           


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

# Get W-er for all families
w_families = []
w_sym_families = []
b_families = []

from joblib import Parallel, delayed                                                                     
import expectation_reflection as ER    

# useful calculation to have mean-adjusted sequences
s_centered = s - s.mean(axis=0)      

# Expectation Reflection                                                                                 
#=========================================================================================#
# Define parallel run for ER-fit method
def predict_w(s,i0,i1i2,niter_max,l2):                                                                   
    #print('i0:',i0)                                                                                     
    i1,i2 = i1i2[i0,0],i1i2[i0,1]                                                                        
    x = np.hstack([s[:,:i1],s[:,i2:]])                                                                   
    y = s[:,i1:i2]                                                                                       
    h01,w1 = ER.fit(x,y,niter_max,l2)                                                                    
    return h01,w1                                                                                        

for msa_id in range(7):
    # we want to use all sequences
    s_train = s[family_indx[msa_id]] 
    
    # Define w matrix with variable for each possible amino acid at each sequence position               
    w_ER = np.zeros((mx.sum(),mx.sum()))                                                                     
    h0 = np.zeros(mx.sum())             
    
    
    # define files to save results to
    w_file = "%s/%s_w.npy" % (data_dir, lp_names[msa_id])   
    wsym_file = "%s/%s_w_sym.npy" % (data_dir, lp_names[msa_id])   
    b_file = "%s/%s_b.npy" % (data_dir, lp_names[msa_id])        
    
    
    # w_ER calculation can be time consuming so it may be better to load it
    create_new = True # set True if you want to compute a new w_ER regardless
    if os.path.exists(w_file) and not create_new:                                                          
        w_ER = np.load(w_file)      
        w_ER_sym = np.load(w_file)                                                                               
        b = np.load(b_file)                                                                                                               
    else:                                                                                                    
        #-------------------------------                                                                     
        # parallel                                                                                           
        start_time = timeit.default_timer()                                                                  
        res = Parallel(n_jobs = n_cpu-2)(delayed(predict_w)                                                   
                (s_train,i0,i1i2,niter_max=10,l2=100.0)                                                          
                for i0 in range(n_var))                                                                      
                                                                                                             
        run_time = timeit.default_timer() - start_time                                                       
        print('run time:',run_time)                                                                          
        #------------------------------- 
        for i0 in range(n_var):
            i1,i2 = i1i2[i0,0],i1i2[i0,1]                                                                    
            
            h01 = res[i0][0]                                                                                 
            w1 = res[i0][1]
            
            h0[i1:i2] = h01                                                                                  
            w_ER[:i1,i1:i2] = w1[:i1,:]                                                                      
            w_ER[i2:,i1:i2] = w1[i1:,:]                                                                      
            
        # make w symmetric                                                                                   
        w_ER_sym = (w_ER + w_ER.T)/2.            
        b = h0
        
        np.save(w_file, w_ER)
        np.save(wsym_file, w_ER_sym)
        np.save(b_file, b)
    

    w_families.append(w_ER)
    w_sym_families.append(w_ER_sym)
    b_families.append(b)

s_families = [] # list of full msa embeded in all family sequences spaces

