import sys,os
import random
import tools
import timeit
# import pydca-ER module
import matplotlib
#matplotlib.use('agg')
#matplotlib.rcParams['text.usetex'] = True
import matplotlib.pyplot as plt
from scipy import linalg
from sklearn.preprocessing import OneHotEncoder
import expectation_reflection as ER
from direct_info import direct_info
from direct_info import sort_di
from joblib import Parallel, delayed
import numpy as np
import pickle


#========================================================================================
data_path = '/data/cresswellclayec/DCA_ER/covid_proteins/'
root_dir = '/data/cresswellclayec/DCA_ER/covid_proteins/'
# TO RUN: 		singularity exec -B /data/cresswellclayec/DCA_ER/biowulf/,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python codon_mapping.py
       
base_pairs = ['A','T','G','C']

table = { 
'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
'NNN':'X'
} 
nsp12_range = [(13442-1,13468-1),(13468-1,16236-1)]
def translate_weird_nucleotide(nuc):
	if nuc == 'N':
		nucleotide = random.choice(base_pairs)
		print('N --> ',nucleotide)
	if nuc == 'Y':
		nucleotide = random.choice(['T','C'])
	else:
		nucleotide = nuc
	return nucleotide

def translate_sequence(seq,indices=None,subject_index=0): 

    if indices is None:
        indices = np.array([i for i in len(seq)])

    if '-' in seq:
       tab_indx = seq=='-'
       seq = seq[~tab_indx]
       indices = indices[~tab_indx]
    protein =[]
    index_mapping = {}
    if len(seq)%3 == 0: 
        amino_index = 0
        for i in range(0, len(seq), 3): 
            #print('i1 %d  i2 %d i3 %d'%(indices[i],indices[i+1],indices[i+2]))
            # add index mapping from gene to amino acid array
            #if subject_index ==indices[i] or subject_index == indices[i+1] or subject_index == indices[i+2]:
            #    index_tuple = (i,i+1,i+2)
            index_mapping[indices[i]] = amino_index
            index_mapping[indices[i+1]] = amino_index
            index_mapping[indices[i+2]] = amino_index
            codon = seq[i:i + 3] 
            protein.append( table[''.join(codon).upper()] )
            amino_index += 1
    return protein, index_mapping#, index_tuple

# convert nucleotide sequence to amino acid sequences
# default parameters are for SARS CoV-2 genome
def convert_codon(subject_index=14407, subject_encoding_region='NSP12', gene_range=nsp12_range , aligned_file = root_dir+"covid_genome_full_aligned.fasta", ref_file= root_dir+"wuhan_ref.fasta"):
	from Bio import SeqIO
	with open(aligned_file,"r") as handle:


		#subject_index = 14407
		#subject_encoding_region = 'nsp12'
		column_aa = []
		column_bp = []

		n_codon = 0
		dash_codon = 0
		for i,record in enumerate(SeqIO.parse(handle, "fasta")):
			seq_array = [char for char in ''.join(record.seq).upper()]
			seq_indices = [i for i,char in enumerate(''.join(record.seq))]
			seq_range_array = []
			seq_range_indices = []
			for start,end in gene_range:
				seq_range_indices.extend( seq_indices[start:end+1])
				seq_range_array.extend(seq_array[start:end+1])

			if i==0:
				if 0: # to compare against ncbi genes and proteins
					# compared against nsp12 data @ https://www.ncbi.nlm.nih.gov/protein/1802476815
					# WORKS
					print('record: ',record.id)
					print(seq_range_indices)
					print(seq_range_array)
					print('length: ', len(seq_range_array))
					print('length % 3 ', len(seq_range_array)%3)

				# get index mapping with reference sequence to use for the entire alignment
				protein_seq,codon_index_map,subject_codon_indices = translate_sequence(''.join(seq_range_array),seq_range_indices,subject_index)
				#print(protein_seq)
				print('\n\namino acid array len:', len(protein_seq))
				print('bp to amino acid mapping len: ',len(codon_index_map))

				print('\n#------------------------ %d Mapping -----------------------------#'%subject_index)
				print('#-----------------------  Reference Seq -----------------------------#')
				#test_index = 14407
				#print( '	14408 index in array: ',test_index)
				i1,i2,i3 = subject_codon_indices
				subject_codon_indices = (seq_range_indices[i1],seq_range_indices[i2],seq_range_indices[i3])
				bp1,bp2,bp3 = subject_codon_indices
				subject_codon = [seq_range_array[i1],seq_range_array[i2],seq_range_array[i3]]
				print( '	%d codon indices: '%subject_index,subject_codon_indices)
				print( '	%d codon nucleotieds: '%subject_index,subject_codon,' --> ',table[''.join(subject_codon)])
				subject_amino_index = codon_index_map[subject_index]
				print(	'	corresponding amino acid index and letter: %d, %s'%(subject_amino_index, protein_seq[subject_amino_index]))
				print('#--------------------------------------------------------------------#\n\n')
			# add aa corresponding to gene subject index
			subject_codon = [seq_range_array[i1],seq_range_array[i2],seq_range_array[i3]]
			n_nuc = False
			dash_nuc = False
			for ii,nucleotide in enumerate(subject_codon):
				if nucleotide not in base_pairs:
					print('index %d has abnormal nucleotide: '%i,nucleotide)
					subject_codon[ii] = translate_weird_nucleotide(nucleotide)
				if nucleotide == 'N':
					n_nuc= True 
				if nucleotide == '-':
					dash_nuc=True
			if n_nuc:
				n_codon += 1
			if dash_nuc:
				dash_codon += 1
				
			try:
				subject_codon_aa = table[''.join(subject_codon)]
			except:
				print('index %d could not convert codon: '%i,''.join(subject_codon))
				pass
			amino_index = codon_index_map[subject_index]
			column_aa.append(subject_codon_aa)	
			column_bp.append(seq_array[subject_index])	
			#print(	'	seq %d aa letter:  %s'%(i, subject_codon_aa))

	print('\n\nThere were %d codons with N nucleotides and %d codons with dashes...\n\n'%(n_codon,dash_codon))
	print('#--------------------------------------------------------------------#')
	print('\n\nSaving...')
	np.save('%d_aa_column.npy'%subject_index,column_aa)
	np.save('%d_bp_column.npy'%subject_index,column_bp)
	with open('%s_codon_index_map.pkl'%subject_encoding_region, 'wb') as f:
		pickle.dump(codon_index_map, f)
	f.close()
	print('...Done\n')


def get_aa_pair_counts(pos1,pos1_gene_range, pos2, pos2_gene_range, aligned_file = root_dir+"covid_genome_full_aligned.fasta", ref_file= root_dir+"wuhan_ref.fasta"):

	from Bio import SeqIO
	with open(aligned_file,"r") as handle:
		aa_pairs = []

		for i,record in enumerate(SeqIO.parse(handle, "fasta")):
			seq_array = [char for char in ''.join(record.seq).upper()]
			seq_indices = [i for i,char in enumerate(''.join(record.seq))]

			# get position 1 sequence array
			pos1_seq_range_array = []
			pos1_seq_range_indices = []
			for start,end in pos1_gene_range:
				pos1_seq_range_indices.extend( seq_indices[start:end+1])
				pos1_seq_range_array.extend(seq_array[start:end+1])

			# get position 2 sequence array
			pos2_seq_range_array = []
			pos2_seq_range_indices = []
			for start,end in pos2_gene_range:
				pos2_seq_range_indices.extend(seq_indices[start:end+1])
				pos2_seq_range_array.extend(seq_array[start:end+1])


			if i==0:
				# get position 1 index mapping with reference sequence to use for the entire alignment
				pos1_protein_seq,pos1_codon_index_map,pos1_codon_indices = translate_sequence(''.join(pos1_seq_range_array),pos1_seq_range_indices,pos1)
				# get position 2 index mapping with reference sequence to use for the entire alignment
				pos2_protein_seq,pos2_codon_index_map,pos2_codon_indices = translate_sequence(''.join(pos2_seq_range_array),pos2_seq_range_indices,pos2)

				print('#----- Position 1 -----#')
				print('\n\namino acid array len:', len(pos1_protein_seq))
				print('bp to amino acid mapping len: ',len(pos1_codon_index_map))
				print('#----------------------#')
				print('#----- Position 2 -----#')
				print('\n\namino acid array len:', len(pos1_protein_seq))
				print('bp to amino acid mapping len: ',len(pos1_codon_index_map))
				print('#----------------------#')

				print('\n#------------------------Mapping Pair:  %d, %d -----------------------------#'%(pos1,pos2) )
				print('#---------------------------- Pos 1 Seq -----------------------------#')
				#test_index = 14407
				#print( '	14408 index in array: ',test_index)
				pos1_i1,pos1_i2,pos1_i3 = pos1_codon_indices
				pos1_codon_indices = (pos1_seq_range_indices[pos1_i1],pos1_seq_range_indices[pos1_i2],pos1_seq_range_indices[pos1_i3])
				pos1_bp1,pos1_bp2,pos1_bp3 = pos1_codon_indices
				pos1_codon = [pos1_seq_range_array[pos1_i1],pos1_seq_range_array[pos1_i2],pos1_seq_range_array[pos1_i3]]
				print( '	%d codon indices: '%pos1,pos1_codon_indices)
				print( '	%d codon nucleotieds: '%pos1,pos1_codon,' --> ',table[''.join(pos1_codon)])
				pos1_amino_index = pos1_codon_index_map[pos1]
				print(	'	corresponding amino acid index and letter: %d, %s'%(pos1_amino_index, pos1_protein_seq[pos1_amino_index]))
				print('#--------------------------------------------------------------------#')
				print('#---------------------------- Pos 2 Seq -----------------------------#')
				#test_index = 14407
				#print( '	14408 index in array: ',test_index)
				pos2_i1,pos2_i2,pos2_i3 = pos2_codon_indices
				pos2_codon_indices = (pos2_seq_range_indices[pos2_i1],pos2_seq_range_indices[pos2_i2],pos2_seq_range_indices[pos2_i3])
				pos2_bp1,pos2_bp2,pos2_bp3 = pos2_codon_indices
				pos2_codon = [pos2_seq_range_array[pos2_i1],pos2_seq_range_array[pos2_i2],pos2_seq_range_array[pos2_i3]]
				print( '	%d codon indices: '%pos2,pos2_codon_indices)
				print( '	%d codon nucleotieds: '%pos2,pos2_codon,' --> ',table[''.join(pos2_codon)])
				pos2_amino_index = pos2_codon_index_map[pos2]
				print(	'	corresponding amino acid index and letter: %d, %s'%(pos2_amino_index, pos2_protein_seq[pos2_amino_index]))
				print('#--------------------------------------------------------------------#\n\n')

			# add aa pair corresponding to gene index pos1 and pos2
			pos1_codon = [pos1_seq_range_array[pos1_i1],pos1_seq_range_array[pos1_i2],pos1_seq_range_array[pos1_i3]]
			pos2_codon = [pos2_seq_range_array[pos2_i1],pos2_seq_range_array[pos2_i2],pos2_seq_range_array[pos2_i3]]
			# Position 1
			for ii,nucleotide in enumerate(pos1_codon):
				if nucleotide not in base_pairs:
					print('index %d has abnormal nucleotide: '%i,nucleotide)
					pos1_codon[ii] = translate_weird_nucleotide(nucleotide)
			try:
				pos1_codon_aa = table[''.join(pos1_codon)]
			except:
				print('index %d could not convert codon: '%i,''.join(pos1_codon))
				pass
			# Position 2
			for ii,nucleotide in enumerate(pos2_codon):
				if nucleotide not in base_pairs:
					print('index %d has abnormal nucleotide: '%i,nucleotide)
					pos2_codon[ii] = translate_weird_nucleotide(nucleotide)
			try:
				pos2_codon_aa = table[''.join(pos2_codon)]
			except:
				print('index %d could not convert codon: '%i,''.join(pos2_codon))
				pass

			aa_pairs.append((pos1_codon_aa,pos2_codon_aa))	

	print('#--------------------------------------------------------------------#')
	print('\n\nSaving Pairs Array...')
	np.save('%d_%d_aa_pairs.npy'%(pos1,pos2),aa_pairs)
	print('...Done\n')


# If you want to run codon_mapping as a driver (this was done for the SARS CoV-2 genome analysis
if __name__ == "__main__":
    # Swarm aligned file 
    msa_file = root_dir+"covid_genome_full_aligned.fasta"
    ref_file = root_dir+"wuhan_ref.fasta"
    
    encoding_ranges =	{
    			'NSP2' : [(806-1,2719-1)],
    			'NSP3' : [(2720-1,8554-1)],
    			'NSP12' : [(13442-1,13468-1),(13468-1,16236-1)],
    			'NSP13' : [(16237-1,18039-1)],
    			'NSP14a2' : [(18040-1,19620-1)],
    			'NSP15' : [(19621-1,20658-1)],
    			'NSP16' : [(20659-1,21552-1)],
    			'S' : [(21563-1,25384-1)],
    			'ORF3a' : [(25393-1,26220-1)],
    			'E' : [(26245-1,26472-1)],
    			'ORF7a' : [(27394-1,27759-1)],
    			'ORF7b' : [(27756-1,27887-1)],
    			'N' : [(28274-1,29533-1)],
    			'Full' : [(266-1,29674-1)]
    			}
    
    # RUN COMMAND:
    # singularity exec -B /data/cresswellclayec/DCA_ER/biowulf/,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python codon_mapping.py 14407 NSP12
    # ---------------------------------- Find Subject AA-BP ------------------------------------------ #
    
    if len(sys.argv) == 3:
    	subject_index = int(sys.argv[1])
    	subject_encoding_region = sys.argv[2] 
    	gene_range = encoding_ranges[subject_encoding_region]
    
    
    	convert_codon(subject_index=subject_index, subject_encoding_region=subject_encoding_region, gene_range = gene_range)
    
    	aa = np.load('%d_aa_column.npy'%subject_index)
    	bp = np.load('%d_bp_column.npy'%subject_index)
    
    	unique, counts = np.unique(aa, return_counts=True)
    	bp_unique, bp_counts = np.unique(bp, return_counts=True)
    
    	print(unique)	
    	print(counts)	
    	print(bp_unique)	
    	print(bp_counts)	
    # ------------------------------------------------------------------------------------------------ #
    
    # RUN COMMAND:
    # singularity exec -B /data/cresswellclayec/DCA_ER/biowulf/,/data/cresswellclayec/DCA_ER/covid_proteins /data/cresswellclayec/DCA_ER/LADER.simg python codon_mapping.py 3036 NSP3 14407 NSP12
    # ----------------------------------- Get Pair AA Counts ----------------------------------------- #
    elif len(sys.argv) == 5:
    	pos1 = int(sys.argv[1])
    	pos1_encoding_region = sys.argv[2] 
    	pos1_gene_range = encoding_ranges[pos1_encoding_region]
    
    	pos2 = int(sys.argv[3])
    	pos2_encoding_region = sys.argv[4] 
    	pos2_gene_range = encoding_ranges[pos2_encoding_region]
    
    
    
    	get_aa_pair_counts(pos1,pos1_gene_range, pos2, pos2_gene_range)
    
    # ------------------------------------------------------------------------------------------------ #
    
    else: 
    	print('Incorrect parameters')
    
