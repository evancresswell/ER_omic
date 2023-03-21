# Data Processing Added Steps:
## 2018.12.24: replace 'Z', 'X', and gap by elements in the same columns and with probability
## 2018.12.26: separate remove gaps (first) and remove conserved positions (last)



import numpy as np
import pandas as pd
import pickle
from pypdb import Query
from time import sleep
import os, sys
from pathlib import Path

from Bio import SeqIO
import Bio.PDB, warnings
from Bio.PDB import *
from scipy.spatial import distance_matrix
from Bio import BiopythonWarning
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from joblib import Parallel, delayed
pdb_parser = Bio.PDB.PDBParser()

from scipy.spatial import distance_matrix
from urllib.error import HTTPError
from prody import *
from prody import searchPfam, fetchPfamMSA
#from ProDy import *
#from ProDy.prody import searchPfam, fetchPfamMSA


# ============================================================================================================================================ #
# -------------------------------------------- Connecting PDB with MSA ----------------------------------------------------------------------- #
# ============================================================================================================================================ #
def pdb2msa(pdb_file, pdb_dir, create_new=True):
    pdb_id = os.path.basename(pdb_file)[3:7]

    # if you just want to load ProDy Dataframe (if it exists)
    if os.path.exists('%s/%s_pdb_df.csv' % (pdb_dir, pdb_id)) and not create_new:
        prody_df = pd.read_csv('%s/%s_pdb_df.csv' % (pdb_dir, pdb_id))
        return [prody_df]

    # Load pdb structure
    # print(pdb_id)
    chain_matches = {}
    for record in SeqIO.parse(pdb_file, "pdb-seqres"):
        print("Record id %s, chain %s" % (record.id, record.annotations["chain"]))
        print(record.dbxrefs)
    pdb_model = pdb_parser.get_structure(str(pdb_id), pdb_file)[0]
    

    # create DataFram of Pfam ID and metadata matching the PDB strucutre by ProDy-BLAST search 
    prody_columns =  ['PDB ID' , 'Chain', 'Polypeptide Index', 'Pfam', 'accession', 'class', 'id', 'type', 'PDB Sequence']
    float_values = ['bitscore', 'ind_evalue', 'cond_evalue']
    int_values = ['ali_end', 'ali_start', 'end', 'hmm_end', 'hmm_start', 'start']
    prody_df = None
    prody_dict = {}

    # loop through different polypeptide chains to search for matching MSAs
    for chain in pdb_model.get_chains():

        ppb = PPBuilder().build_peptides(chain)
        if len(ppb) == 0: 
            print('\nChain %s has no polypeptide sequence\n' % chain.get_id())
            continue

        for i, pp in enumerate(ppb):
            # get amino acid sequence from chain
            poly_seq = list()
            for char in str(pp.get_sequence()):
                poly_seq.append(char)
            print('\nChain %s polypeptide %d (length %d): ' % (chain.get_id(), i, len(''.join(poly_seq))),''.join(poly_seq))

            # Search for matching Pfam to given sequence with ProDy
            try:
                prody_search = searchPfam(''.join(poly_seq), timeout=300)
                print(prody_search)
            except Exception as e:
                print('Error with prody.searchPfam: ', e, '\n')
                continue

            
            # add search results to ProDy-search DataFrame 
            for pfam_key in prody_search.keys():
                ps = prody_search[pfam_key]
           
                # prody_alignment_values = list(ps['locations'].values())
                prody_alignment_values = []
                for l_key in ps['locations'].keys():
                    if l_key in float_values:
                        prody_alignment_values.append(float(ps['locations'][l_key]))
                    elif l_key in int_values:
                        prody_alignment_values.append(int(ps['locations'][l_key]))
                    else:
                        prody_alignment_values.append(ps['locations'][l_key])

                prody_lst = [pdb_id, chain.get_id(), i, pfam_key, ps['accession'], ps['class'], ps['id'], ps['type'], ''.join(poly_seq)] + prody_alignment_values

                if prody_df is None:
                    prody_alignment_columns = [key for key in ps['locations'].keys()]
                    prody_df = pd.DataFrame([prody_lst], columns = prody_columns + prody_alignment_columns)
                else:
                    prody_df.loc[len(prody_df.index)] = prody_lst


    # reorder df to get best matches on top.
    if prody_df is None:
        print('PDB2MSA ERROR: No pdb matches found using prody.searchPfam for any of the chains/polypeptide sequences. !!')
        sys.exit(42) # its the answer
    print('sorting ProdyDataframe')
    
    prody_df = prody_df.sort_values(by='bitscore', ascending = False).sort_values(by='ind_evalue') # could be cond_evalue??
    prody_df = prody_df.reset_index(drop=True)

    prody_df.to_csv('%s/%s_pdb_df.csv' % (pdb_dir, pdb_id))

    return [prody_df]




def get_tpdb(s, ali_start_indx, ali_end_indx, pfam_start_indx, pfam_end_indx, aligned_pdb_str):
    # Find the reference sequence in the Alignment
    #	-- Requires DataFrame from data_processing.pdb2msa() function as input

    alignment_len =  ali_end_indx - ali_start_indx + 1
    print('alignment length: ', alignment_len)
    print('PDB sequence: (len %d)' % len(aligned_pdb_str), aligned_pdb_str)
    #print(alignment_len, len(aligned_pdb_str), aligned_pdb_str)

    min_ham = alignment_len
    max_pair_score = 0
    min_indx = -1
    best_alignment = None
    for i, seq in enumerate(s):
        gap_seq = seq == '-'  # returns True/False for gaps/no gaps
        subject = seq[~gap_seq]
        seq_str = ''.join(subject).upper()
        aligned_seq_str = seq_str[pfam_start_indx : pfam_end_indx+1]
        alignments = pairwise2.align.globalxx(aligned_pdb_str, aligned_seq_str)
        
        if len(alignments) == 0:
            continue
        try:
            pair_score = alignments[0].score
        except(AttributeError):
            pair_score = alignments[0][2]

        if pair_score > max_pair_score:
            #min_ham = ham_dist
            max_pair_score = pair_score
            min_indx = i
            print(format_alignment(*alignments[0]))
            best_alignment = alignments[0]
            print('match upgrade at' , i)
            #print('%d: hamm dist=%d, pairwise score=%f\n' % (i, ham_dist, pair_score))
            print('%d: pairwise score=%f\n' % (i, pair_score))
            print('lengths: ', len(aligned_pdb_str), len(aligned_seq_str))
            print(best_alignment)

    gap_seq = s[min_indx] == '-'
    print('best match is sequence %d with hamming distance of %d (length %d)' % (min_indx, min_ham, len(s[min_indx][~gap_seq])))
    
    if best_alignment is None:
        print('ERROR could not find alignment with score better than 0..') 
        
    # get matching seq for both sequences (no gaps in either)
    try:
        aligned_pdb_char_array = np.array([char for char in best_alignment.seqA])
        aligned_ref_char_array = np.array([char for char in best_alignment.seqB])
    except(AttributeError):
        aligned_pdb_char_array = np.array([char for char in best_alignment[0]])
        aligned_ref_char_array = np.array([char for char in best_alignment[1]])
       
    
    # get array of gaps for both sequences
    seqA_gaps = aligned_pdb_char_array == '-'
    seqB_gaps = aligned_ref_char_array == '-'
    aligned_gaps = np.logical_or(seqA_gaps, seqB_gaps)

    
    # create index array for reference sequence so we know which msa columns associated with aligned arrays
    pdb_count = 0
    ref_count = 0
    gap_ref_index = -1 * np.ones(len(aligned_ref_char_array), dtype=int)
    gap_pdb_index = -1 * np.ones(len(aligned_pdb_char_array), dtype=int)
    for i, char in enumerate(aligned_ref_char_array):
        if char !='-':
            gap_ref_index[i] = int(ref_count)
            ref_count += 1
        if aligned_pdb_char_array[i] !='-':
            gap_pdb_index[i] = int(pdb_count)
            pdb_count += 1            
    
    # get columns to remove (gap in PDB) in MSA
    pdb_gap_cols_in_ref = gap_ref_index[seqA_gaps]
    print(len(pdb_gap_cols_in_ref), pdb_gap_cols_in_ref)

    # get s_index for mapping msa to pdb sequence.
    pdb_s_index = gap_pdb_index[~aligned_gaps]
    print('PDB index map: ', len(pdb_s_index), pdb_s_index)
    
    # Extract further infor for aligned seqs.
    aligned_pdb_nogap = aligned_pdb_char_array[~aligned_gaps]
    aligned_ref_nogap = aligned_ref_char_array[~aligned_gaps]
    print('\n aligned PDB and ref seq:')
    print(len(aligned_pdb_nogap),aligned_pdb_nogap)
    print(len(aligned_ref_nogap),aligned_ref_nogap)
    
    # Trim By gaps in Ref seq (tbr). Then Trim By gaps in Pdb seq (tpb)
    s_tbr = s[:, ~gap_seq]
    s_tbr = s_tbr[:,pfam_start_indx : pfam_end_indx+1]
    print('s_tbr[tpdb] = ', s_tbr[min_indx])
    print(s_tbr.shape)
    s_tbr_tbp = np.delete(s_tbr, pdb_gap_cols_in_ref, axis=1)
    print(s_tbr_tbp.shape)

    # printed ref seq should be the same as the fully alinged, gapless pdb and ref seqs above.
    print('MSA trimmed by reference and pdb_sequence: ', len(s_tbr_tbp[min_indx]), s_tbr_tbp[min_indx])

    return min_indx, best_alignment, s_tbr_tbp, pdb_s_index


def query_pdb(msa, i):
    # This function allows for parallel ProDy-BLAST query
    # 	-- ISSUE: Biowulf does not like when you make a lot of query/downloads to outside servers
    seq = msa[i]
    gap_seq = seq == '-'  # returns True/False for gaps/no gaps
    subject = seq[~gap_seq]
    seq_str = ''.join(subject)

    try: 
        q = Query(seq_str.upper(),
          query_type="sequence",
          return_type="polymer_entity")
        search_results = q.search()
    
    
        # print(search_results)
        # print('sequence %d: ' % i, search_results['result_set'][0]['identifier'])
        
        # parse BLAST PDB search results
        pdb_id = search_results['result_set'][0]['identifier']
        score = search_results['result_set'][0]['score']
        
        match_dict = search_results['result_set'][0]['services'][0]['nodes'][0]['match_context'][0]
        identity = match_dict['sequence_identity']
        e_value = match_dict['evalue']
        bitscore = match_dict['bitscore']
        alignment_length = match_dict['alignment_length']
        mismatches = match_dict['mismatches']
        gaps_open = match_dict['gaps_opened']
        query_beg = match_dict['query_beg']
        query_end = match_dict['query_end']
        subject_beg = match_dict['subject_beg']
        subject_end = match_dict['subject_end']
        query_len = match_dict['query_length']
        subject_len = match_dict['subject_length']
        query_aligned_seq = match_dict['query_aligned_seq']
        subject_aligned_seq = match_dict['subject_aligned_seq']
        pass
    except(UserWarning):
        print('bad query')
        # fill with empty search result 
        pdb_id = '0000' 
        score = 0.0 
        
        match_dict = {} 
        identity = 0. 
        e_value = 1. 
        bitscore = 0. 
        alignment_length = 0 
        mismatches = 0 
        gaps_open = 0  
        query_beg = 0  
        query_end = 0 
        subject_beg = 0 
        subject_end = 0 
        query_len = 0  
        subject_len = 0 
        query_aligned_seq = 0
        subject_aligned_seq = 0
        pass
 
    # sleep(10) # Don't query too fast or biowulf complains
    return  [i, pdb_id, score, identity, e_value, bitscore, alignment_length, mismatches, gaps_open, \
             query_beg, query_end, subject_beg, subject_end, query_len, subject_len, \
             query_aligned_seq, subject_aligned_seq]

 
def find_best_pdb(pfam_id, data_path, pdb_dir, n_cpu=20):
    # Given a Pfam_id find the best corresponding PDB strucutre 
    # -- This search direction is inherently flawed because of how much time it takes and some MSA's just dont have matching PDB strucutres
    from IPython.display import HTML
    # Set Columns for eventual DataFrame sort.
    columns = ['MSA Index', 'PDB ID', 'Score', 'Identity', 'E-value', 'Bitscore', 'Alignment Length', 'Mismatches', 'Gaps Opened', \
                                  'Query Beg', 'Query End', 'Subject Beg', 'Subject End', 'Query Len', 'Subject Len', \
                                  'Query Aligned Seq', 'Subject Aligned Seq']

    # Import from installed package
    print('looking for pdb_references in %s ' % pdb_dir)
    if os.path.exists('%s/%s_pdb_references.pkl' % (pdb_dir, pfam_id)):
        # if the query has already been done, load the raw query dataframe 
        try:
            pdb_refs_ecc  = pd.read_pickle('%s/%s_pdb_references.pkl' % (pdb_dir, pfam_id))
        except(ValueError): # Most likely unsupported pickle protocol from trying to load up-to-date pickled df with old python version... read csv.
            pdb_refs_ecc  = pd.read_csv('%s/%s_pdb_references.csv' % (pdb_dir, pfam_id))
    else:
        msa = load_msa(data_path, pfam_id)
        dup_rows = []

        print('Finding best PDB match for (Searching %d sequences)' % len(msa), pfam_id)
        #try:  # try parallel run, if it doesent work just do in serial
        if 0:
            print('\n\nQuerying PDB database in Parallel...\n')
            from joblib import Parallel, delayed
            pdb_matches = Parallel(n_jobs=2, verbose=10)(delayed(query_pdb)(msa, i) for i in range(len(msa)))     
            pdb_refs_ecc = pd.DataFrame(pdb_matches, columns=columns)
        #except(AttributeError):
        else:
            # Parallel query didnt work so do it in serial
            print('Parallel query didnt work so we do it in serial')
            pdb_matches = {}
            for i, seq in enumerate(msa):
                print('query sequence %d of %d ' % (i, msa.shape[0]))
                # --- check for duplicate row --- #
                if [a for a in seq] in dup_rows:
                    continue
                else:
                    dup_rows.append([a for a in seq])
                # ------------------------------- #
                try:
                    gap_seq = seq == '-'  # returns True/False for gaps/no gaps
                    subject = seq[~gap_seq]
                    seq_str = ''.join(subject)
                    q = Query(seq_str.upper(),
                      query_type="sequence",
                      return_type="polymer_entity")
                    search_results = q.search()
                
                
                    # print(search_results)
                    # print('sequence %d: ' % i, search_results['result_set'][0]['identifier'])
                    
                    # parse BLAST PDB search results
                    pdb_id = search_results['result_set'][0]['identifier']
                    score = search_results['result_set'][0]['score']
                    
                    match_dict = search_results['result_set'][0]['services'][0]['nodes'][0]['match_context'][0]
                    identity = match_dict['sequence_identity']
                    e_value = match_dict['evalue']
                    bitscore = match_dict['bitscore']
                    alignment_length = match_dict['alignment_length']
                    mismatches = match_dict['mismatches']
                    gaps_open = match_dict['gaps_opened']
                    query_beg = match_dict['query_beg']
                    query_end = match_dict['query_end']
                    subject_beg = match_dict['subject_beg']
                    subject_end = match_dict['subject_end']
                    query_len = match_dict['query_length']
                    subject_len = match_dict['subject_length']
                    query_aligned_seq = match_dict['query_aligned_seq']
                    subject_aligned_seq = match_dict['subject_aligned_seq']
                    
                    pdb_matches[i] = [i, pdb_id, score, identity, e_value, bitscore, alignment_length, mismatches, gaps_open, \
                                      query_beg, query_end, subject_beg, subject_end, query_len, subject_len, \
                                      query_aligned_seq, subject_aligned_seq]
                                
                    #print('query %d: good' % i)
                except(UserWarning):
                    #print('query %d: bad' % i)
                    pass
            pdb_refs_ecc = pd.DataFrame.from_dict(pdb_matches, orient='index', columns=columns)
            if len(pdb_refs_ecc) == 0:
                print('PDB-sequence query yields no results!!')
                sys.exit(23)
     
            pdb_refs_ecc.to_pickle('%s/%s_pdb_references.pkl' % (pdb_dir, pfam_id))
            pdb_refs_ecc.to_csv('%s/%s_pdb_references.csv' % (pdb_dir, pfam_id))

        print('\nPDB query Dictionary length: ', len(pdb_matches)) 
    print('Raw-Query PDB dataframe gives %d matches... \n' % len(pdb_refs_ecc))
    pdb_sorted = pdb_refs_ecc.loc[pdb_refs_ecc['Identity']>.90] 
    print('PDB dataframe with better than .90 Identity: %d' % len(pdb_sorted))
    pdb_sorted = pdb_sorted.loc[pdb_sorted['Gaps Opened']==0]
    print('PDB dataframe with no gaps open: %d' % len(pdb_sorted))
    pdb_sorted = pdb_sorted.sort_values(by=['Score', 'Bitscore'], ascending = False).sort_values(by='E-value')
    pdb_matched = pdb_sorted.loc[pdb_sorted['Alignment Length']==pdb_sorted['Query Len']]
    pdb_matched_sorted = pdb_matched.sort_values(by=['Alignment Length', 'Bitscore', 'Score'], ascending = False).sort_values(by='E-value')
    pdb_matched_matched_sorted = pdb_matched_sorted.reset_index(drop=True)
    print('Sorted PDB matches (%d matches): \n' % len(pdb_matched_sorted), pdb_matched_sorted.head())
    return pdb_matched_sorted


# ============================================================================================================================================ #
# -------------------------------------------------------------------------------------------------------------------------------------------- #
# ============================================================================================================================================ #


# ============================================================================================================================================ #
# ------------------------------------------ MSA Processing ---------------------------------------------------------------------------------- #
# ============================================================================================================================================ #

# ------------------------------ #
def remove_bad_seqs(s, tpdb, fgs=0.3, trimmed_by_refseq=True,return_seq_indx=False):
    # remove sequence rows with too many gaps 
    #     -- update tpdb

    # if trimmed by reference sequence, create a temp matrix to find bad sequences
    if not trimmed_by_refseq:
        s_temp = s.copy()
        gap_pdb = s[tpdb] == '-'  # returns True/False for gaps/no gaps in reference sequence
        s_temp = s_temp[:, ~gap_pdb]  # removes gaps in reference sequence
    else:
        s_temp = s

    # remove bad sequences having a gap fraction of fgs  
    l, n = s_temp.shape

    frequency = [(s_temp[t, :] == '-').sum() / float(n) for t in range(l)]
    bad_seq = [t for t in range(l) if frequency[t] > fgs]
    print(len(bad_seq))
    new_s = np.delete(s, bad_seq, axis=0)
    # Find new sequence index of Reference sequence tpdb
    seq_index = np.arange(s.shape[0])
    seq_index = np.delete(seq_index, bad_seq)
    new_tpdb = np.where(seq_index == tpdb)
    print("After removing bad sequences, tpdb is now ", new_tpdb[0][0])
    if not return_seq_indx:
        return new_s, new_tpdb[0][0]
    else:
        return new_s, new_tpdb[0][0], bad_seq

# ------------------------------ #
def remove_seqs_list(s, tpdb, seqs_to_remove):
    # remove sequence rows in list seqs_to_remove
    #     -- update tpdb
    l, n = s.shape

    new_s = np.delete(s, seqs_to_remove, axis=0)
    # Find new sequence index of Reference sequence tpdb
    seq_index = np.arange(s.shape[0])
    seq_index = np.delete(seq_index, seqs_to_remove)
    new_tpdb = np.where(seq_index == tpdb)
    print("After removing bad sequences, tpdb is now ", new_tpdb[0][0])

    return new_s, new_tpdb[0][0]



# ------------------------------ #
def remove_bad_cols(s, fg=0.3, fc=0.9):
    # remove positions having a fraction fc of converved residues or a fraction fg of gaps 

    l, n = s.shape
    # gap positions:
    frequency = [(s[:, i] == '-').sum() / float(l) for i in range(n)]
    cols_gap = [i for i in range(n) if frequency[i] > fg]

    # conserved positions:
    frequency = [max(np.unique(s[:, i], return_counts=True)[1]) for i in range(n)]
    cols_conserved = [i for i in range(n) if frequency[i] / float(l) > fc]

    cols_remove = cols_gap + cols_conserved

    return np.delete(s, cols_remove, axis=1), cols_remove


# ------------------------------ #
def find_bad_cols(s, fg=0.2):
    # remove positions having a fraction fg of gaps
    l, n = s.shape
    # gap positions:
    frequency = [(s[:, i] == '-').sum() / float(l) for i in range(n)]
    bad_cols = [i for i in range(n) if frequency[i] > fg]

    # return np.delete(s,gap_cols,axis=1),np.array(gap_cols)
    return np.array(bad_cols)


# ------------------------------ #
def find_conserved_cols(s, fc=0.8):
    # remove positions having a fraction fc of converved residues
    l, n = s.shape

    # conserved positions:
    frequency = [max(np.unique(s[:, i], return_counts=True)[1]) for i in range(n)]
    conserved_cols = [i for i in range(n) if frequency[i] / float(l) > fc]

    # return np.delete(s,conserved_cols,axis=1),np.array(conserved_cols)
    return np.array(conserved_cols)


# ------------------------------ #
def number_residues(s):
    # number of residues at each position
    l, n = s.shape
    mi = np.zeros(n)
    for i in range(n):
        s_unique = np.unique(s[:, i])
        mi[i] = len(s_unique)

    return mi


# ------------------------------ #
def convert_letter2number(s):
    letter2number = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, \
                     'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19, '-': 20,
                     'U': 21}
    # ,'B':20, 'Z':21, 'X':22}
    try:
        l, n = s.shape
    except(ValueError):
        n = s.shape[0] # if s is only one row.
        return np.array([letter2number[s[i].upper()]  for i in range(n)])
    # making sure all amino acids are uppercase # # this is done in PYDCA as well. though there are references that say lowercase means no-good
    return np.array([letter2number[s[t, i].upper()] for t in range(l) for i in range(n)]).reshape(l, n)


# ------------------------------ #
def convert_number2letter(s):
    number2letter = {0: 'A', 1: 'C', 2: 'D', 3: 'E', 4: 'F', 5: 'G', 6: 'H', 7: 'I', 8: 'K', 9: 'L', \
                     10: 'M', 11: 'N', 12: 'P', 13: 'Q', 14: 'R', 15: 'S', 16: 'T', 17: 'V', 18: 'W', 19: 'Y', 20: '-',
                     21: 'U'}
    try:
        l, n = s.shape
        return np.array([number2letter[s[t, i]] for t in range(l) for i in range(n)]).reshape(l, n)
    except(ValueError):
        print(s)
        return np.array([number2letter[r] for r in s])


# ------------------------------ #
# 2018.12.24: replace value at a column with probility of elements in that column
def value_with_prob(a, p1):
    """ generate a value (in a) with probability
    input: a = np.array(['A','B','C','D']) and p = np.array([0.4,0.5,0.05,0.05]) 
    output: B or A (likely), C or D (unlikely)
    """
    p = p1.copy()
    # if no-specific prob --> set as uniform distribution
    if p.sum() == 0:
        p[:] = 1. / a.shape[0]  # uniform
    else:
        p[:] /= p.sum()  # normalize

    ia = int((p.cumsum() < np.random.rand()).sum())  # cordinate

    return a[ia]



# ------------------------------ #
def find_and_replace(s, z, a):
    """ find positions of s having z and replace by a with a probality of elements in s column
    input: s = np.array([['A','Q','A'],['A','E','C'],['Z','Q','A'],['A','Z','-']])
           z = 'Z' , a = np.array(['Q','E'])    
    output: s = np.array([['A','Q','A'],['A','E','C'],['E','Q','A'],['A','Q','-']]           
    """
    print('Find and Replace in sequence.\nNOTE: This can be done in Parallel\nLook at data_processing.find_and_replace_parallel()')
    xy = np.argwhere(s == z)

    for it in range(xy.shape[0]):
        t, i = xy[it, 0], xy[it, 1]

        na = a.shape[0]
        p = np.zeros(na)
        for ii in range(na):
            p[ii] = (s[:, i] == a[ii]).sum()

        s[t, i] = value_with_prob(a, p)
    return s

# ------------------------------ #
def find_and_replace_section(s_sec,z, a,i0):
    # parallel implementation
    """ find positions of s having z and replace by a with a probality of elements in s column
    input: s = np.array([['A','Q','A'],['A','E','C'],['Z','Q','A'],['A','Z','-']])
           z = 'Z' , a = np.array(['Q','E'])    
    output: s = np.array([['A','Q','A'],['A','E','C'],['E','Q','A'],['A','Q','-']]           
    """ 
    s = np.copy(s_sec[i0])
    xy = np.argwhere(s == z)

    try:
        for it in range(xy.shape[0]):
            t,i = xy[it,0],xy[it,1]
    
            na = a.shape[0]
            p = np.zeros(na)    
            for ii in range(na):
                p[ii] = (s[:,i] == a[ii]).sum()
    
            s[t,i] = value_with_prob(a, p)
    except(ValueError):
        print('error in section %d' % i0)
        print(s.flags)
        sys.exit(0)

    return s 


# ------------------------------ #
def find_and_replace_parallel(s,z,a, cpus_per_job):
    from math import ceil
    ax = 0
    ncpu = max(cpus_per_job - 4, 2)
    n_sec = ncpu
    n_seq = ceil(float(s.shape[0]) / ncpu)
    print('dealing with %f sequences in each of the %d sections' % (n_seq,n_sec))
    timeout=9999999
    # create sectioned version of s
    print('in find_and_replace_parallel: ', s.shape)
    print(n_sec)
    s_sec = np.array_split(s, n_sec, axis=ax)
    print(s_sec[0].shape)
    for ii in range(len(s_sec)):
        s_sec[ii].setflags(write=1)

    res = Parallel(n_jobs = ncpu, timeout=timeout)(delayed(find_and_replace_section)\
        (s_sec,z, a,i0)\
        for i0 in range(int(n_sec)))
    res = np.concatenate(res,axis=ax)


    print('parallel result: ', res.shape)
    return res
    

# ------------------------------ #
def replace_lower_by_higher_prob(s, p0=0.3):
    # input: s: 1D numpy array ; threshold p0
    # output: s in which element having p < p0 were placed by elements with p > p0, according to prob

    # f = itemfreq(s)  replaced by next line due to warning
    f = np.unique(s, return_counts=True)
    # element and number of occurence
    a, p = f[0], f[1].astype(float)

    # probabilities    
    p /= float(p.sum())

    # find elements having p > p0:
    iapmax = np.argwhere(p > p0).reshape((-1,))  # position

    apmax = a[iapmax].reshape((-1,))  # name of aminoacid
    pmax = p[iapmax].reshape((-1,))  # probability

    # find elements having p < p0
    apmin = a[np.argwhere(p < p0)].reshape((-1,))

    if apmin.shape[0] > 0:
        for a in apmin:
            ia = np.argwhere(s == a).reshape((-1,))
            for iia in ia:
                s[iia] = value_with_prob(apmax, pmax)

    return s


# ------------------------------ #
def write_FASTA(msa, pfam_id, s_ipdb, number_form=True, processed=True, path='./'):
    # Write MSA to fasta format file
    # Processed MSA to file in FASTA format
    msa_outfile = 'MSA_' + pfam_id + '.fa'

    # Reference sequence to file in FASTA format
    ref_outfile = 'ref_' + pfam_id + '.fa'
    ref_seq = s_ipdb

    # print("Reference Sequence (shape=",msa[ref_seq].shape,"):\n",msa[ref_seq])

    if number_form:
        msa_letters = convert_number2letter(msa)
        ref_array = msa_letters[ref_seq]
        if not processed:
            gap_ref = ref_array == '-'  # remove gaps from reference array
            ref_letters = msa_letters[ref_seq][~gap_ref]
        else:
            ref_letters = msa_letters[ref_seq]
    else:
        msa_letters = msa
        ref_array = msa[ref_seq]
        gap_ref = ref_array == '-'  # remove gaps from reference array
        ref_letters = ref_array[~gap_ref]

    printing = False
    if printing:
        print("Reference Sequence number: ", ref_seq)
        print("Reference Sequence (shape=", ref_letters.shape, "):\n", ref_letters)

        print("Writing processed MSA (shape=", msa_letters.shape, ") to FASTA files:\n", msa_letters)

    # First save reference sequence to FASTA file
    ref_str = ''
    ref_list = ref_letters.tolist()
    ref_str = ref_str.join(ref_list)
    if processed:
        with open(path + ref_outfile, 'w') as fh:
            fh.write('>{}\n{}\n'.format(pfam_id + ' | REFERENCE', ref_str))
    else:
        ref_outfile = 'orig_ref_' + pfam_id + '.fa'
        with open(path + ref_outfile, 'w') as fh:
            fh.write('>{}\n{}\n'.format(pfam_id + ' | REFERENCE', ref_str))

    # Next save MSA to FAST file

    with open(path + msa_outfile, 'w') as fh:
        for seq_num, seq in enumerate(msa_letters):
            msa_list = seq.tolist()
            msa_str = ''
            msa_str = msa_str.join(msa_list)
            if seq_num == ref_seq:
                fasta_header = pfam_id + ' | REFERENCE'
            else:
                fasta_header = pfam_id
            if printing:
                print(msa_str)
            fh.write('>{}\n{}\n'.format(fasta_header, msa_str))

# =========================================================================================
# noinspection PyBroadException
def load_msa(data_path, pfam_id):
    s = np.load('%s/%s/msa.npy' % (data_path, pfam_id)).T

    # convert bytes to str
    try:
        s = np.array([s[t, i].decode('UTF-8') for t in range(s.shape[0]) \
                      for i in range(s.shape[1])]).reshape(s.shape[0], s.shape[1])
    # print("shape of s (after UTF-8 decode):\n",s.shape)
    except:
        print("\n\nUTF not decoded, pfam_id: %s\n" % pfam_id, s.shape)
        print(s[0])
        # last ditch effort -- takes much more memory/time to decode
        try:
            s = np.array([str(s[t, i])[2] for t in range(s.shape[0]) \
                      for i in range(s.shape[1])]).reshape(s.shape[0], s.shape[1])
            return s
        except:
            print('trying long way after UTF exception didnt work: ', sys.exc_info()[0])
            pass
        print("Exception: ", sys.exc_info()[0])
        # Create list file for missing pdb structures
        if not os.path.exists('missing_MSA.txt'):
            file_missing_msa = open("missing_MSA.txt", 'w')
            file_missing_msa.write("%s\n" % pfam_id)
            file_missing_msa.close()
        else:
            file_missing_msa = open("missing_MSA.txt", 'a')
            file_missing_msa.write("%s\n" % pfam_id)
            file_missing_msa.close()
        return
    return s


def data_processing_pdb2msa(data_path, pdb_df,gap_seqs=0.2, gap_cols=0.2, prob_low=0.004, 
                        conserved_cols=0.8, printing=True, out_dir='./', pdb_dir='./', letter_format=False, 
                        remove_cols=True, create_new=True):
    pfam_id = pdb_df['Pfam']
    pdb_seq = pdb_df['PDB Sequence']
    pdb_id = pdb_df['PDB ID']
    ali_start_indx = int(pdb_df['ali_start'])-1
    ali_end_indx = int(pdb_df['ali_end'])-1
    pfam_start_indx = int(pdb_df['hmm_start'])-1
    pfam_end_indx = int(pdb_df['hmm_end'])-1

    aligned_pdb_str  = pdb_df['PDB Sequence'][ali_start_indx:ali_end_indx+1]


    print('PDB ID: %s, Pfam ID: %s' % (pdb_id, pfam_id))

    np.random.seed(123456789)
    #if not create_new and os.path.exists("%s/%s_pdb_df.csv" % (out_dir, pfam_id)):
    if 0:
        print('Because create_new is False and files exist we will load preprocessed data:')
        if remove_cols:
            s = np.load("%s/%s_%s_preproc_msa.npy" % (out_dir, pfam_id, pdb_id))
            s_index = np.load("%s/%s_%s_preproc_sindex.npy" % (out_dir, pfam_id, pdb_id))
            removed_cols = np.load("%s/%s_%s_removed_cols.npy" % (out_dir, pfam_id, pdb_id))
            ref_seq = np.load("%s/%s_%s_preproc_refseq.npy" % (out_dir, pfam_id, pdb_id))
        else: 
            s = np.load("%s/%s_%s_allCols_msa.npy" % (out_dir, pfam_id, pdb_id))
            s_index = np.load("%s/%s_%s_allCols_sindex.npy" % (out_dir, pfam_id, pdb_id)) 
            removed_cols = np.load("%s/%s_%s_removed_cols.npy" % (out_dir, pfam_id, pdb_id))
            ref_seq = np.load("%s/%s_%s_allCols_refseq.npy" % (out_dir, pfam_id, pdb_id))

        if not letter_format and isinstance(s[0][0], str):
            s = convert_letter2number(s)

      
    
    # Load MSA
    s = load_msa(data_path, pfam_id)
    print('s: ', s)
    orig_seq_len = s.shape[1]
    print('Original Sequence length: ', orig_seq_len)

    
    # Using given MSA find best matching PDB structure from all available MSA sequences.
    
    if printing:
        print("\n\n#--------------------- Find PDB Sequence in MSA ---------------#")
       

       
    # Find PDB seq in MSA current
    tpdb, alignment, s, pdb_s_index = get_tpdb(s, ali_start_indx, ali_end_indx, pfam_start_indx, pfam_end_indx, aligned_pdb_str) # requires prody.searchPfam DF from pdb2msa as input

    if printing:
        print('MSA index %d matches PDB' % tpdb)
        print("#--------------------------------------------------------------#\n\n")

    if printing:
        print("#\n\n-------------------------Remove Gaps--------------------------#")
        print('Shape of s is : ', s.shape)
        print("s = \n", s)
    


 
    # remove gaps to allow alignment with PDB sequence..

    # remove columns not in alignment range
    in_range_indices = np.arange(ali_start_indx, ali_end_indx) 
    print('in range indices: ', in_range_indices)


    s_index = np.arange(s.shape[1])

    if printing:
        print("s[tpdb] shape is ", s[tpdb].shape)
        print("s = \n", s)
        print("though s still has gaps, s[%d] does not:\n" % (tpdb), s[tpdb])
        print("s shape is ", s.shape)
        print("Saving indices of reference sequence s[%d](length=%d):\n" % (tpdb, s_index.shape[0]), s_index)
        print("#--------------------------------------------------------------#\n\n")

    lower_cols = np.array([i for i in range(s.shape[1]) if s[tpdb, i].islower()])
    if printing:
        print("removing non aligned (lower case) columns in subject sequence:\n ", lower_cols, '\n')
        
        
    # lower case removal reference: https://onlinelibrary.wiley.com/doi/full/10.1002/1097-0134%2820001101%2941%3A2%3C224%3A%3AAID-PROT70%3E3.0.CO%3B2-Z
    # upper = np.array([x.isupper() for x in s[tpdb]])
    # print('select only column presenting as uppercase at the first row')
    # upper = np.array([x.isupper() for x in s[0]])
    # s = s[:,upper]


    # --- remove duplicates before processing (as done in pydca) --- #
    if 0:
        dup_rows = []
        s_no_dup = []
        for i, row in enumerate(s):
            if [a for a in row] in s_no_dup:
                if i != tpdb:   # do not want to remove reference sequence
                    dup_rows.append(i)
                else:           # we need to add the reference sequence back in even if its a duplicate row.    
                    s_no_dup.append([a for a in row])
            else:
                s_no_dup.append([a for a in row])
        if printing:
            print('found %d duplicates! (Removing...)' % len(dup_rows))
        s, tpdb = remove_seqs_list(s, tpdb, dup_rows)
    # -------------------------------------------------------------- #

    
    # - Removing bad sequences (>gap_seqs gaps) -------------------- #    
    s, tpdb, bad_seq_indx = remove_bad_seqs(s, tpdb, gap_seqs, return_seq_indx=True)  # removes all sequences (rows) with >gap_seqs gap %
    
    if printing:
        print('\nAfter removing bad sequences...\ntpdb (s_ipdb) is : ', tpdb)
        print(s.shape)
    # -------------------------------------------------------------- #

    # - Finding bad columns (>gap_cols) -------------------- #
    bad_cols = find_bad_cols(s, gap_cols)
    if printing:
        print('found bad columns :=', bad_cols)
    # ------------------------------------------------------ #

    # ------------ Replace Bad Amino Acid Letters if valid ones, two strategies --------- #
    if 1: # replace aa with potential correct aa
        # 2018.12.24:
        # replace 'Z' by 'Q' or 'E' with prob
        # print('replace Z by Q or E')
        s = find_and_replace(s, 'Z', np.array(['Q', 'E']))

        # replace 'B' by Asparagine (N) or Aspartic (D)
        # print('replace B by N or D')
        s = find_and_replace(s, 'B', np.array(['N', 'D']))

        # replace 'X' as amino acids with prob
        # print('replace X by other aminoacids')
        amino_acids = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', \
                                'T', 'V', 'W', 'Y', 'U'])
        s = find_and_replace(s, 'X', amino_acids)
    else: # replace aa with gap '-' (like in pydca)
        if printing:
            print('\n\nUsing PYDCA\'s method of replacing Z, B, and X with - instead of likely alternates...\n\n')
        s = find_and_replace(s, 'Z', np.array(['-']))
        s = find_and_replace(s, 'B', np.array(['-']))
        s = find_and_replace(s, 'X', np.array(['-']))
    # ----------------------------------------------------------------------------------- #

    
    # --------------------- Find Conserved Columns ------------------------------------ #
    conserved_cols = find_conserved_cols(s, conserved_cols)
    if printing:
        print("found conserved columns (80% repetition):\n", conserved_cols)
    # --------------------------------------------------------------------------------- #

    
    # ----------------------- Remove all Bad Columns ---------------------------------- #

    removed_cols = np.array(list(set(bad_cols) | set(conserved_cols)))

    removed_cols = np.array(list(set(removed_cols) | set(lower_cols)))
        
    # remove columns in which references sequence is different from matched PDB sequence
    #mismatch_cols = [i for i,a in enumerate(ref_seq) if a != pdb_seq[i]]
    #removed_cols = np.array(list(set(removed_cols) | set(mismatch_cols)))

    if printing:
        print("We remove conserved and bad columns with, at the following indices (len %d):\n" % len(removed_cols), removed_cols)
    # ----------------------- Remove all Bad Columns ---------------------------------- #


    # info still pased through removed_cols but this way we can interact with full msa if remove_cols is False
    if remove_cols:
        s = np.delete(s, removed_cols, axis=1)
        s_index = np.delete(s_index, removed_cols)
        pdb_s_index = np.delete(pdb_s_index, removed_cols)

    if printing and remove_cols:
        print("Removed Columns...")
        print("s now has shape: ", s.shape)
        print("s_index (length=%d) = \n" % s_index.shape[0], s_index)
        print("pdb_s_index (length=%d) = \n" % pdb_s_index.shape[0], pdb_s_index)
        print("In data processing Ref Seq (shape=", s[tpdb].shape, "): \n", s[tpdb])

    # Convert S to number format (number representation of amino acids)
    if not letter_format:
        # convert letter to number:
        s = convert_letter2number(s)
    print(s)

    # replace lower probs by higher probs 
    for i in range(s.shape[1]):
        s[:, i] = replace_lower_by_higher_prob(s[:, i], prob_low)

    np.save("%s/%s_%s_removed_cols.npy" % (out_dir, pfam_id, pdb_id), removed_cols)
    if remove_cols:
        np.save("%s/%s_%s_preproc_msa.npy" % (out_dir, pfam_id, pdb_id), s)
        np.save("%s/%s_%s_preproc_sindex.npy" % (out_dir, pfam_id, pdb_id), s_index)
        np.save("%s/%s_%s_preproc_pdb_sindex.npy" % (out_dir, pfam_id, pdb_id), pdb_s_index)
        np.save("%s/%s_%s_preproc_refseq.npy" % (out_dir, pfam_id, pdb_id), s[tpdb])
    else:
        np.save("%s/%s_%s_allCols_msa.npy" % (out_dir, pfam_id, pdb_id), s)
        np.save("%s/%s_%s_allCols_sindex.npy" % (out_dir, pfam_id, pdb_id), s_index)
        np.save("%s/%s_%s_allCols_refseq.npy" % (out_dir, pfam_id, pdb_id), s[tpdb])



    # - Removing bad sequences (>gap_seqs gaps) -------------------- #
    if printing:
        print(s.shape)
        print("In Data Processing Final Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])

    return [s, removed_cols, s_index, tpdb, pdb_s_index, bad_seq_indx]


# =========================================================================================



def data_processing_msa2pdb(data_path, pfam_id, index_pdb=0, gap_seqs=0.2, gap_cols=0.2, prob_low=0.004, 
                        conserved_cols=0.8, printing=True, out_dir='./', pdb_dir='./', letter_format=False, 
                        remove_cols=True, create_new=True, n_cpu=4):
    n_cpu = min(10, n_cpu-2) 
    np.random.seed(123456789)
    if not create_new and os.path.exists("%s/%s_pdb_query.npy" % (out_dir, pfam_id)):
        print('Because create_new is False and files exist we will load preprocessed data:')
        if remove_cols:
            s = np.load("%s/%s_preproc_msa.npy" % (out_dir, pfam_id))
            s_index = np.load("%s/%s_preproc_sindex.npy" % (out_dir, pfam_id))
            removed_cols = np.load("%s/%s_removed_cols.npy" % (out_dir, pfam_id))
            ref_seq = np.load("%s/%s_preproc_refseq.npy" % (out_dir, pfam_id))
            pdb_seq = np.load("%s/%s_preproc_pdb_seq.npy" % (out_dir, pfam_id))
        else: 
            s = np.load("%s/%s_allCols_msa.npy" % (out_dir, pfam_id))
            s_index = np.load("%s/%s_allCols_sindex.npy" % (out_dir, pfam_id))
            removed_cols = np.load("%s/%s_removed_cols.npy" % (out_dir, pfam_id))
            ref_seq = np.load("%s/%s_allCols_refseq.npy" % (out_dir, pfam_id))
            pdb_seq = np.load("%s/%s_allCols_pdb_seq.npy" % (out_dir, pfam_id))

        if not letter_format and isinstance(s[0][0], str):
            s = convert_letter2number(s)

        pdb_matches= find_best_pdb(pfam_id, data_path, pdb_dir, n_cpu=n_cpu)
      
    
        # --------------------------------------------------------------------------------------------------------- #
        # Load PDB reference from Uniprot/Pfam to cross reference Queried PDB matches from find_best_pdb()
        try:
            # Load PDB reference from Uniprot/Pfam to cross reference Queried PDB matches from find_best_pdb()
            individual_pdb_ref_file = Path(data_path, pfam_id, 'pdb_refs.npy')
            pdb = np.load(individual_pdb_ref_file)
            # Load PDB reference file and delete 'b' in front of letters (python 2 --> python 3)
            pdb = np.array([pdb[t,i].decode('UTF-8') for t in range(pdb.shape[0]) \
                     for i in range(pdb.shape[1])]).reshape(pdb.shape[0],pdb.shape[1])
            # Create pandas dataframe for protein structure
            pdb_df = pd.DataFrame(pdb,columns = ['PF','seq','id','uniprot_start','uniprot_end',\
                                             'pdb_id','chain','pdb_start','pdb_end'])
            # Extract PDB IDs from PDB reference dataframe
            pdb_reference_ids = np.unique(pdb_df['pdb_id'].to_numpy())
            # Find PDB matches who's PDB ID also shows up in Uniprot/Pfam PDB reference file
            matches_in_reference = []
            for i, pid in enumerate(pdb_matches['PDB ID'].to_numpy()):
                if pid[:4] in pdb_reference_ids:
                    matches_in_reference.append(i)
            if len(pdb_matches.iloc[matches_in_reference]) > 0:
                pdb_reference_matches = pdb_matches.iloc[matches_in_reference]
            else:
                pdb_reference_matches = pdb_matches
        except(IndexError):
            print('Loaded pdb: ', pdb)
            print('\n\nEMPTY PDB REFERENCE!!!\n(in data_processing. so we wont have a pdb_ref match..)\n\n')
            pdb_reference_ids = None
            pdb_reference_matches = pdb_matches

         # --------------------------------------------------------------------------------------------------------- #

        
        pdb_select = pdb_reference_matches.iloc[index_pdb]
       #  pdb_select = pdb_reference_matches.loc[pdb_matches['MSA Index']==69].iloc[0]
        print(pdb_matches.head())

        pdb_id = pdb_select['PDB ID']
        original_tpdb = pdb_select['MSA Index']
        # print(original_tpdb)
        # print(ref_seq)
        for i, seq in enumerate(s):
            if not letter_format:
                # print(i, ''.join(convert_number2letter(seq)).upper())
                # print(''.join(ref_seq).upper())
                try:
                    if ''.join(convert_number2letter(seq)).upper() == ''.join(ref_seq).upper():
                        current_tpdb = i
                        break
                except(TypeError):
                    if ''.join(convert_number2letter(seq)).upper() == ''.join(convert_number2letter(ref_seq)).upper():
                        current_tpdb = i
                        break

            else: 
                if ''.join(seq) == ''.join(ref_seq):
                    current_tpdb = i
                    break
     
        print('Pre-processed MSA (shape:', s.shape, '):')
        print(s)
        print('s_index:')
        print(s_index)
        print('Matched following sequences from MSA (sequence %d) and PDB (id %s):' % (current_tpdb, pdb_id))
        print(ref_seq)
        print(pdb_seq)
        print('Reference sequence index: %d --> %d (after processing)' % (original_tpdb, current_tpdb))
        return s, removed_cols, s_index, current_tpdb, pdb_select

    # Load MSA
    s = load_msa(data_path, pfam_id)
    orig_seq_len = s.shape[1]
    print('Original Sequence length: ', orig_seq_len)

    
    # Using given MSA find best matching PDB structure from all available MSA sequences.
    
    if printing:
        print("#\n\n--------------------- Find Matching PDB Strucutre for MSA ----#")
        
    pdb_matches= find_best_pdb(pfam_id, data_path, pdb_dir, n_cpu=n_cpu)
    print(pdb_matches.head())

    
    
    # --------------------------------------------------------------------------------------------------------- #
    # re-do PDB reference file match in case it wasn't loaded ie, no if file exisists pass-through
    try:
        # Load PDB reference from Uniprot/Pfam to cross reference Queried PDB matches from find_best_pdb()
        individual_pdb_ref_file = Path(data_path, pfam_id, 'pdb_refs.npy')
        pdb = np.load(individual_pdb_ref_file)
        # Load PDB reference file and delete 'b' in front of letters (python 2 --> python 3)
        pdb = np.array([pdb[t,i].decode('UTF-8') for t in range(pdb.shape[0]) \
                 for i in range(pdb.shape[1])]).reshape(pdb.shape[0],pdb.shape[1])
        # Create pandas dataframe for protein structure
        pdb_df = pd.DataFrame(pdb,columns = ['PF','seq','id','uniprot_start','uniprot_end',\
                                         'pdb_id','chain','pdb_start','pdb_end'])
        # Extract PDB IDs from PDB reference dataframe
        pdb_reference_ids = np.unique(pdb_df['pdb_id'].to_numpy())
        # Find PDB matches who's PDB ID also shows up in Uniprot/Pfam PDB reference file
        matches_in_reference = []
        for i, pid in enumerate(pdb_matches['PDB ID'].to_numpy()):
            if pid[:4] in pdb_reference_ids:
                matches_in_reference.append(i)
        if len(pdb_matches.iloc[matches_in_reference]) > 0:
            pdb_reference_matches = pdb_matches.iloc[matches_in_reference]
        else:
            pdb_reference_matches = pdb_matches
    except(IndexError):
        print('Loaded pdb: ', pdb)
        print('\n\nEMPTY PDB REFERENCE!!!\n(in data_processing. so we wont have a pdb_ref match..)\n\n')
        pdb_reference_ids = None
        pdb_reference_matches = pdb_matches

 

    # --------------------------------------------------------------------------------------------------------- #

    # Since PDB matches are ordered by best matching.. choose first (0th) one ie index_pdb
    if len(pdb_reference_matches) == 0:
        print('No PDB matches found in MSA')
        print('should resort to pdb_refs to try as last ditch...\nIncomplete...\n\n\n\n')
        sys.exit(23)


    pdb_select = pdb_reference_matches.iloc[index_pdb]
#     # enforce old PDB refs structure for PF00186
#    pdb_select = pdb_matches.loc[pdb_matches['MSA Index']==69].iloc[0]

    ref_seq = pdb_select['Query Aligned Seq']
    pdb_seq = pdb_select['Subject Aligned Seq']
    query_seq = pdb_select['Query Aligned Seq']
    pdb_id = pdb_select['PDB ID']
    tpdb = pdb_select['MSA Index']
    if printing:
        print('Using MSA sequence %d (length %d)' % (tpdb, len(query_seq)))

    if printing:
        print("#--------------------------------------------------------------#\n\n")
    
    
    if printing:
        print("#\n\n-------------------------Remove Gaps--------------------------#")
        print('Shape of s is : ', s.shape)
        print("s = \n", s)


    gap_pdb = s[tpdb] == '-'  # returns True/False for gaps/no gaps
    s = s[:, ~gap_pdb]        # removes gaps

    s_index = np.arange(s.shape[1])

    if printing:
        print("s[tpdb] shape is ", s[tpdb].shape)
        print("s = \n", s)
        print("though s still has gaps, s[%d] does not:\n" % (tpdb), s[tpdb])
        print("s shape is ", s.shape)
        print("Saving indices of reference sequence s[%d](length=%d):\n" % (tpdb, s_index.shape[0]), s_index)
        print("#--------------------------------------------------------------#\n\n")

    lower_cols = np.array([i for i in range(s.shape[1]) if s[tpdb, i].islower()])
    if printing:
        print("removing non aligned (lower case) columns in subject sequence:\n ", lower_cols, '\n')
        
        
    # lower case removal reference: https://onlinelibrary.wiley.com/doi/full/10.1002/1097-0134%2820001101%2941%3A2%3C224%3A%3AAID-PROT70%3E3.0.CO%3B2-Z
    # upper = np.array([x.isupper() for x in s[tpdb]])
    # print('select only column presenting as uppercase at the first row')
    # upper = np.array([x.isupper() for x in s[0]])
    # s = s[:,upper]


    # --- remove duplicates before processing (as done in pydca) --- #
    if 1:
        dup_rows = []
        s_no_dup = []
        for i, row in enumerate(s):
            if [a for a in row] in s_no_dup:
                if i != tpdb:   # do not want to remove reference sequence
                    dup_rows.append(i)
                else:           # we need to add the reference sequence back in even if its a duplicate row.    
                    s_no_dup.append([a for a in row])
            else:
                s_no_dup.append([a for a in row])
        if printing:
            print('found %d duplicates! (Removing...)' % len(dup_rows))
        s, tpdb = remove_seqs_list(s, tpdb, dup_rows)
    # -------------------------------------------------------------- #

    
    # - Removing bad sequences (>gap_seqs gaps) -------------------- #    
    s, tpdb, removed_seq_indx = remove_bad_seqs(s, tpdb, gap_seqs, return_indx=True)  # removes all sequences (rows) with >gap_seqs gap %
    
    if printing:
        print('\nAfter removing bad sequences...\ntpdb (s_ipdb) is : ', tpdb)
        print(s.shape)
    # -------------------------------------------------------------- #

    # - Finding bad columns (>gap_cols) -------------------- #
    bad_cols = find_bad_cols(s, gap_cols)
    if printing:
        print('found bad columns :=', bad_cols)
    # ------------------------------------------------------ #

    # ------------ Replace Bad Amino Acid Letters if valid ones, two strategies --------- #
    if 1: # replace aa with potential correct aa
        # 2018.12.24:
        # replace 'Z' by 'Q' or 'E' with prob
        # print('replace Z by Q or E')
        s = find_and_replace(s, 'Z', np.array(['Q', 'E']))

        # replace 'B' by Asparagine (N) or Aspartic (D)
        # print('replace B by N or D')
        s = find_and_replace(s, 'B', np.array(['N', 'D']))

        # replace 'X' as amino acids with prob
        # print('replace X by other aminoacids')
        amino_acids = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', \
                                'T', 'V', 'W', 'Y', 'U'])
        s = find_and_replace(s, 'X', amino_acids)
    else: # replace aa with gap '-' (like in pydca)
        if printing:
            print('\n\nUsing PYDCA\'s method of replacing Z, B, and X with - instead of likely alternates...\n\n')
        s = find_and_replace(s, 'Z', np.array(['-']))
        s = find_and_replace(s, 'B', np.array(['-']))
        s = find_and_replace(s, 'X', np.array(['-']))
    # ----------------------------------------------------------------------------------- #

    
    # --------------------- Find Conserved Columns ------------------------------------ #
    conserved_cols = find_conserved_cols(s, conserved_cols)
    if printing:
        print("found conserved columns (80% repetition):\n", conserved_cols)
    # --------------------------------------------------------------------------------- #

    
    # ----------------------- Remove all Bad Columns ---------------------------------- #

    removed_cols = np.array(list(set(bad_cols) | set(conserved_cols)))

    removed_cols = np.array(list(set(removed_cols) | set(lower_cols)))
        
    # remove columns in which references sequence is different from matched PDB sequence
    #mismatch_cols = [i for i,a in enumerate(ref_seq) if a != pdb_seq[i]]
    #removed_cols = np.array(list(set(removed_cols) | set(mismatch_cols)))

    if printing:
        print("We remove conserved and bad columns with, at the following indices (len %d):\n" % len(removed_cols), removed_cols)
    # ----------------------- Remove all Bad Columns ---------------------------------- #


    # info still pased through removed_cols but this way we can interact with full msa if remove_cols is False
    if remove_cols:
        s = np.delete(s, removed_cols, axis=1)
        s_index = np.delete(s_index, removed_cols)

    if printing and remove_cols:
        print("Removed Columns...")
        print("s now has shape: ", s.shape)
        print("s_index (length=%d) = \n" % s_index.shape[0], s_index)
        print("In data processing Ref Seq (shape=", s[tpdb].shape, "): \n", s[tpdb])

    # Convert S to number format (number representation of amino acids)
    if not letter_format:
        # convert letter to number:
        s = convert_letter2number(s)
    print(s)

    # replace lower probs by higher probs 
    for i in range(s.shape[1]):
        s[:, i] = replace_lower_by_higher_prob(s[:, i], prob_low)

    np.save("%s/%s_pdb_query.npy" % (out_dir, pfam_id), pdb_matches)
    np.save("%s/%s_removed_cols.npy" % (out_dir, pfam_id), removed_cols)
    np.save("%s/%s_pdb_seq.npy" % (out_dir, pfam_id), pdb_seq)
    if remove_cols:
        np.save("%s/%s_preproc_msa.npy" % (out_dir, pfam_id), s)
        np.save("%s/%s_preproc_sindex.npy" % (out_dir, pfam_id), s_index)
        np.save("%s/%s_preproc_refseq.npy" % (out_dir, pfam_id), s[tpdb])
        np.save("%s/%s_preproc_pdb_seq.npy" % (out_dir, pfam_id), pdb_seq)
    else:
        np.save("%s/%s_allCols_msa.npy" % (out_dir, pfam_id), s)
        np.save("%s/%s_allCols_sindex.npy" % (out_dir, pfam_id), s_index)
        np.save("%s/%s_allCols_refseq.npy" % (out_dir, pfam_id), s[tpdb])
        np.save("%s/%s_allCols_pdb_seq.npy" % (out_dir, pfam_id), pdb_seq)



    # - Removing bad sequences (>gap_seqs gaps) -------------------- #
    if printing:
        print(s.shape)
        print("In Data Processing Final Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])

    return s, removed_cols, s_index, tpdb, pdb_select


# =========================================================================================
# process data without given pdb structure or pfam_id
def data_processing_experiment(s, pfam_id, ipdb=0, gap_seqs=0.2, gap_cols=0.2, prob_low=0.004, conserved_cols_thresh=0.8, printing=True, out_dir='./', letter_format=False, remove_cols=True,n_cpus=2):
    # def data_processing(data_path,pfam_id,ipdb=0,gap_seqs=0.2,gap_cols=0.2,prob_low=0.004):
    # read parse_pfam data:
    # print('read original aligned pfam data')
    # s = np.load('../%s/msa.npy'%pfam_id).T

    np.random.seed(123456789)

    orig_seq_len = s.shape[1]
    print('Original Sequence length: ', orig_seq_len)


    tpdb = 0 
    if printing:
        print('tpdb (s_ipdb) is : ', tpdb)

    if printing:
        print("#\n\n-------------------------Remove Gaps--------------------------#")
        print('Shape of s is : ', s.shape)
        print("s = \n", s)


    gap_pdb = s[tpdb] == '-'   # returns True/False for gaps/no gaps
    # print("removing gaps...")
    print(s[0])
    print(s[-1])
    s = s[:, ~gap_pdb]  # removes gaps
    print(s.shape)


    # print('shape of s without reference sequence gaps: ', s.shape)
    s_index = np.arange(s.shape[1])

    if printing:
        print("s[tpdb] shape is ", s[tpdb].shape)
        print("s = \n", s)
        print("though s still has gaps, s[%d] does not:\n" % (tpdb), s[tpdb])
        print("s shape is ", s.shape)
        print("Saving indices of reference sequence s[%d](length=%d):\n" % (tpdb, s_index.shape[0]), s_index)
        print("#--------------------------------------------------------------#\n\n")

    lower_cols = np.array([i for i in range(s.shape[1]) if s[tpdb, i].islower()])
    if printing:
        print("removing non aligned (lower case) columns in subject sequence:\n ", lower_cols, '\n')
    # lower case removal reference: https://onlinelibrary.wiley.com/doi/full/10.1002/1097-0134%2820001101%2941%3A2%3C224%3A%3AAID-PROT70%3E3.0.CO%3B2-Z

    # --- remove duplicates before processing (as done in pydca) --- #

    if 0:
        dup_rows = []
        s_no_dup = []
        for i, row in enumerate(s):
            if [a for a in row] in s_no_dup:
                if i != tpdb: 	# do not want to remove reference sequence
                    dup_rows.append(i)
                    print('found duplicate, row %d' % i)
                else:    	# we need to add the reference sequence back in even if its a duplicate row.	
                    s_no_dup.append([a for a in row])
             
            else:
                s_no_dup.append([a for a in row])
        np.save('%s_dup.npy' % pfam_id ,dup_rows)
        if printing:
            print('found %d duplicates! (Removing...)' % len(dup_rows))
        s, tpdb = remove_seqs_list(s, tpdb, dup_rows)
    else:
        print('bypassing duplicate search...')
    # -------------------------------------------------------------- #



    # - Removing bad sequences (>gap_seqs gaps) -------------------- #
    if printing:
        print(s.shape)
        print("In Data Processing Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])
    # print('remove sequences containing too many gaps')
    s, tpdb, bad_seq_indx = remove_bad_seqs(s, tpdb, gap_seqs, return_seq_indx=True)  # removes all sequences (rows) with >gap_seqs gap %
    if printing:
        print('\nAfter removing bad sequences...\ntpdb (s_ipdb) is : ', tpdb)
        print(s.shape)
    # -------------------------------------------------------------- #



    bad_cols = find_bad_cols(s, gap_cols)
    if printing:
        print('found bad columns :=', bad_cols)

    if 1: # replace aa with potential correct aa
        # # 2021.12.24: -- replace 'Z' by 'Q' or 'E' with prob
        # # 2023.1.6: -- using parallel find and replace (instead of serial)
        # # print('replace Z by Q or E')
        # s = find_and_replace(s, 'Z', np.array(['Q', 'E']))
        s = find_and_replace_parallel(s,'Z',np.array(['Q', 'E']), n_cpus)



        # # replace 'B' by Asparagine (N) or Aspartic (D)
        # # print('replace B by N or D')
        # s = find_and_replace(s, 'B', np.array(['N', 'D']))
        s = find_and_replace_parallel(s,'B', np.array(['N', 'D']), n_cpus)

        # # replace 'X' as amino acids with prob
        # # print('replace X by other aminoacids')
        amino_acids = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', \
                                'T', 'V', 'W', 'Y', 'U'])
        # s = find_and_replace(s, 'X', amino_acids)
        s = find_and_replace_parallel(s,'X', amino_acids, n_cpus)

    else: # replace aa with gap '-' (like in pydca)
        if printing:
            print('\n\nUsing PYDCA\'s method of replacing Z, B, and X with - instead of likely alternates...\n\n')
        s = find_and_replace(s, 'Z', np.array(['-']))
        s = find_and_replace(s, 'B', np.array(['-']))
        s = find_and_replace(s, 'X', np.array(['-']))

    # remove conserved cols
    conserved_cols = find_conserved_cols(s, fc=conserved_cols_thresh)
    if printing:
        print("found conserved columns (%f repetition):\n" % conserved_cols_thresh, conserved_cols)

    # print(s.shape)
    # print('number of conserved columns removed:',conserved_cols.shape[0])

    removed_cols = np.array(list(set(bad_cols) | set(conserved_cols)))

    removed_cols = np.array(list(set(removed_cols) | set(lower_cols)))
    if printing:
        print("We remove conserved and bad columns with, at the following indices (len %d):\n" % len(removed_cols), removed_cols)


    # info still pased through removed_cols but this way we can interact with full msa if remove_cols is False
    if remove_cols:
        s = np.delete(s, removed_cols, axis=1)
        s_index = np.delete(s_index, removed_cols)

    if printing:
        print("Removed Columns...")
        print("s now has shape: ", s.shape)
        print("s_index (length=%d) = \n" % s_index.shape[0], s_index)

    # print('replace gap(-) by other aminoacids')
    # s = find_and_replace(s,'-',amino_acids)

    if printing:
        print("In Data Processing Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])

    # Convert S to number format (number representation of amino acids)
    if not letter_format:
        # convert letter to number:
        s = convert_letter2number(s)

    # print(s.shape)
    # replace lower probs by higher probs 
    # print('replace lower probs by higher probs')
    for i in range(s.shape[1]):
        s[:, i] = replace_lower_by_higher_prob(s[:, i], prob_low)

    # print("s[tpdb] (shape=",s[tpdb].shape,"):\n",s[tpdb])
    # remove_cols = np.hstack([gap_cols,conserved_cols])
    # remove_cols = np.hstack([remove_cols,lower_cols]) ## 2019.01.22

    # np.savetxt('s0.txt',s,fmt='%i')
    # np.savetxt('cols_remove.txt',remove_cols,fmt='%i')

    # f = open('n_pos.txt','w')
    # f.write('%i'%(s.shape[1]))
    # f.close()

    # mi = number_residues(s)
    # print(mi.mean())
    np.save("%s/%s_removed_cols.npy" % (out_dir, pfam_id), removed_cols)

    # - Removing bad sequences (>gap_seqs gaps) -------------------- #
    if printing:
        print(s.shape)
        print("In Data Processing Final Reference Sequence (shape=", s[tpdb].shape, "): \n", convert_number2letter(s[tpdb]))

    return s, removed_cols, s_index, tpdb, orig_seq_len, bad_seq_indx


# =========================================================================================



def data_processing(data_path, pfam_id, ipdb=0, gap_seqs=0.2, gap_cols=0.2, prob_low=0.004, conserved_cols=0.8, printing=True, out_dir='./', letter_format=False, remove_cols=True):
    # def data_processing(data_path,pfam_id,ipdb=0,gap_seqs=0.2,gap_cols=0.2,prob_low=0.004):
    # read parse_pfam data:
    # print('read original aligned pfam data')
    # s = np.load('../%s/msa.npy'%pfam_id).T

    np.random.seed(123456789)
    s = load_msa(data_path, pfam_id)
    orig_seq_len = s.shape[1]
    print('Original Sequence length: ', orig_seq_len)


    # print('select only column presenting as uppercase at PDB sequence')
    # pdb = np.load('../%s/pdb_refs.npy'%pfam_id)
    pdb = np.load('%s/%s/pdb_refs.npy' % (data_path, pfam_id))
    # ipdb = 0

    # convert bytes to str (python 2 to python 3)
    pdb = np.array([pdb[t, i].decode('UTF-8') for t in range(pdb.shape[0]) \
                    for i in range(pdb.shape[1])]).reshape(pdb.shape[0], pdb.shape[1])
    if printing:
        print("pdb ref example (pdb[0])  (after UTF-8 decode, removing 'b'):\n", pdb[0])

    tpdb = int(pdb[ipdb, 1])
    if printing:
        print('tpdb (s_ipdb) is : ', tpdb)
    # tpdb is the sequence #
    # print(tpdb)

    if printing:
        print("#\n\n-------------------------Remove Gaps--------------------------#")
        print('Shape of s is : ', s.shape)
        print("s = \n", s)


    gap_pdb = s[tpdb] == '-'  # returns True/False for gaps/no gaps
    # print("removing gaps...")
    s = s[:, ~gap_pdb]  # removes gaps

    # print('shape of s without reference sequence gaps: ', s.shape)
    s_index = np.arange(s.shape[1])

    if printing:
        print("s[tpdb] shape is ", s[tpdb].shape)
        print("s = \n", s)
        print("though s still has gaps, s[%d] does not:\n" % (tpdb), s[tpdb])
        print("s shape is ", s.shape)
        print("Saving indices of reference sequence s[%d](length=%d):\n" % (tpdb, s_index.shape[0]), s_index)
        print("#--------------------------------------------------------------#\n\n")

    lower_cols = np.array([i for i in range(s.shape[1]) if s[tpdb, i].islower()])
    if printing:
        print("removing non aligned (lower case) columns in subject sequence:\n ", lower_cols, '\n')
    # lower case removal reference: https://onlinelibrary.wiley.com/doi/full/10.1002/1097-0134%2820001101%2941%3A2%3C224%3A%3AAID-PROT70%3E3.0.CO%3B2-Z

    # upper = np.array([x.isupper() for x in s[tpdb]])
    # print('select only column presenting as uppercase at the first row')
    # upper = np.array([x.isupper() for x in s[0]])
    # s = s[:,upper]


    # --- remove duplicates before processing (as done in pydca) --- #
    if 1:
        dup_rows = []
        s_no_dup = []
        for i, row in enumerate(s):
            if [a for a in row] in s_no_dup:
                if i != tpdb: 	# do not want to remove reference sequence
                    dup_rows.append(i)
                else:    	# we need to add the reference sequence back in even if its a duplicate row.	
                    s_no_dup.append([a for a in row])
            else:
                s_no_dup.append([a for a in row])
        if printing:
            print('found %d duplicates! (Removing...)' % len(dup_rows))
        s, tpdb = remove_seqs_list(s, tpdb, dup_rows)
    # -------------------------------------------------------------- #



    # - Removing bad sequences (>gap_seqs gaps) -------------------- #
    if printing:
        print(s.shape)
        print("In Data Processing Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])
    # print('remove sequences containing too many gaps')
    s, tpdb = remove_bad_seqs(s, tpdb, gap_seqs)  # removes all sequences (rows) with >gap_seqs gap %
    # s, tpdb, bad_seq_indx = remove_bad_seqs(s, tpdb, gap_seqs, return_seq_indx=True)  # returns original MSA indices/rows of sequences removed 
    if printing:
        print('\nAfter removing bad sequences...\ntpdb (s_ipdb) is : ', tpdb)
        print(s.shape)
    # -------------------------------------------------------------- #



    bad_cols = find_bad_cols(s, gap_cols)
    if printing:
        print('found bad columns :=', bad_cols)

    if 1: # replace aa with potential correct aa
        # 2018.12.24:
        # replace 'Z' by 'Q' or 'E' with prob
        # print('replace Z by Q or E')
        s = find_and_replace(s, 'Z', np.array(['Q', 'E']))

        # replace 'B' by Asparagine (N) or Aspartic (D)
        # print('replace B by N or D')
        s = find_and_replace(s, 'B', np.array(['N', 'D']))

        # replace 'X' as amino acids with prob
        # print('replace X by other aminoacids')
        amino_acids = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', \
                                'T', 'V', 'W', 'Y', 'U'])
        s = find_and_replace(s, 'X', amino_acids)
    else: # replace aa with gap '-' (like in pydca)
        if printing:
            print('\n\nUsing PYDCA\'s method of replacing Z, B, and X with - instead of likely alternates...\n\n')
        s = find_and_replace(s, 'Z', np.array(['-']))
        s = find_and_replace(s, 'B', np.array(['-']))
        s = find_and_replace(s, 'X', np.array(['-']))

    # remove conserved cols
    conserved_cols = find_conserved_cols(s, conserved_cols)
    if printing:
        print("found conserved columns (80% repetition):\n", conserved_cols)

    # print(s.shape)
    # print('number of conserved columns removed:',conserved_cols.shape[0])

    removed_cols = np.array(list(set(bad_cols) | set(conserved_cols)))

    removed_cols = np.array(list(set(removed_cols) | set(lower_cols)))
    if printing:
        print("We remove conserved and bad columns with, at the following indices (len %d):\n" % len(removed_cols), removed_cols)


    # info still pased through removed_cols but this way we can interact with full msa if remove_cols is False
    if remove_cols:
        s = np.delete(s, removed_cols, axis=1)
        s_index = np.delete(s_index, removed_cols)

    if printing:
        print("Removed Columns...")
        print("s now has shape: ", s.shape)
        print("s_index (length=%d) = \n" % s_index.shape[0], s_index)

    # print('replace gap(-) by other aminoacids')
    # s = find_and_replace(s,'-',amino_acids)

    if printing:
        print("In Data Processing Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])

    # Convert S to number format (number representation of amino acids)
    if not letter_format:
        # convert letter to number:
        s = convert_letter2number(s)

    # print(s.shape)
    # replace lower probs by higher probs 
    # print('replace lower probs by higher probs')
    for i in range(s.shape[1]):
        s[:, i] = replace_lower_by_higher_prob(s[:, i], prob_low)

    # print("s[tpdb] (shape=",s[tpdb].shape,"):\n",s[tpdb])
    # remove_cols = np.hstack([gap_cols,conserved_cols])
    # remove_cols = np.hstack([remove_cols,lower_cols]) ## 2019.01.22

    # np.savetxt('s0.txt',s,fmt='%i')
    # np.savetxt('cols_remove.txt',remove_cols,fmt='%i')

    # f = open('n_pos.txt','w')
    # f.write('%i'%(s.shape[1]))
    # f.close()

    # mi = number_residues(s)
    # print(mi.mean())
    np.save("%s/%s_removed_cols.npy" % (out_dir, pfam_id), removed_cols)

    # - Removing bad sequences (>gap_seqs gaps) -------------------- #
    if printing:
        print(s.shape)
        print("In Data Processing Final Reference Sequence (shape=", s[tpdb].shape, "): \n", convert_number2letter(s[tpdb]))

    return s, removed_cols, s_index, tpdb, orig_seq_len


# =========================================================================================

def data_processing_covid(data_path, pfam_id, ipdb=0, gap_seqs=0.2, gap_cols=0.2, prob_low=0.004, conserved_cols=0.8):
    # def data_processing(data_path,pfam_id,ipdb=0,gap_seqs=0.2,gap_cols=0.2,prob_low=0.004):

    printing = True
    printing = False

    # read parse_pfam data:
    # print('read original aligned pfam data')
    # s = np.load('../%s/msa.npy'%pfam_id).T
    s = np.load('%s/%s/msa.npy' % (data_path, pfam_id))
    # print(type(s))
    # print(type(s[0]))
    # print(s[0])
    if printing:
        print("shape of s (import from msa.npy):\n", s.shape)

    # convert bytes to str
    """
    try:
        s = np.array([s[t,i].decode('UTF-8') for t in range(s.shape[0]) \
             for i in range(s.shape[1])]).reshape(s.shape[0],s.shape[1])
        if printing:
    	    print("shape of s (after UTF-8 decode):\n",s.shape)
    except:
        print("\n\nUTF not decoded, pfam_id: %s \n\n"%pfam_id,s.shape)
        print("Exception: ",sys.exc_info()[0])
        # Create list file for missing pdb structures
        if not os.path.exists('missing_MSA.txt'):
            file_missing_msa = open("missing_MSA.txt",'w')
            file_missing_msa.write("%s\n"% pfam_id)
            file_missing_msa.close()
        else:
            file_missing_msa = open("missing_MSA.txt",'a')
            file_missing_msa.write("%s\n"% pfam_id)
            file_missing_msa.close()
        return
    """
    if printing:
        print('\n\nstarting shape: ', s.shape)

    if printing:
        print("#\n\n-------------------------Remove Gaps--------------------------#")
        print("s = \n", s)

    # no pdb_ref structure for covid proteins, ref strucutre is always s[0]
    tpdb = ipdb
    gap_pdb = s[tpdb] == '-'  # returns True/False for gaps/no gaps

    # print("removing gaps...")
    s = s[:, ~gap_pdb]  # removes gaps
    if printing:
        print(s.shape)
    s_index = np.arange(s.shape[1])
    print(s_index)

    if printing:
        print("s[tpdb] shape is ", s[tpdb].shape)
        print("s = \n", s)
        print("though s still has gaps, s[%d] does not:\n" % (tpdb), s[tpdb])
        print("s shape is ", s.shape)
        print("Saving indices of reference sequence s[%d](length=%d):\n" % (tpdb, s_index.shape[0]), s_index)
        print("#--------------------------------------------------------------#\n\n")

    # print(s.shape)
    # print(s)

    lower_cols = np.array([i for i in range(s.shape[1]) if s[tpdb, i].islower()])
    # print("lower_cols: ",lower_cols)

    # upper = np.array([x.isupper() for x in s[tpdb]])

    # print('select only column presenting as uppercase at the first row')
    # upper = np.array([x.isupper() for x in s[0]])
    # s = s[:,upper]
    if printing:
        print(s.shape)

    if printing:
        print("In Data Processing Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])
    # print('remove sequences containing too many gaps')
    s, tpdb = remove_bad_seqs(s, tpdb, gap_seqs)  # removes all sequences (rows) with >gap_seqs gap %
    # s, tpdb, bad_seq_indx = remove_bad_seqs(s, tpdb, gap_seqs, return_seq_indx=True)  # returns original MSA indices/rows of sequences removed 
    if printing:
        print(s.shape)

    bad_cols = find_bad_cols(s, gap_cols)
    if printing:
        print('found bad columns :=', bad_cols)

    # 2018.12.24:
    # replace 'Z' by 'Q' or 'E' with prob
    # print('replace Z by Q or E')
    s = find_and_replace(s, 'Z', np.array(['Q', 'E']))

    # replace 'B' by Asparagine (N) or Aspartic (D)
    # print('replace B by N or D')
    s = find_and_replace(s, 'B', np.array(['N', 'D']))

    # replace 'X' as amino acids with prob
    # print('replace X by other aminoacids')
    amino_acids = np.array(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', \
                            'T', 'V', 'W', 'Y', 'U'])
    s = find_and_replace(s, 'X', amino_acids)

    # remove conserved cols
    conserved_cols = find_conserved_cols(s, conserved_cols)
    if printing:
        print("found conserved columns (80% repetition):\n", conserved_cols)

    # print(s.shape)
    # print('number of conserved columns removed:',conserved_cols.shape[0])

    removed_cols = np.array(list(set(bad_cols) | set(conserved_cols)))

    removed_cols = np.array(list(set(removed_cols) | set(lower_cols)))
    if printing:
        print("We remove conserved and bad columns with, at the following indices:\n", removed_cols)

    # 2019.09.17: excluse conserved cols
    # removed_cols = np.array(list(set(bad_cols) | set(lower_cols)))

    s = np.delete(s, removed_cols, axis=1)
    s_index = np.delete(s_index, removed_cols)

    print(s_index)
    if printing:
        print("Removed Columns...")
        print("s now has shape: ", s.shape)
        print("s_index (length=%d) = \n" % s_index.shape[0], s_index)

    # print('replace gap(-) by other aminoacids')
    # s = find_and_replace(s,'-',amino_acids)

    if printing:
        print("In Data Processing Reference Sequence (shape=", s[tpdb].shape, "): \n", s[tpdb])

    # convert letter to number:
    s = convert_letter2number(s)
    # print(s.shape)
    # replace lower probs by higher probs 
    # print('replace lower probs by higher probs')
    for i in range(s.shape[1]):
        s[:, i] = replace_lower_by_higher_prob(s[:, i], prob_low)

    # print("s[tpdb] (shape=",s[tpdb].shape,"):\n",s[tpdb])
    # min_res = min_res(s)
    # print(min_res)

    # remove_cols = np.hstack([gap_cols,conserved_cols])
    # remove_cols = np.hstack([remove_cols,lower_cols]) ## 2019.01.22

    # np.savetxt('s0.txt',s,fmt='%i')
    # np.savetxt('cols_remove.txt',remove_cols,fmt='%i')

    # f = open('n_pos.txt','w')
    # f.write('%i'%(s.shape[1]))
    # f.close()

    # mi = number_residues(s)
    # print(mi.mean())
    np.save("%s/removed_cols.npy" % pfam_id, removed_cols)
    return s, removed_cols, s_index, tpdb


# =========================================================================================

def generate_pfam_data(data_path, pfam_id, ipdb):
    pdb = np.load('%s/%s/pdb_refs.npy' % (data_path, pfam_id))

    # Pre-Process Structure Data
    # delete 'b' in front of letters (python 2 --> python 3)
    print(pdb.shape)
    print(pdb)
    if len(pdb) == 0:
        print("Missing PDB structure")
        file_missing_pdb = open("missing_PDB.txt", 'a')
        file_missing_pdb.write("%s\n" % pfam_id)
        file_missing_pdb.close()
    else:
        pdb = np.array([pdb[t, i].decode('UTF-8') for t in range(pdb.shape[0]) \
                        for i in range(pdb.shape[1])]).reshape(pdb.shape[0], pdb.shape[1])

        # Create pandas dataframe for protein structure
        df = pd.DataFrame(pdb, columns=['PF', 'seq', 'id', 'uniprot_start', 'uniprot_start', \
                                        'pdb_id', 'chain', 'pdb_start', 'pdb_end'])
        print(df.head())

        # print('seq:',int(pdb[ipdb,1]))

        # try:
        # data processing
        s0, cols_removed, s_index, s_ipdb = data_processing(data_path, pfam_id, ipdb, \
                                                            gap_seqs=0.2, gap_cols=0.2, prob_low=0.004,
                                                            conserved_cols=0.9)

        # Save processed data
        msa_outfile, ref_outfile = write_FASTA(s0, pfam_id, s_ipdb, path='pfam_ecc/')
        pf_dict = {}
        pf_dict['s0'] = s0
        pf_dict['s_index'] = s_index
        pf_dict['s_ipdb'] = s_ipdb
        pf_dict['cols_removed'] = cols_removed

        with open('pfam_ecc/%s_DP.pickle' % (pfam_id), 'wb') as f:
            pickle.dump(pf_dict, f)
        f.close()
        """
        except:
            print("Could not generate data for %s: "%(pfam_id),sys.exc_info())

            if not os.path.exists('data_not_generated.txt'):
                data_fail = open("data_not_generated.txt",'w')
                data_fail.write("%s\n"% pfam_id)
                data_fail.close()
            else:
                data_fail = open("data_not_generated.txt",'a')
                data_fail.write("%s\n"% pfam_id)
                data_fail.close()
        return
        """
    return pf_dict, pdb
# -------------------------------
