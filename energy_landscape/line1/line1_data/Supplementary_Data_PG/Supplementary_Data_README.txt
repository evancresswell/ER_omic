FILES AND DIRECTORIES
Supplementary_Data

 Supplementary_Data_README.txt  <---this file
 orf1_FL
 orf1_cc_names
 orf1_coiled_coil
 orf1_cterm
 orf2
 retrotransposition_vector

l1pa1 and l1pa2 orf1 were retrieved from hg19, the rest from hg18 as described in the text

all files are in fasta format regardless of their file extenstion 

file name annotations:

	seqs = nucleotides; peps = amino acids
	
	cn or cns = 50% consensus sequence unless otherwise indicated

	null indicates the CG-affected codons were translated to "O", replaced by "-", see text

	CG indicates that CG-containing positions were restored to CG - see text
--------------------------------------------------------------------------------------------

special family name desingations
--------------------------------
l1pa2:
2a (416 seqs) and 2b (527 seqs), coiled coils are identical
	2a & 2b differ by one amino acid at position 14
	  l1pa2a codon 14, AAG = K
	  l1pa2b codon 14, ACG = T

	l1pa2ab files contain the combined 2a and 2b variants in one file
    ral2ab files contain 500 sequences randomly selected from l1pa2ab files
    
------------------------------------------------------------------------------
l1pa5:

5anc, 5mod - coiled coils are identical 
	see text for how these L1Pa5 subfamilies are distinguished

	5anc_cc_peps.fst 372 peps 
	5mod_cc_peps.fst 426 peps
	5am_cc_peps.fst  748 peps (combined, 5anc and 5mod sequencs)
================================================================================

FILE TREE

avf@supplementary_data$ ls -d */  top directories
	orf1_FL/	orf1_cc_names/	orf1_coiled_coil/	orf1_cterm/		orf2/

----------------------------------------------------------
+++++++
orf1_FL/
avf@orf1_FL$ ls -d */
	l1pa1/	l1pa2/	l1pa3/	l1pa4/	l1pa5/	l1pa6/	l1pa7/
	
l1pa1/
	l1hs_orf1_peps-cg-null.1.fst <- cg-null
	l1hs_orf1_peps_w_cns.fst
	l1hs_orf1_seqs_w_cns.fst

l1pa2/
	l1pa2a_orf1_peps-cg-null.1.fst
	l1pa2a_orf1_peps_w_cns.fst
	l1pa2a_orf1_seqs_w_cns.fst
	l1pa2ab_orf1_peps-cg-null.1.fst
	l1pa2b_orf1_peps-cg-null.1.fst
	l1pa2b_orf1_peps_w_cns.fst
	l1pa2b_orf1_seqs_w_cns.fst

l1pa3/
	3_FLorf1_seqs_w_cn_cg.fst
 
l1pa4/
	4_FLorf1_seqs_w_cn-cg.fst
 
l1pa5/
	5_FLorf1anc_seqs_w_cn.fst
	5_FLorf1mod_seqs_w_cn.fst

l1pa6/
	6_FLorf1_seqs_w_cns.fst
 
l1pa7/
	7_FLorf1_seqs_w_cns.fst

------------------------------------------------------------
+++++++++++++
orf1_cc_names/ 

 these files contain the chromosome coordinates for each L1 element from which the coiled coil sequence
 was extracted...most sequences are indexed by a trailing number preceded by "_"

 except for  L1Pa1, each sequence name begins with the family number (i.e. 2a, 2b, 3, 4, 5anc, 5mod, 6 or 7)
 followed by either just the bare chromosome number (e.g., the first entry in 1_cc_names.1, 8, for chr8), or by chr
 and its number, the chromsomal coordinates for the L1 element (and optionally + or -, indicating the DNA 
 strand), and usually followed by the _index number if the chromosmal coordinates are not given...in some cases the
 cluster number (e.g., cL1) is appended to the index number

 therefore the sequence names can take various forms such as:
 L1HS::chr8:105971290-105977058(-)_2; >6__chr10_15954505_15960680; >6__chr10_15954505_15960680_2;
 >5mod_chrY_405-cL1 

 the chromosomal coordinates are not included in sequence names for some files in the Supplementary Data set

 to retrieve these coordinates search the below files using the unix grep command with the family number and chr number
 and _index number 

L1 family file name			sequence name
--------- ---------			---------------------
L1Pa1     1_cc_names		>L1HS::chr8:105971290-105977058(-)_2
          1_cc_names.1		>8_105971290-105977058(-)_2

L1Pa2     2ab_cc_names.1	>2a_chr19_29247788-29253876_2
          2ab_cc_names		>2a_chr19_29247788-29253876_2
	
L1Pa3     3_cc_names.1		>3_chr1_80126885_80133004_2
          3_cc_names		>3__chr1_80126885_80133004_2

L1Pa4     4_cc_names.1		>4_chr10_20620352_20626507_2
          4_cc_names		>4__chr10_20620352_20626507_2

L1Pa5     5anc_cc_names.1	>5_chr4_110915651_110921757_2
          5anc_cc_names		>5__chr4_110915651_110921757_2

          5mod_cc_names.1	>5_chr3_100060365_100066400_2
          5mod_cc_names		>5__chr3_100060365_100066400_2

L1pa6     6_cc_names.1		>6_chr10_15954505_15960680_2
          6_cc_names		>6__chr10_15954505_15960680_2

L1Pa7     7_cc_names.1		>7_chr11_48426178_48432210_2
          7_cc_names		>7__chr11_48426178_48432210_2

----------------------------------------------------------------
++++++++++++++++
orf1_coiled_coil/

 avf@orf1_coiled_coil$ ls -d */
  l1pa1/	l1pa2/	l1pa3/	l1pa4/	l1pa5/	l1pa6/	l1pa7/

  l1pa1/
  1_cc_peps-cg-null.1.fst
  1_cc_peps.1.fst
  1_cc_seqs.1.fst
  cls1L.1.fasta

 l1pa2/
  2ab_cc_peps-cg-null.1.fst
  2ab_cc_seqs.fst
  ral2ab_cc_peps.fst
  ral2ab_cc_peps_in-cg-null.1.fst
  ral2ab_cc_seqs.fst

 l1pa3/
  3_cc_peps.fst
  3_cc_peps_a-cg-null.1.fst
  3_cc_seqs.fst
  cls123L.3.fasta
  cls1L.3.fst
  
 l1pa4/
  4_cc_peps-cg-null.1.fst
  4_cc_peps.fst
  4_cc_seqs.fst
  cls1234L.4.fasta

 l1pa5/
  5am_cc_peps.fst
  5anc_cc_peps-cg-null.1.fst
  5anc_cc_peps.fst
  5anc_cc_seqs.fst
  5mod_cc_peps-cg-null.1.fst
  5mod_cc_peps.fst
  5mod_cc_seqs.fst
  cls123L.5am.fasta
	   
 l1pa6/
  6_cc_peps-cg-null.1.fst
  6_cc_peps.fst
  6_cc_seqs.fst
  cls1234L.6.fasta

 l1pa7/
  7_cc_peps-cg-null.1.fst
  7_cc_peps.fst
  7_cc_seqs.fst
  cls123L.7.fasta
----------------------------------------------------------------
++++++++++
orf1_cterm/

  1-7_cterm_cn_peps.fasta
  3ct_seq.fst
  4ct_seq.fst
  5anc_ct.fst
  5mod_ct.fst
  6ct.fst
  7_cterm_names
  7_cterm_seqs.fst
  C_term_l1pa1.fst
  C_term_l1pa2a.fst
  C_term_l1pa2b.fst
----------------------------------------------------------------
++++
orf2/
 Jul 10 13:07 ORF2_h_orth_l1pa2-l1pa7.xlsx
 Aug  1  2018 h_1.7_orf2_cn_pep.fasta
