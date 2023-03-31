import numpy as np
from joblib import Parallel, delayed
import random
random.seed(42)  



# Proplerly defined Energy.. 1/24/2023

def col_H(i, i1i2, w, b, s):
    for ii, (i1,i2) in enumerate(i1i2):
        if i in range(i1,i2):
            break
       
    H_i = 0
    for j in range(len(s)):
        if j in range(i1,i2):
            continue
        H_i +=  w[i,j] * s[j] 
    H_i += b[i]
    return H_i

def E_parallel(i1i2, s, w, b, s2=None):
    s_len = len(s)
    
    # Caluculate Hi for every column of H, includeing bias.
    resH = Parallel(n_jobs = 20-2)(delayed(col_H)                                                   
        (i0, i1i2, w, b, s)
        for i0 in range(s_len)) 
    H_array = resH

    coshH = 0
    for i, H_i in enumerate(H_array):
            coshH += np.log(2 * np.cosh(H_i))
    
    if s2 is not None:
        sum_sigH = np.sum(np.array([sig_i * H_array[i] for i, sig_i in enumerate(s2) ]))
        E_array  = np.array([ np.log(2 * np.cosh(H_i)) for H_i in H_array]) 
        - np.array([sig_i * H_array[i] for i, sig_i in enumerate(s2) ])
    else:
        sum_sigH = np.sum(np.array([sig_i * H_array[i] for i, sig_i in enumerate(s) ]))
        E_array  = np.array([ np.log(2 * np.cosh(H_i)) for H_i in H_array]) 
        - np.array([sig_i * H_array[i] for i, sig_i in enumerate(s) ])
    E = coshH - sum_sigH
    
    return E, E_array
   
def energy_diff(i1i2, s1, s2, w, b, return_pos_array=False):
    e_diff = 0.
    s_len = len(s1)
    if s_len==1:
        s1 = s1.transpose()
    s_len = len(s1)
    


    if s_len != len(s2):
        print('sequences not comparable!')
        print(s_len)
        print(len(s2))
        sys.exit(24)
   
    E_diff, E_diff_array = E_parallel(i1i2, s1, w, b, s2=s2)
    if return_pos_array:
        return E_diff, E_diff_array
    else:
        return E_diff

def prob_mut_LD(seq, i1i2, w, b, ncpu=2):
    S = (w + w.T)/2. # notation from LD notes                                                                           
    A = (w - w.T)/2.
    # Caluculate Hi for every column of H, includeing bias.
    resH = Parallel(n_jobs = ncpu)(delayed(col_H)                                                   
        (i0, i1i2, S, b, seq)
        for i0 in range(len(seq))) 
    H_array_S = resH
    # Caluculate Hi for every column of H, includeing bias.
    resH = Parallel(n_jobs = ncpu)(delayed(col_H)                                                   
        (i0, i1i2, A, b, seq)
        for i0 in range(len(seq))) 
    H_array_A = resH
    mut_probs = []
    prob_tot = 0.
    for i in range(len(seq)):
        H_A_i = H_array_A[i]
        H_S_i = H_array_S[i]
        seq_i = seq[i]
        mut_prob =  np.exp( seq_i * H_A_i) / ( np.exp(seq_i * H_S_i) * (2 * np.cosh(H_A_i)) )
        # mut_prob =  seq_i * np.exp( H_i) / (2 * np.cosh(H_i))
        
        prob_tot *= mut_prob
        mut_probs.append(mut_prob)
        
    return prob_tot, np.array(mut_probs)

def w_seq_LD(seq, i1i2, w_in, b_in, n_iter=1000,seed=42,ncpu=2):
    random.seed(seed)
    
    if w_in.ndim > 2:
        if w_in.ndim-1 != b_in.ndim or w_in[0].shape[1] != b_in[0].shape[0]:
            sys.exit(42)
        print('alternating walk between %d w/bs' % len(w_in))
        directions = []

    else:
        print('random walk using a cluster\'s w/b')
        w = w_in
        b = b_in
         
    seq_walk = [seq]
    n_var = len(i1i2)
    for itr in range(n_iter):
        if w_in.ndim > 2:
            indx = np.random.choice(range(len(w_in)))
            w = w_in[indx]
            b = b_in[indx]
            directions.append(indx)
            
        # Randomly choose sequence mutation (weighted for LARGEST DEVIATION)
        prob_gp1, prob_array_gp1 = prob_mut_LD(seq_walk[-1], i1i2, w, b, ncpu=ncpu)
        prob_mut = prob_array_gp1/np.sum(prob_array_gp1)
        prob_mut = prob_mut.reshape(len(prob_mut))
        mut_pos_aa = np.random.choice(range(len(prob_array_gp1)), p=prob_mut)  

        # find position mut_pos_aa is in
        found = False
        for i0 in range(n_var):
            i1,i2 = i1i2[i0][0], i1i2[i0][1]
            for i, ii0 in enumerate(range(i1,i2)):
                if mut_pos_aa == ii0:
                    found = True
                    break
            if found:
                break
        
        # apply mutation
        temp_sequence = np.copy(seq_walk[-1])
        sig_section = np.zeros(i2-i1)
        sig_section[i] = 1.
        sig_section = sig_section.reshape(len(sig_section),1)
        temp_sequence[i1:i2] = sig_section
        seq_walk.append(temp_sequence)
        
    if w_in.ndim > 2:
        return seq_walk, directions
    else:
        return seq_walk
    


def prob_mut_LD_balanced(i1i2, w, b, s_init, s_cons, ncpu=2):
    small = np.exp(-20)
    # Get the probability of seeing an amino acid at each position
    # --- use LD by scaling probability by ratio of initial vs ratio of final to keep 
    #     sequence energy from going down
    # --- balance mutations by incorporating liklihood of final sequence in context of initial
    #     sequence and global (consensus, sequence background) sequence
    # --- given an initial sequence and a global consensus (w consensus or sequence background)

    # Calculate local field for initial sequence
    # Caluculate Hi for every column of H, includeing bias.
    resH = Parallel(n_jobs = 20-2)(delayed(col_H)
        (i0, i1i2, w, b, s_init)
        for i0 in range(len(s_init)))
    H_init = resH

    # Calculate local field for consensus sequence
    # Caluculate Hi for every column of H, includeing bias.
    resH = Parallel(n_jobs = 20-2)(delayed(col_H)
        (i0, i1i2, w, b, s_cons)
        for i0 in range(len(s_cons)))
    H_cons = resH

    s_ones = np.ones(s_init.shape)

    assert len(s_ones) == len(H_cons) and len(s_ones) == len(H_init)

    # get the left hand term -- local probability of aa at position #
    num_array  = np.array([np.exp(sig_i * H_init[i]) / (2 * np.cosh(H_init[i]))
                           for i, sig_i in enumerate(s_init) ])
    # we want to probability of all possible states (ie np.ones instead of actual sequnce)
    denom_array  = np.array([np.exp(sig_i * H_init[i]) / (2 * np.cosh(H_init[i]))
                           for i, sig_i in enumerate(s_ones) ])
    denom_array[denom_array ==0.] = small  # dont want to divide by zero for zero probability
    lh_array = num_array / denom_array
    # ------------------------------------------------------------- #

    # get the right hand term -- global probability of aa at position
    num_array  = np.array([sig_i * H_init[i] / (2 * np.cosh(H_cons[i]))
                           for i, sig_i in enumerate(s_init) ])
    # we want to probability of all possible states (ie np.ones instead of actual sequnce)
    denom_array  = np.array([sig_i * H_cons[i] / (2 * np.cosh(H_cons[i]))
                           for i, sig_i in enumerate(s_ones) ])
    denom_array[denom_array ==0.] = small  # dont want to divide by zero for zero probability
    rh_array = num_array / denom_array
    # ------------------------------------------------------------- #

    assert len(lh_array) == len(rh_array)
    prob_array = lh_array * rh_array
    assert len(prob_array) == len(lh_array)

    return prob_array


def w_seq_LD_balanced(s_init, s_cons, i1i2, w_in, b_in, n_iter=1000,seed=42,ncpu=2):
    random.seed(seed)

    if w_in.ndim > 2:
        if w_in.ndim-1 != b_in.ndim or w_in[0].shape[1] != b_in[0].shape[0]:
            sys.exit(42)
        print('alternating walk between %d w/bs' % len(w_in))
        directions = []

    else:
        print('random walk using a cluster\'s w/b')
        w = w_in
        b = b_in

    seq_walk = [s_init]
    n_var = len(i1i2)
    for itr in range(n_iter):
        if w_in.ndim > 2:
            indx = np.random.choice(range(len(w_in)))
            w = w_in[indx]
            b = b_in[indx]
            directions.append(indx)

        # Randomly choose sequence mutation (weighted for LARGEST DEVIATION)
        prob_array = prob_mut_LD_balanced(i1i2, w, b, s_init, s_cons, ncpu=2)
        prob_mut = prob_array/np.sum(prob_array)
        prob_mut = prob_mut.reshape(len(prob_mut))

        # keep choosing mutations until you get new sequence (or reach max choice iterations)
        mutated = False
        counter = 0
        while not mutated or counter > 1000:
            mut_pos_aa = np.random.choice(range(len(prob_array)), p=prob_mut)
            if seq_walk[-1][mut_pos_aa] != 1:
                mutated = True
            else:
                counter += 1

        # find position mut_pos_aa is in
        found = False
        for i0 in range(n_var):
            i1,i2 = i1i2[i0][0], i1i2[i0][1]
            for i, ii0 in enumerate(range(i1,i2)):
                if mut_pos_aa == ii0:
                    found = True
                    break
            if found:
                break

        # apply mutation (i.e. swap previous onehot positive for new one)
        temp_sequence = np.copy(seq_walk[-1])
        sig_section = np.zeros(i2-i1)
        sig_section[i] = 1.
        sig_section = sig_section.reshape(len(sig_section),1)
        temp_sequence[i1:i2] = sig_section.reshape((len(sig_section),))
        seq_walk.append(temp_sequence)

    if w_in.ndim > 2:
        return seq_walk, directions
    else:
        return seq_walk       
 
