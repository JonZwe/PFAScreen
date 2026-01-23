# %%
# Function to annotate peaks that have a known mass difference to annotated diagnostic fragement
# Data processing could also be done by the networkx package in the future
import numpy as np

def combine_dias_diffs(mz_array,
                       diffs, 
                       m_frags, 
                       idx_prec_diffs, 
                       idx_prec_dia_frags, 
                       idx_dia_frags, 
                       idx_differences, 
                       dist_logical, 
                       spec_formula_dias
                       ):

    # find spectra that have at least one fragment with diagnostic and difference evidence
    dia_diff_idx = np.zeros((len(m_frags), len(mz_array)))
    for idx in np.where(np.logical_and(idx_prec_diffs, idx_prec_dia_frags))[0]:
        for diff in range(len(m_frags)):
            dia_diff_idx[diff][idx] = np.sum(np.logical_and(idx_dia_frags[idx], idx_differences[diff][idx]))


    spec_idx_dia_diff = []
    new_formula_list = []
    frag_idx_list = []
    for spec in np.where(np.sum(dia_diff_idx, axis = 0))[0]:

        # array with indices of fragments that have both dias and diffs, row 1 specifies the fragment difference # ZEROS CAN BE GET IF SECOND ENTRY IS SAVED IN A VARIABLE!!!
        diff_dia_evid = np.vstack( (np.zeros(len( np.where(np.logical_and(idx_dia_frags[spec], idx_differences[0][spec]))[0])), np.where(np.logical_and(idx_dia_frags[spec], idx_differences[0][spec]))[0]) )
        frag_diff_pairs = np.array((np.where(dist_logical[0][spec])[0], np.where(dist_logical[0][spec])[1]))
        frag_diff_pairs = np.vstack((np.zeros(frag_diff_pairs.shape[1]), frag_diff_pairs))
        for i in range(1, len(m_frags)):
            diff_and_dia_evid_i = np.vstack( (i*np.ones(len(np.where(np.logical_and(idx_dia_frags[spec], idx_differences[i][spec]))[0])), np.where(np.logical_and(idx_dia_frags[spec], idx_differences[i][spec]))[0]) ) 
            diff_dia_evid = np.concatenate( (diff_dia_evid, diff_and_dia_evid_i), axis = 1)

            frag_diff_pairs_i = np.array((np.where(dist_logical[i][spec])[0], np.where(dist_logical[i][spec])[1]))
            frag_diff_pairs_i = np.vstack((i*np.ones(frag_diff_pairs_i.shape[1]), frag_diff_pairs_i))
            frag_diff_pairs = np.concatenate((frag_diff_pairs, frag_diff_pairs_i), axis = 1)

        # fragment difference pairs that also have diagnostic evidence
        frag_diff_pairs_dias = frag_diff_pairs[:, np.in1d(frag_diff_pairs[1,:], diff_dia_evid[1,:])]

        # collect neighbours, multipicator of mass difference and whether the neighbour is below or above
        # first entries are diff_and_dia_evidence and then neighours are added sequencially
        frag_idx = frag_diff_pairs_dias[2,:].astype('int')             # idx of first diffs
        partner_idx = frag_diff_pairs_dias[1,:].astype('int')          # respective partner that has diagnostic evidence
        diff_type = frag_diff_pairs_dias[0,:].astype('int')            # type of fragment difference
        diff_multiplication = np.ones(frag_diff_pairs_dias.shape[1])   # multiplication of mass diff (times 1 for next neighbour)

        # MAYBE TAKE np.unique() to increase performance
        checked_frags = frag_diff_pairs_dias[1:,:].flatten().astype('int')               # array with already checked fragments (for while loop)
        frag_diffs_comparis = frag_diff_pairs[:,np.in1d(frag_diff_pairs[1,:], frag_idx)] # find neigbours of frags which are diff*1 away from dias

        factor = 2  # start with two since differences from before have already a distance of diff*1
        while np.all( np.in1d(frag_diffs_comparis[2,:], checked_frags) ) == False:

            # find new neighbours of found diffs
            neighbour_idx = frag_diffs_comparis[2, np.in1d(frag_diffs_comparis[2,:] , checked_frags, invert = True)]
            neighbour_neighbour_idx = frag_diffs_comparis[1, np.in1d(frag_diffs_comparis[2,:] , checked_frags, invert = True)]
            diff_type_i = frag_diffs_comparis[0, np.in1d(frag_diffs_comparis[2,:] , checked_frags, invert = True)]

            p = np.zeros(len(neighbour_neighbour_idx))
            for s in range(len(neighbour_neighbour_idx)):
                p[s] = partner_idx[neighbour_neighbour_idx[s] == frag_idx][0]
            partner_idx = np.append(partner_idx, p).astype('int')

            # append data
            frag_idx = np.append(frag_idx, neighbour_idx).astype('int')
            diff_multiplication = np.append(diff_multiplication, factor * np.ones(len(neighbour_idx)))
            diff_type = np.append(diff_type, diff_type_i).astype('int')
            
            frag_diffs_comparis = frag_diff_pairs[:, np.in1d(frag_diff_pairs[1,:], frag_idx)]
            checked_frags = np.append(checked_frags, neighbour_idx)
            factor += 1

        # only if diffs without dias are present calculate new formulas 
        if np.all(np.in1d(np.unique(frag_idx), np.unique(frag_diff_pairs_dias[1,:]))) == False:
        
            formula_idx = np.vstack( ( np.where(idx_dia_frags[spec])[0] , np.arange(0, np.sum(idx_dia_frags[spec])) ) ) # indices of dia peaks with corresponding formula indices
            new_formula = [None]*len(frag_idx)
            for n in range(len(frag_idx)):

                if frag_idx[n] > partner_idx[n]:
                    if diff_multiplication[n] == 1:
                        new_formula[n] = spec_formula_dias[spec][formula_idx[1,formula_idx[0,:] == partner_idx[n]][0]] + '+' + diffs[diff_type[n]] 
                    else:
                        new_formula[n] = spec_formula_dias[spec][formula_idx[1,formula_idx[0,:] == partner_idx[n]][0]] + '+' + str(int(diff_multiplication[n])) + '(' + diffs[diff_type[n]] + ')' 
                else:
                    if diff_multiplication[n] == 1:
                        new_formula[n] = spec_formula_dias[spec][formula_idx[1,formula_idx[0,:] == partner_idx[n]][0]] + '-' + diffs[diff_type[n]] 
                    else:
                        new_formula[n] = spec_formula_dias[spec][formula_idx[1,formula_idx[0,:] == partner_idx[n]][0]] + '-' + str(int(diff_multiplication[n])) + '(' + diffs[diff_type[n]] + ')' 
        
            # Duplicates should be remove in the future to make data visualization easier!
            spec_idx_dia_diff.append(spec) 
            new_formula_list.append(new_formula)
            frag_idx_list.append(frag_idx)

    return spec_idx_dia_diff, frag_idx_list, new_formula_list