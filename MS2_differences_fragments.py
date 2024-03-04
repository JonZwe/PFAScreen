'''
FindPFAS: search for fragment mass differences, diagnostic fragments and combinations of both in centroid raw MS/MS spectra
'''
import numpy as np
import pandas as pd
from itertools import compress
from pyteomics import mass
from ismembertol import ismembertol                 # Array comparison within an absolute tolerance
from combine_dias_diffs import combine_dias_diffs   # Function to annotate MS/MS neighbour peaks from diagnostic fragments

def MS2_differences_fragments(
        mz_array,
        RT_array,
        intensity_array,
        spec_mz_list,
        spec_intens_list,
        idx_MS1,
        diffs = ['CF2'],
        number_of_fragments = 1,
        mass_tolerance = 0.002,
        intensity_threshold = 5,
        adducts = 1,
        ):

    def frags_to_mass(diffs): # NOTE: MUSS AUF GUI ANGEPASST WERDEN!
        # convert fragments specified as numbers to float
        for n in range(len(diffs)):
            if (type(diffs[n]) == float) or (type(diffs[n]) == int): # ERSETZEN DURCH ISINSTANCE?
                diffs[n] = float(diffs[n])
        # create mass of difference array from chemical formulae
        m_diffs = np.array(np.zeros((len(diffs))))
        for n in range(len(diffs)):
            if isinstance(diffs[n], str):
                m_diffs[n] = mass.calculate_mass(formula = diffs[n])
            else:
                m_diffs[n] = diffs[n]
        return m_diffs

    m_diffs = frags_to_mass(diffs)

    def spectra_filtering(spec_mz_list, spec_intens_list, intensity_threshold):
        # normalize spectra to basebeak if relative intensity is desired
        # if rel_intens == True:
        #     for n in range(len(spec_intens_list)):
        #         spec_intens_list[n] = spec_intens_list[n]/max(spec_intens_list[n])*100

        # intensity filtering: remove intensity below set threshold
        spec_intens_list_fil = [None]*len(spec_intens_list)
        spec_mz_list_fil = [None]*len(spec_mz_list)
        for n in range(len(spec_intens_list)):
            keep_intens = spec_intens_list[n] > intensity_threshold
            spec_intens_list_fil[n] = spec_intens_list[n][keep_intens]
            spec_mz_list_fil[n] = spec_mz_list[n][keep_intens]

        return spec_mz_list_fil, spec_intens_list_fil

    spec_mz_list_corr_fil, spec_intens_list_corr_fil = spectra_filtering(
                                                            spec_mz_list, 
                                                            spec_intens_list,
                                                            intensity_threshold)

    def calc_dist_matrix(spec_mz_list):
        # create a list of matrices with n x n dimensions of each intensity array
        diff_matrix = [None]*len(spec_mz_list)
        for n in range(len(spec_mz_list)):
            diff_matrix[n] = np.array(np.zeros((len(spec_mz_list[n]),len(spec_mz_list[n]))))
            XX, YY = np.meshgrid(spec_mz_list[n], spec_mz_list[n])
            # fill list with the distances between fragments (n x n matrix)
            diff_matrix[n] = abs(XX-YY)
        return diff_matrix

    def calc_intens_min(spec_intens_list):
        # create empty array for intens minimum matrices
        intens_min = [None]*len(spec_intens_list)
        for n in range(len(spec_intens_list)):
            intens_min[n] = np.array(np.zeros((len(spec_intens_list[n]),len(spec_intens_list[n]))))
            # fill array with the min intensity of two corresponding fragments
            NN, MM = np.meshgrid(spec_intens_list[n], spec_intens_list[n])
            intens_min[n] = np.minimum(np.absolute(NN),np.absolute(MM))
        return intens_min

    diff_matrix = calc_dist_matrix(spec_mz_list_corr_fil)


    def find_mass_differences(diff_matrix, spec_mz_list_corr_fil, spec_intens_list_corr_fil, m_diffs, mass_tolerance):
        # lower and upper boundary for difference between fragments
        l_bound = m_diffs - mass_tolerance
        u_bound = m_diffs + mass_tolerance

        # create empty list for logicals
        dist_logical = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_diffs))]
        idx_differences = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_diffs))]

        # loop to find index of fragment with mass differences
        for n in range(len(m_diffs)):
            for m in range(len(spec_mz_list_corr_fil)):
                # find desired differences in difference matrices 
                dist_logical[n][m] = np.logical_and(diff_matrix[m] >= l_bound[n], diff_matrix[m] <= u_bound[n]).astype(int)
                # create empty index arrays and fill with indices 
                idx_differences[n][m] = np.array(np.zeros(((len(spec_intens_list_corr_fil[m])))))
                idx_differences[n][m] = np.array(np.sum(dist_logical[n][m], axis = 0))
                idx_differences[n][m] = idx_differences[n][m] > 0
        return idx_differences, dist_logical

    idx_differences, dist_logical = find_mass_differences(diff_matrix, 
                                                          spec_mz_list_corr_fil,
                                                          spec_intens_list_corr_fil,
                                                          m_diffs,
                                                          mass_tolerance)


    def get_diff_peaks(spec_mz_list_corr_fil, spec_intens_list_corr_fil, m_diffs, idx_differences):
        # create empty lists: len(m_diffs) x len(spectra_list)
        spec_mz_diffs = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_diffs))]
        spec_intens_diffs = [[None for x in range(len(spec_mz_list_corr_fil))] for x in range(len(m_diffs))]
        n_diffs = np.zeros((len(m_diffs), len(spec_mz_list_corr_fil)))
        sum_int_weighted = np.zeros((len(n_diffs), len(spec_mz_list_corr_fil)))

        for n in range(len(m_diffs)):
            for m in range(len(spec_mz_list_corr_fil)):
                # zeros array with the length of the positive peaks at the respective positions
                spec_mz_diffs[n][m] =  np.array(np.zeros((sum(idx_differences[n][m]))))
                spec_intens_diffs[n][m] = np.array(np.zeros((sum(idx_differences[n][m]))))
                # read out the spectra and save positive m/z's and intensities in new lists
                spec_mz_diffs[n][m] = spec_mz_list_corr_fil[m][idx_differences[n][m]]
                spec_intens_diffs[n][m] = spec_intens_list_corr_fil[m][idx_differences[n][m]]
                # create list with number of fragments for each precursor m/z
                n_diffs[n][m] = len(spec_mz_diffs[n][m])
                # sum of weighted intensities for all fragment differences (NOT IN USE YET)
                sum_int_weighted[n][m] = sum(spec_intens_diffs[n][m]/len(spec_intens_diffs[n][m]))
        
        return spec_mz_diffs, spec_intens_diffs, n_diffs, sum_int_weighted

    spec_mz_diffs, spec_intens_diffs, n_diffs, _ = get_diff_peaks(spec_mz_list_corr_fil, 
                                                                 spec_intens_list_corr_fil, 
                                                                 m_diffs, 
                                                                 idx_differences)

    # sum up the total number of fragments for each precursor m/z
    n_diff_tot = np.sum(n_diffs, axis = 0)
    idx_prec_diffs = n_diff_tot >= number_of_fragments  # find all precursors which have n positive fragments
    print(str(sum(idx_prec_diffs)) + ' of ' + str(len(mz_array)) + ' MSMS spectra have fragment differences')


    def get_diagnostic_fragments(spec_mz_list_corr_fil, adducts):
        # read in DF xlsx
        dia_frags = pd.read_csv('diagnostic_fragments.csv')

        if adducts == 1: # check which polarity is desired
            dia_frags = dia_frags[dia_frags['polarity'] == 'neg']
        else:
            dia_frags = dia_frags[dia_frags['polarity'] == 'pos']

        dia_frags['F_containing'] = dia_frags['fragment_formula'].str.contains('F')
        dia_frags['fragment_formula'] = dia_frags['fragment_formula'].str.replace(' ', '', regex=False)
        dia_frags['fragment_formula'] = dia_frags['fragment_formula'].str.replace('-', '', regex=False)
        dia_frags['fragment_formula'] = dia_frags['fragment_formula'].str.replace('+', '', regex=False)
        dia_frags = dia_frags.drop_duplicates(subset = 'fragment_formula').reset_index(drop = True)
        dia_frags = dia_frags.sort_values(by=['fragment_mass']).reset_index(drop = True)
        
        # NOTE: dia_frags has to be sorted according to increasing mass!
        # NOTE: REMOVE MASS DUPLICATES!!! OR CLUSTER TOGETHER   
        remove_duplicates = True
        if remove_duplicates == True:
            decimals_round = 2 # Extremly important, depends on sample
            _, idx_unique, = np.unique(np.round(dia_frags['fragment_mass'], decimals_round), return_index=True)

            uniques = np.in1d(np.arange(0, len(dia_frags)), idx_unique)
            dia_frags = dia_frags[uniques]
        # NOTE: Non fluorine containing DFs should be marked so that precursors are only found if at least n F-containing fragments are present

        # create None lists
        idx_dia_frags = [None]*len(spec_mz_list_corr_fil)
        idx_dia_formula = [None]*len(spec_mz_list_corr_fil)
        spec_mz_dias = [None]*len(spec_mz_list_corr_fil)
        spec_intens_dias = [None]*len(spec_mz_list_corr_fil)
        spec_formula_dias = [None]*len(spec_mz_list_corr_fil)
        spec_formula_mass = [None]*len(spec_mz_list_corr_fil)
        n_dia_frags_f_contain = np.array(np.zeros((len(spec_mz_dias))))
        n_dia_frags = np.array(np.zeros((len(spec_mz_dias))))

        for n in range(len(spec_mz_list_corr_fil)):
            # find diagnostic fragments in spectra
            idx_dia_frags[n], idx_dia_formula[n], _ = ismembertol(spec_mz_list_corr_fil[n], dia_frags['fragment_mass'], mass_tolerance)
            spec_mz_dias[n] = spec_mz_list_corr_fil[n][idx_dia_frags[n]]
            spec_intens_dias[n] = spec_intens_list_corr_fil[n][idx_dia_frags[n]]
            spec_formula_dias[n] = dia_frags['fragment_formula'].iloc[idx_dia_formula[n]].reset_index(drop = True)
            spec_formula_mass[n] = dia_frags['fragment_mass'].iloc[idx_dia_formula[n]].reset_index(drop = True)
            n_dia_frags_f_contain[n] = np.sum(dia_frags['F_containing'].iloc[idx_dia_formula[n]].reset_index(drop = True))
            n_dia_frags[n] = sum(idx_dia_frags[n])

        return idx_dia_frags, spec_mz_dias, spec_intens_dias, spec_formula_dias, n_dia_frags_f_contain, n_dia_frags
    idx_dia_frags, spec_mz_dias, spec_intens_dias, spec_formula_dias, n_dia_frags_f_contain, n_dia_frags = get_diagnostic_fragments(spec_mz_list_corr_fil, adducts)

    # find indices of m/z's with diagnostic fragments
    idx_prec_dia_frags = n_dia_frags_f_contain >= number_of_fragments
    n_dia_frag_prec_tot = sum(idx_prec_dia_frags)
    print(str(n_dia_frag_prec_tot) + ' of ' + str(len(mz_array)) + " MSMS spectra have diagostic fragments")


    if all(isinstance(item, str) for item in diffs) == True:
        # annotate neighbours of diagnostic fragments characterized by mass differences
        spec_idx_dia_diff, frag_idx_list, new_formula_list = combine_dias_diffs(mz_array,
                                                                                diffs, 
                                                                                m_diffs, 
                                                                                idx_prec_diffs, 
                                                                                idx_prec_dia_frags, 
                                                                                idx_dia_frags, 
                                                                                idx_differences, 
                                                                                dist_logical, 
                                                                                spec_formula_dias)

    
    # ========== Save data for positive precursors in Pandas Dataframe ======================
    idx_diff_or_dia_prec = np.logical_or(idx_prec_diffs, idx_prec_dia_frags) # either mass difference or diagnostic evidence
     
    # save positive precursors in new arrays
    mz_hit = mz_array[idx_diff_or_dia_prec]
    RT_hit = RT_array[idx_diff_or_dia_prec]
    intensity_hit = intensity_array[idx_diff_or_dia_prec]
    idx_MS1_hit = idx_MS1[idx_diff_or_dia_prec]
    n_diffs_tot_hit = n_diff_tot[idx_diff_or_dia_prec]
    n_dias_hit = n_dia_frags[idx_diff_or_dia_prec]

    # save new list with all precursors and their number of positive mass differences
    n_diff_hit = [[None] for x in range(len(m_diffs))]
    for n in range(len(m_diffs)):
        n_diff_hit[n] = list(compress(n_diffs[n], idx_diff_or_dia_prec))

    # peaks with mass differences: delete other precursor with no positive fragments
    spec_mz_diffs_hit = [[None] for x in range(len(m_diffs))]
    spec_intens_diffs_hit = [[None] for x in range(len(m_diffs))]
    for n in range(len(m_diffs)):
        spec_mz_diffs_hit[n] = list(compress(spec_mz_diffs[n], idx_diff_or_dia_prec))
        spec_intens_diffs_hit[n] = list(compress(spec_intens_diffs[n], idx_diff_or_dia_prec))

    # save spectra with diagnostic fragments in new arrays
    spec_mz_dias_hit = list(compress(spec_mz_dias, idx_diff_or_dia_prec))
    spec_intens_dias_hit =  list(compress(spec_intens_dias, idx_diff_or_dia_prec))
    spec_formula_dias_hit = list(compress(spec_formula_dias, idx_diff_or_dia_prec))

    # all MS2 peaks: delete all precursor with no positive fragments
    spec_mz_hit = list(compress(spec_mz_list_corr_fil, idx_diff_or_dia_prec))
    spec_intens_hit = list(compress(spec_intens_list_corr_fil, idx_diff_or_dia_prec))

    # create pandas DataFrame
    Df = pd.DataFrame(data = {'mz_msms': mz_hit,
                              'rt_msms': RT_hit,
                              'intensity': intensity_hit,
                              'n_diffs': n_diffs_tot_hit,
                              'n_dias': n_dias_hit,
                              'mz_peaks': spec_mz_hit, 
                              'intens_peaks': spec_intens_hit,
                              'mz_peaks_diagnostic': spec_mz_dias_hit,
                              'intens_peaks_diagnostic': spec_intens_dias_hit,
                              'formula_diagnostic': spec_formula_dias_hit})

    # convert fragment masses to string for column names, important if mass is given instead of formula
    frags_str = [None]*len(diffs)
    for n in range(len(diffs)):
        frags_str[n] = str(diffs[n])

    # create pandas DataFrame for n_diff_hit
    Df_diffs = pd.DataFrame(n_diff_hit).T # transpose
    Df_diffs.columns = frags_str

    # add prefix to fragment string
    append_str_frag = 'mz_'
    append_str_intens = 'intens_'
    spec_frag_str = [append_str_frag + sub for sub in frags_str]
    spec_intens_str = [append_str_intens + sub for sub in frags_str]

    # create pandas DataFrames for positive spectra peaks
    Df_spec_mz_diffs = pd.DataFrame(spec_mz_diffs_hit).T
    Df_spec_mz_diffs.columns = spec_frag_str
    Df_spec_intens_diffs = pd.DataFrame(spec_intens_diffs_hit).T
    Df_spec_intens_diffs.columns = spec_intens_str

    # merge DataFrames together
    Df_fin = pd.concat([Df, Df_diffs.reindex(Df.index),
                        Df_spec_mz_diffs.reindex(Df.index),
                        Df_spec_intens_diffs.reindex(Df.index)],
                        axis = 1)

    # NOTE: DAS SOLLTE NOCHMAL AUF KORREKTHEIT ÜBERPRÜFT WERDEN
    idx_dia_diff_evid = np.where(np.in1d(np.where(idx_diff_or_dia_prec)[0], spec_idx_dia_diff))[0] # index where dia evidence was gained for diffs

    D_combine_dia_diffs = pd.DataFrame(data = {'frag_idx':frag_idx_list, 'new_formulas':new_formula_list}).set_index(idx_dia_diff_evid)
    Df_fin = pd.concat([Df_fin, D_combine_dia_diffs], axis=1).set_index(idx_MS1_hit.astype('int'))

    return Df_fin
