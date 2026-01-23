
import numpy as np
import pyopenms as oms

def match_isotope_patters_raw(exp, 
                              formula, 
                              element, 
                              element_bool, 
                              rt, 
                              md_deviation = 0.02):

    """
    Function to match a given formula to a measured isotope 
    pattern in a oms.MSExperiment at a given retention time
    """

    if element_bool == '+':
        formula = formula + element
    elif element_bool == '-':
        element_counts = oms.EmpiricalFormula(formula).getElementalComposition()
        # make a compatable binary string for openms b'H == 'H'.encode()
        if element.encode() in element_counts:
            element_counts[element.encode()] -= 1
            formula = ''.join([key.decode() + str(value) for key, value in element_counts.items()])


    # Theoretical isotopes
    formula = oms.EmpiricalFormula(formula)
    isotopes = formula.getIsotopeDistribution( oms.CoarseIsotopePatternGenerator(4) )
    mz_arr_theo = np.array([iso.getMZ() for iso in isotopes.getContainer()])
    ints_arr_theo = np.array([iso.getIntensity() for iso in isotopes.getContainer()])
    ints_arr_theo = ints_arr_theo/np.max(ints_arr_theo)*100

    # Extraction from measured raw data
    rt_array = np.array([spec.getRT() for spec in exp])
    # closest index to the target retention time
    idx = np.argmin(np.abs(rt_array - rt))
    # array with MS levels
    ms_level_arr = np.array([spec.getMSLevel() for spec in exp])
    # find indices where the MS1 level is 1
    indices_ms1 = np.where(ms_level_arr == 1)[0]
    # compute the distances to the target index
    distances = np.abs(indices_ms1 - idx)   
    # find the closest 1 (ties are automatically resolved)
    idx = indices_ms1[np.argmin(distances)]

    mz_array = exp[int(idx)].get_peaks()[0]
    ints_array = exp[int(idx)].get_peaks()[1]

    idx_cut = np.logical_and(mz_array > mz_arr_theo[0]-0.1, mz_array < mz_arr_theo[-1]+0.1)

    # check if there are m/z present after removing below and above the isotopes (can be empty if scan is sparse)
    if np.sum(idx_cut) > 0:
        
        mz_arr_exp = mz_array[idx_cut]
        ints_arr_exp = ints_array[idx_cut]/np.max(ints_array[idx_cut])*100

        # get closest isotopes to perform matching
        md_precursor = mz_arr_theo[0] - np.round(mz_arr_theo[0])
        md_mz_arr_exp = mz_arr_exp - np.round(mz_arr_exp)

        # max MD deviations: 3*13C = 0.01008, 37Cl = -0.0029, 81Br = -0.0020, 34S = -0.0042
        # assume 0.02 Da deviation to be conservative enough
        idx_md_match = np.abs(md_mz_arr_exp - md_precursor) < md_deviation
        mz_arr_exp_md_cleaned = mz_arr_exp[idx_md_match]
        ints_arr_exp_md_cleaned = ints_arr_exp[idx_md_match]

        # check if cleaned up array is empty
        if len(mz_arr_exp_md_cleaned) > 0:

            # find closest mz
            mz_precursor_exp = mz_arr_exp_md_cleaned[np.argmin(abs(mz_arr_exp_md_cleaned - mz_arr_theo[0]))]
            ppm_error = ((mz_precursor_exp - mz_arr_theo[0])/mz_arr_theo[0])*1e6

            # find closest isotopes: remove mz's that deviate > 0.1 Da to ensure that wrong ones are removed
            # NOTE: consider to improve: could be done like in MZmine, where intensities are summed (than normalized again), 
            # this minimizes the risk that very low abundant ions are erroneously considered
            mz_arr_exp_score = []
            ints_arr_exp_score = []
            for mz_theo in mz_arr_theo:
                
                idx_closest = np.argmin(abs(mz_arr_exp_md_cleaned - mz_theo))

                if abs(mz_arr_exp_md_cleaned[idx_closest] - mz_theo) < 0.1:
                    mz_arr_exp_score.append(mz_arr_exp_md_cleaned[idx_closest])
                    ints_arr_exp_score.append(ints_arr_exp_md_cleaned[idx_closest])
                else:
                    mz_arr_exp_score.append(np.nan)
                    ints_arr_exp_score.append(0)

            # NOTE: Consider also normalizing the intensities here again!
            # Can be a bug!

            mz_arr_exp_score = np.array(mz_arr_exp_score)
            ints_arr_exp_score = np.array(ints_arr_exp_score)


            # score similar to the isotope score in MZMine https://doi.org/10.1021/ac3000418
            def isotope_score(intens_array_theoretical, intens_array_measured):

                intens_difference = abs(intens_array_theoretical - intens_array_measured)
                score = 1
                for n in range(len(intens_difference)):
                    score *= (1 - intens_difference[n])
                    
                return score
            
            iso_score = isotope_score(ints_arr_theo/np.max(ints_arr_theo), ints_arr_exp_score/np.max(ints_arr_exp_score))
        else:
            print('Isotope matching failed due to no m/z found after cleaning up!')

            mz_arr_exp = np.nan
            ints_arr_exp = np.nan
            mz_arr_exp_score = np.nan
            ints_arr_exp_score = np.nan
            iso_score = np.nan
            ppm_error = np.nan
    else:
        print('Isotope matching failed due to missing data!')
        
        mz_arr_exp = np.nan
        ints_arr_exp = np.nan
        mz_arr_exp_score = np.nan
        ints_arr_exp_score = np.nan
        iso_score = np.nan
        ppm_error = np.nan

    isotope_dict = {'mz_arr_theo': mz_arr_theo, 
                    'ints_arr_theo': ints_arr_theo, 
                    'mz_arr_exp': mz_arr_exp, 
                    'ints_arr_exp': ints_arr_exp,
                    'mz_arr_exp_score': mz_arr_exp_score,
                    'ints_arr_exp_score': ints_arr_exp_score,
                    'ppm_error': ppm_error, 
                    'iso_score': iso_score}

    return isotope_dict