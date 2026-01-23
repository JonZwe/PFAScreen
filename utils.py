import os
from itertools import compress
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyopenms as oms
from chemformula import ChemFormula
from pyteomics import mass
from rdkit import Chem
from rdkit.Chem import Draw
from chromatogram_utils import get_eic
from pyopenms.Constants import C13C12_MASSDIFF_U

"""
Collection of various functions
classes: Chemical (dev)
Jonathan Zweigle, 06/2025
"""

def get_all_mass_differences(arr):
    """
    Function to return all differences in an array
    """
    x, y = np.meshgrid(arr, arr)
    diff_matrix = abs(x-y)

    return diff_matrix.flatten()


def get_top_n_peaks(mzs_arr, 
                    ints_arr, 
                    n):
    """
    Returns the top n m/z values by intensity and their original indices.
    
    Parameters:
        mzs_arr (array-like): m/z values.
        ints_arr (array-like): intensity values.
        n (int): number of top peaks to return.
    
    Returns:
        top_mzs (np.ndarray): top n m/z values sorted by intensity.
        top_indices (np.ndarray): corresponding indices in original arrays.
    """
    mzs_arr = np.asarray(mzs_arr)
    ints_arr = np.asarray(ints_arr)

    if len(mzs_arr) != len(ints_arr):
        raise ValueError("mzs_arr and ints_arr must be the same length.")
    
    if n == 0:
        return np.array([]), np.array([], dtype=int)

    n_actual = min(n, len(ints_arr))
    idx = np.argsort(ints_arr)[-n_actual:][::-1]
    
    return mzs_arr[idx], idx


def ismembertol_orig(vec1, 
                     vec2, 
                     tol) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    '''
    Backup of original function to find identicals within an absolute tolerance
    input: two arrays (vec_1 and vec_2) and absolute tolerance (tol)
    returns booleans with indices of identical values 
    e.g.: vec_1[bool_1] = values in vec_1 that are identical with vec_2 within tol
    and boolean matrix diff_bool
    '''
    vec1 = np.array(vec1)
    vec2 = np.array(vec2)
    
    diff = np.abs(np.subtract(vec1, vec2[:, np.newaxis]))
    diff_bool = np.less_equal(diff, tol)
    bool1 = np.any(diff_bool, axis=0)
    bool2 = np.any(diff_bool, axis=1)

    return bool1, bool2, diff_bool


def ismembertol(vec1, 
                vec2, 
                tol, 
                chunk_size=1000):
    '''
    Memory-efficient function to find identicals within an absolute tolerance
    input: two arrays (vec_1 and vec_2) and absolute tolerance (tol)
    returns booleans with indices of identical values and match indices
    e.g.: vec_1[bool_1] = values in vec_1 that are identical with vec_2 within tol
    NOTE: TESTING REQUIRED! Is currently used for suspect screening!
    '''

    vec1 = np.array(vec1, dtype=np.float64)
    vec2 = np.array(vec2, dtype=np.float64)
    
    # Initialize result arrays
    bool1 = np.zeros(len(vec1), dtype=bool)
    bool2 = np.zeros(len(vec2), dtype=bool)
    idx_list_1 = []
    idx_list_2 = []
    
    # Process vec1 in chunks to be memory efficient (note: processing vec1, not vec2!)
    for start_idx in range(0, len(vec1), chunk_size):
        end_idx = min(start_idx + chunk_size, len(vec1))
        vec1_chunk = vec1[start_idx:end_idx]
        
        # Same computation as original: vec2[:, np.newaxis] - vec1
        # But now with chunks: vec2[:, np.newaxis] - vec1_chunk
        diff_chunk = np.abs(np.subtract(vec1_chunk, vec2[:, np.newaxis]))
        diff_bool_chunk = np.less_equal(diff_chunk, tol)
        
        # Update bool1 and bool2 for this chunk
        bool1[start_idx:end_idx] |= np.any(diff_bool_chunk, axis=0)
        bool2 |= np.any(diff_bool_chunk, axis=1)
        
        # Collect indices - these match the original matrix orientation
        chunk_idx_2, chunk_idx_1 = np.where(diff_bool_chunk)  # Note: swapped order to match original
        # Adjust chunk_idx_1 to global indices
        chunk_idx_1 += start_idx
        idx_list_1.extend(chunk_idx_1)
        idx_list_2.extend(chunk_idx_2)

    return bool1, bool2, np.array(idx_list_2), np.array(idx_list_1)  # Note: swapped order to match original


def mass_match(mass_vec_1, 
               mass_vec_2, 
               tol
            ):
    # create meshgrid
    xxmz, yymz = np.meshgrid(mass_vec_1, mass_vec_2)

    # calculate difference
    diff_mz = abs(xxmz - yymz)
    
    mass_bool = diff_mz <= tol

    # find indizes of corresponding masses
    idx_in_1 = np.where(mass_bool)[0]
    idx_in_2 = np.where(mass_bool)[1]

    return idx_in_1, idx_in_2


def calculate_mdc_mc(mz, 
                     intens_C12, 
                     intens_C13):
    """
    Function to calculate carbon number (C), mass defect (MD), MD/C, and m/C 
    """
    # estimate number of carbons per molecule
    C = intens_C13/intens_C12/0.011145
    # calculate mass defect
    MD = mz - np.round(mz.astype(float))
    # calculate m/C and MD/C
    MDC = MD/C
    mC = mz/C

    return C, MD, MDC, mC


def repeat_mass(mz, 
                repeating_unit='CF2', 
                n=10) -> list:
    """
    Returns a list of masses for the repeating unit added to a given mass.
    """
    mass_repeating_unit = oms.EmpiricalFormula(repeating_unit).getMonoWeight()
    mass_arr = [mass_repeating_unit * i + mz for i in range(1, n+1)]

    return mass_arr


def count_similar_masses(arr, 
                         tol, 
                         plotting=False, 
                         compute_pairwise=True):
    '''
    Cluster values within tolerance tol.
    If compute_pairwise=True, first compute all pairwise absolute differences
    from the input masses and cluster those differences (useful to find repeated
    mass differences).
    Returns a DataFrame with cluster_mean, cluster_std, n_cluster.
    '''
    a = np.array(arr, dtype=float)
    # drop nans
    a = a[~np.isnan(a)]

    if compute_pairwise:
        if a.size < 2:
            return pd.DataFrame(columns=['cluster_mean', 'cluster_std', 'n_cluster'])
        # compute upper-triangle pairwise absolute differences
        diffs = np.abs(a[:, None] - a[None, :])
        iu = np.triu_indices_from(diffs, k=1)
        values = diffs[iu]
    else:
        values = a

    if values.size == 0:
        return pd.DataFrame(columns=['cluster_mean', 'cluster_std', 'n_cluster'])

    # Sort the array of values/differences
    sorted_arr = np.sort(values)

    # Build clusters by adjacency within tol
    clusters = []
    current_cluster = [sorted_arr[0]]
    for val in sorted_arr[1:]:
        if abs(val - current_cluster[-1]) <= tol:
            current_cluster.append(val)
        else:
            clusters.append(current_cluster)
            current_cluster = [val]
    clusters.append(current_cluster)

    mean_cluster = np.array([np.mean(c) for c in clusters])
    std_cluster = np.array([np.std(c) for c in clusters])
    n_cluster = np.array([len(c) for c in clusters])

    df = pd.DataFrame({'cluster_mean': mean_cluster,
                       'cluster_std': std_cluster,
                       'n_cluster': n_cluster})
    df = df.sort_values('n_cluster', ascending=False).reset_index(drop=True)

    if plotting:
        df_plot = df.copy()[:30]
        df_plot['cluster_mean_str'] = df_plot['cluster_mean'].round(4).astype(str)
        plt.figure(figsize=(4, 8))
        plt.barh(df_plot['cluster_mean_str'], df_plot['n_cluster'])
        plt.xlabel('n')
        plt.ylabel('Cluster Mean')
        plt.title('Top 30 Clusters of Similar Masses/Diffs')
        plt.tight_layout()
        plt.show()

    return df


def md_cone_filter(spec_mzs,
                   ints_mzs, 
                   prec_mz, 
                   slope1 = -0.0001, 
                   intercept1 = 0.05, 
                   slope2 = 0.0008, 
                   intercept2 = -0.05):
    """
    Filter points (mz, md) between two lines defined by slope and offset,
    where both lines start from (prec_mz, prec_md).

    Line equations:
        y1 = slope1 * (x - prec_mz) + prec_md + intercept1
        y2 = slope2 * (x - prec_mz) + prec_md + intercept2
    """

    spec_mzs = np.array(spec_mzs)
    ints_mzs = np.array(ints_mzs)
    md_spec = spec_mzs - np.round(spec_mzs)
    prec_md = prec_mz - np.round(prec_mz)

    x = spec_mzs
    y = md_spec

    y1 = slope1 * (x - prec_mz) + prec_md + intercept1
    y2 = slope2 * (x - prec_mz) + prec_md + intercept2

    lower = np.minimum(y1, y2)
    upper = np.maximum(y1, y2)

    mask = (y >= lower) & (y <= upper)

    return spec_mzs[mask], ints_mzs[mask], md_spec[mask]

# np.random.seed(42)
# mz = np.random.uniform(80, 400, 2000)
# md = np.random.uniform(-0.5, 0.5, 2000)
# plt.figure()
# mz_f, md_f = md_cone_filter(mz, md, 400, 0.1, -0.0001, 0.05, 0.001, -0.05)
# plt.scatter(mz, md, s=10, color='gray', alpha=0.5)
# plt.scatter(mz_f, md_f, s=10, color='blue', alpha=0.7)
# plt.xlabel('m/z')
# plt.ylabel('MD')
# plt.show()

def chemical_formula_substractor(formula_str_1, 
                                 formula_str_2, 
                                 elements):

    """
    Calculates the difference between two chemical formulas in terms of present elements.
    Returns a string representing the difference in the formulas and a boolean indicating 
    whether the difference is a combination of elements from both formulas or not.
    NOTE: Not tested!
    """

    formula_1 = ChemFormula(formula_str_1)
    formula_2 = ChemFormula(formula_str_2)

    def element_counter(formula, elements):
        elem_count = np.zeros(len(elements))
        for n in range(len(formula)):
            for m in range(len(elements)):
                if formula[n][0] == elements[m]:
                    elem_count[m] = formula[n][1]
        return elem_count

    elem_counts_1 = element_counter(list(formula_1.element.items()), elements)
    elem_counts_2 = element_counter(list(formula_2.element.items()), elements)

    if np.sum((elem_counts_1-elem_counts_2) < 0) > 0 and np.sum((elem_counts_2-elem_counts_1) < 0) > 0:
        mixed = True
        diff_counts_1 = (elem_counts_1-elem_counts_2).astype('int')
        diff_counts_2 = (elem_counts_2-elem_counts_1).astype('int')
        elems_1 = list(compress(elements, diff_counts_1 < 0))
        elems_2 = list(compress(elements, diff_counts_2 < 0))
        counts_1 = abs(diff_counts_1[diff_counts_1 < 0])
        counts_2 = abs(diff_counts_2[diff_counts_2 < 0])
        formula_diff_1 = ''.join([n+str(m) for n,m in zip(elems_1,counts_1)])
        formula_diff_2 = ''.join([n+str(m) for n,m in zip(elems_2,counts_2)])
        if mass.calculate_mass(formula_str_1) > mass.calculate_mass(formula_str_2):
            formula_diff_tot = '-'.join([formula_diff_2, formula_diff_1])
        else:
            formula_diff_tot = '-'.join([formula_diff_1, formula_diff_2])

    else:
        mixed = False
        if mass.calculate_mass(formula_str_1) > mass.calculate_mass(formula_str_2):
            diff_counts = abs((elem_counts_1-elem_counts_2).astype('int'))
        else:
            diff_counts = abs((elem_counts_2-elem_counts_1).astype('int'))

        elems = list(compress(elements, diff_counts > 0))
        counts = diff_counts[diff_counts > 0]
        formula_diff = ''.join([n+str(m) for n,m in zip(elems,counts)])

    if np.sum((elem_counts_1-elem_counts_2) < 0) > 0 and np.sum((elem_counts_2-elem_counts_1) < 0) > 0:

        return formula_diff_tot, mixed 
    else:
        return formula_diff, mixed


def estimate_noise_raw(path, 
                       percent_spectra = 10,
                       percent_max = 25, 
                       min_mz = 100, 
                       n_bins = 100,
                       plot = False):
    
    """
    Function to estimate a suitable noise level (e.g., as a threshold for feature detection) from a raw mzML file.
    It should be noted, that the histogram stongly depends on the raw data itselve, meaning that is only gives a first idea
    of the noise level.
    """

    # load mzML file
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)

    # get the number of spectra in the file
    n_spectra_tot = exp.size()

    # get n% of spectra (from the middle of the run)
    n_spectra = int(n_spectra_tot * percent_spectra/100)
    start = int(n_spectra_tot * (0.5 - percent_spectra/100/2))
    end = int(n_spectra_tot * (0.5 + percent_spectra/100/2))

    # concatenate the intensity values of all those spectra
    mz_values_all = [np.array(exp[n].get_peaks()[0]) for n in range(start, end)]
    intensity_values_all = [np.array(exp[n].get_peaks()[1]) for n in range(start, end)]
    # remove m/z values below min_mz
    intensity_values_all = [intensity_values_all[n][mz_values_all[n] > min_mz] for n in range(len(intensity_values_all))]
    intensity_values_all = np.concatenate(intensity_values_all)

    # Sort all intensity values in ascending order
    sorted_indices = np.argsort(intensity_values_all)
    sorted_intensities = intensity_values_all[sorted_indices]

    log_intensities = np.ma.log10(sorted_intensities).compressed()

    # Histogram
    hist_y, bin_edges = np.histogram(log_intensities, bins=n_bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Get the global max of the histogram
    peak_height = np.max(hist_y)
    threshold_height = percent_max/100 * peak_height

    # Walk from left to right and find the first bin that drops below 10% of the peak
    cut_idx = np.argmax(hist_y)  # Start from the peak
    while cut_idx < len(hist_y) and hist_y[cut_idx] > threshold_height:
        cut_idx += 1

    # Use the bin center at the cutoff index
    threshold_pos = bin_centers[cut_idx] if cut_idx < len(bin_centers) else bin_centers[-1]

    noise_estimate = 10**threshold_pos

    if plot:

        plt.figure()
        # Log-transform the intensities (log10)
        log_intensities = np.ma.log10(sorted_intensities)

        # Create histogram (log-transformed intensities)
        _, _, _ = plt.hist(log_intensities, bins=n_bins, alpha=0.7, color='blue', density=True)

        plt.axvline(threshold_pos, color='black', linestyle='--', label='10% drop from peak')
        plt.scatter([threshold_pos], [hist_y[cut_idx]], color='red', zorder=5)
        # Add labels and title
        plt.title(f"Sample: {os.path.basename(path)}\nNoise Estimate: {int(noise_estimate)}")
        plt.xlabel("Log10(Intensity)")
        plt.ylabel("Density")

        plt.tight_layout()
        plt.show()
    
    return noise_estimate


def match_mz_rt_intens(mz_array_1, rt_array_1, intens_array_1, 
                       mz_array_2, rt_array_2, intens_array_2, 
                       mz_tol, 
                       rt_tol, 
                       rel_int_tol):
    """
    Matches features in two lists by m/z, RT, and intensity (with thresholds)
    Parameters
    ----------
    mz_array_1, rt_array_1, intens_array_1 : np.arrays 
        m/z, retention times, and intensities of the first feature list
    The same for the second feature list
    
    mz_tol : float
        m/z tolerance (Da)
    rt_tol : float
        retention time tolerance (min)
    rel_int_tol : float
        relative intensity tolerance (unitless, e.g. 0.1 means 10%)
    
    Returns
    -------
    idx_in_a : array-like
        indices of features in the first list that were matched to the second list
    idx_in_b : array-like
        indices of features in the second list that were matched to the first list
    total_bool : array-like
        boolean array of shape (n x n) with True where the match was successful and False otherwise
    """

    # save memory by converting to float32
    # NOTE: WARNING: THIS ONLY ALLOWS FOR A MAXIMUM FOUR DIGIT PRECISION FOR MZ up to 9999.9999
    # ALSO APPLIES TO INTENSITY! WILL BE ROUNDED!!! CHECK IF ENOUGH FOR INTENSITY COMPARISON
    mz_array_1 = mz_array_1.astype(np.float32)  
    mz_array_2 = mz_array_2.astype(np.float32)
    rt_array_1 = rt_array_1.astype(np.float32)
    rt_array_2 = rt_array_2.astype(np.float32)
    intens_array_1 = intens_array_1.astype(np.float32)
    intens_array_2 = intens_array_2.astype(np.float32)

    # create a difference and ratio matrices with n x n dimensions for mz, RT, and intens
    # create meshgrids
    xxmz, yymz = np.meshgrid(mz_array_1, mz_array_2)
    xxrt, yyrt = np.meshgrid(rt_array_1, rt_array_2)
    xxintens, yyintens = np.meshgrid(intens_array_1, intens_array_2)
    
    # calculate difference and intensity ratio
    diff_mz = abs(xxmz - yymz)
    diff_rt = abs(xxrt - yyrt)
    rel_diff_int = abs( (xxintens - yyintens)/ ((xxintens + yyintens) / 2) )
    
    mz_bool = diff_mz <= mz_tol
    rt_bool = diff_rt <= rt_tol
    intens_bool = rel_diff_int <= rel_int_tol
    
    total_bool = mz_bool.astype('int') * rt_bool.astype('int') * intens_bool.astype('int')

    idx_in_a =np.where(total_bool)[0]
    idx_in_b = np.where(total_bool)[1]
    
    return idx_in_a, idx_in_b, total_bool


def estimate_c(exp, 
               mz_arr, 
               rt_arr, 
               peak_width_arr=None, 
               rt_width=10, 
               extraction_width=0.005):

    print(f'{len(mz_arr)*2} EICs to extract!')

    mz_13c_arr = mz_arr + C13C12_MASSDIFF_U  # add 13C mass difference to the m/z values

    area_mz_arr = np.zeros(len(mz_arr))
    area_mz_13c_arr = np.zeros(len(mz_arr))
    for n in range(len(mz_arr)):

        if peak_width_arr is not None:

            _, ints_m = get_eic(exp, mz_arr[n], rt_arr[n], peak_width_arr[n], extraction_width)
            _, ints_m1 = get_eic(exp, mz_13c_arr[n], rt_arr[n], peak_width_arr[n], extraction_width)
        else:
            _, ints_m,  _ = get_eic(exp, mz_arr[n], rt_arr[n], rt_width, extraction_width)
            _, ints_m1, _ = get_eic(exp, mz_13c_arr[n], rt_arr[n], rt_width, extraction_width)
        
        area_mz_arr[n] = np.trapz(ints_m)
        area_mz_13c_arr[n] = np.trapz(ints_m1)

        # estimate carbon number
        c = area_mz_13c_arr/area_mz_arr/0.011145
            
    return c, area_mz_arr, area_mz_13c_arr


def fraction_with_close_neighbors(arr, 
                                  threshold=0.005):

    """
    Function to count close m/z values (with a given tolerance) 
    """

    arr = np.sort(arr)
    n = len(arr)
    has_neighbor = np.zeros(n, dtype=bool)
    
    i = 0
    for j in range(1, n):
        while arr[j] - arr[i] > threshold:
            i += 1
        if j > i:
            has_neighbor[i:j+1] = True  # mark all in range as having close neighbors

    return np.mean(has_neighbor)  # fraction of elements with close neighbors


def blank_correc(mass_vec, 
                 RT_vec, 
                 intens_vec, 
                 mass_vec_B, 
                 RT_vec_B, 
                 intens_vec_B, 
                 m_tol, RT_tol, 
                 fold_change
                 ):
    '''
    Basic function to perform a blank correction based on accurate mass- and RT-comparision and fold change
    Input data: np.array of mass, RT, intensity for both sample and blank
    m_tol: absolute mass tolerance, RT_tol: retention time tolerance, fold_change: desired fold change
    Output data: boolean indices of overlapping features (not in blank and inverted: in blank)
    '''

    # create a difference and ratio matrices with n x n dimensions for mz, RT, and intens
    # create meshgrids
    xxmz, yymz = np.meshgrid(mass_vec, mass_vec_B)
    xxRT, yyRT = np.meshgrid(RT_vec, RT_vec_B)
    xxintens, yyintens = np.meshgrid(intens_vec, intens_vec_B)
    
    # calculate difference and intensity ratio
    diff_mz = abs(xxmz - yymz)
    diff_RT = abs(xxRT - yyRT)
    intens_ratio = xxintens/yyintens
    
    mass_bool = diff_mz <= m_tol
    RT_bool = diff_RT <= RT_tol
    
    if type(fold_change) == int or type(fold_change) == float:
        intens_bool = intens_ratio <= fold_change
        total_bool = mass_bool.astype('int') * RT_bool.astype('int') * intens_bool.astype('int')
    else:
        total_bool = mass_bool.astype('int') * RT_bool.astype('int')

    # find indizes of features in blank
    idx_in_blank =  np.sum(total_bool, axis = 0) > 0 # index of features that occur in the blank
    idx_not_in_blank = idx_in_blank == 0             # features that are unique to the sample 
    
    return idx_not_in_blank, idx_in_blank, total_bool


def bin_ms2_to_matrix(mz_vec_corr, 
                      TIC_vec_corr, 
                      spec_mz_list_corr_fil, 
                      spec_intens_list_corr_fil, 
                      idx_prec, 
                      TIC_min, 
                      MD_filter, 
                      MD_low, 
                      MD_upp, 
                      tol_mat):

    """
    Function to bin MS2 specta to one homogenous matrix for further data
    evaluation such as e.g. cosine similarity or PCA
    J.Z., May 2022

    # Parameters:
    # tol_mat: absolute tolerance for bin size
    # TIC_min: threshold of precursor intensity below which MS2 specs should be removed
    # MD_filter: True/False: choose whether data should be filtered base on mass defect of precursor m/z
    # MD_upp, MD_low: upper and lower mass defect limit

    """

    # =========================================================================
    # if desired set TIC limit to further reduce data
    idx_TIC = TIC_vec_corr > TIC_min

    # delete entries below set TIC value
    idx_prec = idx_prec[idx_TIC]
    mz_vec_corr = mz_vec_corr[idx_TIC]
    TIC_vec_corr = TIC_vec_corr[idx_TIC]
    spec_mz_list_corr_fil = list(compress(spec_mz_list_corr_fil, idx_TIC))
    spec_intens_list_corr_fil = list(compress(spec_intens_list_corr_fil, idx_TIC))

    # check if spectra are normalized, if not normalize
    all_intens = np.concatenate(spec_intens_list_corr_fil, axis = 0)
    if max(all_intens) > 100:
        for n in range(len(spec_intens_list_corr_fil)):
            try:
                spec_intens_list_corr_fil[n] = spec_intens_list_corr_fil[n]/max(spec_intens_list_corr_fil[n])*100
            except ValueError:
                spec_intens_list_corr_fil[n] = spec_intens_list_corr_fil[n]


    # =========================================================================
    if MD_filter == True:

        # calculate mass defect
        MD = mz_vec_corr - np.round_(mz_vec_corr, decimals=0)
        # find index of mass defect within tolerance
        idx_MD = np.logical_and(MD >= MD_low, MD <= MD_upp)

        # read out arrays
        mz_vec_corr = mz_vec_corr[idx_MD]
        TIC_vec_corr = TIC_vec_corr[idx_MD]
        spec_mz_list_corr_fil = list(compress(spec_mz_list_corr_fil, idx_MD))
        spec_intens_list_corr_fil = list(compress(spec_intens_list_corr_fil, idx_MD))
        idx_prec = idx_prec[idx_MD]

    # =========================================================================
    # create array with all fragment masses
    all_frags = np.concatenate(spec_mz_list_corr_fil, axis = 0)

    # find min and max fragments add certain mass to guarantee that masses are within bins
    max_fragment = max(all_frags)+0.01
    min_fragment = min(all_frags)-0.01

    # =========================================================================
    # calculate number of bins with certain mass tolerance
    n_bins = int((max_fragment-min_fragment)/tol_mat)

    # create bins with masses
    bins = np.linspace(min_fragment, max_fragment, int(n_bins))

    # create matrix with rows (number of specta) and columns (number of bins)
    X = np.array(np.zeros((len(spec_intens_list_corr_fil), len(bins))))
    for n in range(len(spec_mz_list_corr_fil)):
        vec = np.array(np.zeros(len(bins))) # create empty vector with length of bins
        digitized = np.digitize(spec_mz_list_corr_fil[n], bins) # get indized of bin positions where the fragment mass belongs
        vec[digitized] = spec_intens_list_corr_fil[n] # write normalized intensity in vector at the position of the mass bin
        X[n,:] = vec # write vector at position in X matrix


    mass_idx = ~np.all(X == 0, axis=0)
    mass_array = bins[mass_idx]
    # remove all bins that contain no fragments (= 0)
    X = X[:,~np.all(X == 0, axis=0)]

    return X, mass_array


class Chemical():
    """
    A class to represent a chemical compound.
    """
    def __init__(self, formula, smiles):
        self.formula = formula
        self.smiles = smiles
        self.mass = self.get_accurate_mass()
        self.C = self.get_carbon_number()
        self.MD = self.get_MD()
    def get_carbon_number(self):
        return oms.EmpiricalFormula(self.formula).getElementalComposition()[b'C']
    def get_accurate_mass(self):
        isotopes = oms.EmpiricalFormula(self.formula).getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(4))
        return isotopes.getContainer()[0].getMZ()
    def get_MD(self):
        return self.mass - round(self.mass)
    def get_mC(self):
        return self.mass/self.C
    def show_structure(self):
        mol = Chem.MolFromSmiles(self.smiles)
        im = Draw.MolToImage(mol)
        plt.figure()
        plt.imshow(im)
        plt.axis('off')
    def show_isotope_pattern(self):
        isotopes = oms.EmpiricalFormula(self.formula).getIsotopeDistribution(oms.CoarseIsotopePatternGenerator(4))
        mzs = [iso.getMZ() for iso in isotopes.getContainer()]
        ints = [iso.getIntensity() for iso in isotopes.getContainer()]
        plt.figure()
        plt.stem(mzs, ints, 'Black', markerfmt=" ", basefmt=" ")
        _, stemlines, _ = plt.stem(mzs, ints, 'Black',markerfmt=" ", basefmt=" ")
        plt.setp(stemlines, color = 'Black', linewidth= 2)
        for i, txt in enumerate(np.round(mzs, 4)):
            plt.annotate(txt, (mzs[i],ints[i]), color = 'Black', rotation = 20, fontsize = 10)
        plt.ylim(ymin = 0)
        plt.xlabel('m/z')
        plt.ylabel('Relative intensity')

# glyphosate = Chemical('C3H8NO5P','C(C(=O)O)NCP(=O)(O)O')