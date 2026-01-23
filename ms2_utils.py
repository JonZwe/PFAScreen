import numpy as np
from matchms import Spectrum
from matchms.similarity import CosineGreedy
import numpy as np

def merge_spectra(mz_list, 
                  intensity_list, 
                  tolerance=0.005):
    """
    Merge multiple spectra by averaging m/z values within a given tolerance,
    using intensity-weighted averaging for m/z and arithmetic mean for intensities.

    Parameters:
    - mz_list: List of 1D NumPy arrays, each containing m/z values.
    - intensity_list: List of 1D NumPy arrays, each containing corresponding intensities.
    - tolerance: m/z tolerance for merging peaks.

    Returns:
    - (merged_mz, merged_intensity): Two lists with merged values.
    """
    if len(mz_list) != len(intensity_list):
        raise ValueError("The number of m/z and intensity arrays must match.")

    # ✅ Minimal fix: make sure all inputs are at least 1D so len() works
    mz_list = [np.atleast_1d(mz) for mz in mz_list]
    intensity_list = [np.atleast_1d(ints) for ints in intensity_list]

    if not mz_list or not intensity_list or not any(len(mz) for mz in mz_list):
        return [], []

    # Flatten all spectra into a single list of (m/z, intensity)
    all_mz = np.concatenate(mz_list)
    all_intensity = np.concatenate(intensity_list)

    # Sort peaks by m/z
    sorted_indices = np.argsort(all_mz)
    mz_sorted = all_mz[sorted_indices]
    intensity_sorted = all_intensity[sorted_indices]

    merged_mz = []
    merged_intensity = []

    current_group_mz = [mz_sorted[0]]
    current_group_intensity = [intensity_sorted[0]]

    for i in range(1, len(mz_sorted)):
        if mz_sorted[i] - current_group_mz[-1] <= tolerance:
            current_group_mz.append(mz_sorted[i])
            current_group_intensity.append(intensity_sorted[i])
        else:
            # Intensity-weighted m/z
            weighted_mz = np.average(current_group_mz, weights=current_group_intensity)
            avg_intensity = np.mean(current_group_intensity)

            merged_mz.append(weighted_mz)
            merged_intensity.append(avg_intensity)

            # Start new group
            current_group_mz = [mz_sorted[i]]
            current_group_intensity = [intensity_sorted[i]]

    # Add the final group
    if current_group_mz:
        weighted_mz = np.average(current_group_mz, weights=current_group_intensity)
        avg_intensity = np.mean(current_group_intensity)

        merged_mz.append(weighted_mz)
        merged_intensity.append(avg_intensity)

    return merged_mz, merged_intensity


def append_ms2_spec_to_df_align(df_alignment, 
                                df_ms2_all, 
                                mz_tol=0.005, 
                                rt_tol=10, 
                                tol_average=0.005):
    """
    Average all MS2 spectra into one and put into list corresponding to mz and rt of df_alignment (within tolerance)
    """

    ms2_specs_mz = []
    ms2_specs_ints = []
    for idx in df_alignment.index:

        idxs = np.where((df_ms2_all['prec_mz'] > df_alignment['mz'][idx] - mz_tol) & 
                        (df_ms2_all['prec_mz'] < df_alignment['mz'][idx] + mz_tol) &
                        (df_ms2_all['prec_rt'] > df_alignment['rt'][idx] - rt_tol) & 
                        (df_ms2_all['prec_rt'] < df_alignment['rt'][idx] + rt_tol))[0]

        if len(idxs) > 0:
            mzs_ms2, ints_ms2 = merge_spectra(df_ms2_all['mzs_arr_ms2'][idxs].to_list(), 
                                            df_ms2_all['ints_arr_ms2'][idxs].to_list(), 
                                            tolerance=tol_average)
            ms2_specs_mz.append(mzs_ms2)
            ms2_specs_ints.append(ints_ms2)
        else:
            ms2_specs_mz.append(np.nan)
            ms2_specs_ints.append(np.nan)

    return ms2_specs_mz, ms2_specs_ints


def match_ms2_spectra(mzs_arr_exp, ints_arr_exp,
                      mzs_arr_lib, ints_arr_lib,
                      mz_tol=0.01):
    """
    Function to return similarity match between two mass spectra
    """

    spectrum_exp = Spectrum(mz=mzs_arr_exp, intensities=ints_arr_exp)
    spectrum_lib = Spectrum(mz=mzs_arr_lib, intensities=ints_arr_lib)
    
    # Use factory to construct a similarity function
    cosine_greedy = CosineGreedy(tolerance=mz_tol)
    score = cosine_greedy.pair(spectrum_exp, spectrum_lib)

    return score['score'], score['matches']