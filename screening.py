
import os
import time
import numpy as np
import pandas as pd
from tqdm import tqdm

from ms2_differences_fragments import ms2_differences_fragments # FindPFAS: search in MS/MS raw spectra for PFAS
from kmd_analysis import kmd_analysis                           # Kendrick mass defect analysis
from suspect_screening import suspect_screening                 # Basic suspect screening
from get_adduct_data import get_adduct_data
from utils import get_top_n_peaks
from match_subformulas_to_peaks import match_subformulas_to_peaks

def pfascreen(df,
              sample_names,
              params
              ):

    start = time.time()

    # Convert config object to dict if needed (for backward compatibility)
    if hasattr(params, 'to_dict'):
        params = params.to_dict()

    # read out params
    df_suspects = pd.read_csv(params['path_suspect_list'])
    df_diagnostic_fragments = pd.read_csv(params['path_diagnostic_fragments'])
    df_fragment_differences = pd.read_csv(params['path_mass_differences'])
    
    number_of_fragments = params['number_of_fragments']
    mass_tolerance = params['mass_tolerance']
    intensity_threshold = params['intensity_threshold_ms2']
    diff_kmd = params['diff_kmd']
    n_homologues = params['n_homologues']
    polarity = params['polarity']

    # calculate mean of sample
    df['intens_mean'] = df[sample_names].mean(axis = 1)
    
    # write in adduct mass difference
    # NOTE: here is a bug in case formulas like [M+2H] are present!
    adduct_mass_diff = {add: get_adduct_data(add)[0] for add in df['adduct'].unique()}
    df['adduct_mass_diff'] = df['adduct'].map(adduct_mass_diff)

    # NOTE: CAUTION -> This is only a temporary solution to avoid having 
    # wrong isotope pattern matches when searching for M+ or M- during suspect screening
    #if adducts == 3:
    #    df['adduct_mass_diff'] = 0

    # considers only MS2 data
    df_findpfas = ms2_differences_fragments(
                    df['mz'].to_numpy(),
                    df['rt'].to_numpy(),
                    np.ones(len(df)),
                    np.array(df.index),
                    df['mzs_ms2'].to_list(),
                    df['ints_ms2'].to_list(),
                    df_fragment_differences,
                    df_diagnostic_fragments,
                    number_of_fragments=number_of_fragments, # factor of two for diagnostic fragments
                    mass_tolerance=mass_tolerance, 
                    intensity_threshold=intensity_threshold,
                    polarity=polarity
                    )

    df = pd.concat([df, df_findpfas], axis = 1) # NOTE: Works only of only unique MS1 data is present!!! (see line before FindPFAS)

    # write zeros at columns where MSMS spectra are present but not diffs or dias were found
    idx_MSMS_all = np.where(df['mzs_ms2'].notna())[0] # indices with MSMS
    idx_no_hit = df.index[df['n_diffs'].isnull()].to_numpy()
    idx_no_MSMS_hit = idx_no_hit[np.in1d(idx_no_hit, idx_MSMS_all)]
    df.loc[idx_no_MSMS_hit, 'n_diffs'] = 0
    df.loc[idx_no_MSMS_hit, 'n_dias'] = 0

    #%%
    # MD/C-m/C, and MD
    # ==============================================================================================
    def calc_MDC_mC(mz, intens_C12, intens_C13):
        # Avoid divide-by-zero or NaN by replacing 0 or very small values
        intens_C12 = np.where(intens_C12 <= 0, np.nan, intens_C12)
        intens_C13 = np.where(intens_C13 < 0, np.nan, intens_C13)  # Negative intensity makes no sense

        # estimate number of carbons per molecule
        C = intens_C13 / intens_C12 / 0.011145

        # calculate mass defect
        MD = mz - np.round(mz, decimals=0)

        # Avoid divide-by-zero on C
        C = np.where((C == 0) | np.isnan(C), np.nan, C)

        # calculate m/C and MD/C
        mC = mz / C
        MDC = MD / C

        # specify constants
        mC_CF2 = 49.996806
        MDC_CF2 = -0.003194
        m_CF = -8.40596e-05
        m_CHF = -0.0005237

        # transform data to origin and rotate by the slope of the CHF line
        cos_theta = np.cos(-m_CHF)
        sin_theta = np.sin(-m_CHF)

        mC_sr =  (mC - mC_CF2) * cos_theta - (MDC - MDC_CF2) * sin_theta
        MDC_sr = (mC - mC_CF2) * sin_theta + (MDC - MDC_CF2) * cos_theta

        # calculate radial distance from the CF2 location (with ellipse factor lambda)
        r_CF2 = np.sqrt((mC_sr / 3000) ** 2 + (MDC_sr) ** 2)

        return C, MDC, mC, mC_sr, MDC_sr, r_CF2

    # Create masks for valid isotopic intensity rows
    # NOTE: Needs to be checked!
    valid_mask = df['ints_isotopes'].apply(
        lambda x: isinstance(x, (list, tuple)) 
            and len(x) >= 2 
            and all(pd.notna(x)) 
            and x[1] > 0)

    # Initialize new columns with NaNs
    df['C'] = np.nan
    df['MD/C'] = np.nan
    df['m/C'] = np.nan
    df['r_CF2'] = np.nan

    # Only calculate for valid rows
    if valid_mask.any():
        valid_df = df[valid_mask]
        intens_C12 = np.array([n[0] for n in valid_df['ints_isotopes']])
        intens_C13 = np.array([n[1] for n in valid_df['ints_isotopes']])

        # Run your calculation
        C, MDC, mC, _, _, r_CF2 = calc_MDC_mC(valid_df['mz'].values, intens_C12, intens_C13)

        # Assign results back to the original DataFrame
        df.loc[valid_mask, 'C'] = C
        df.loc[valid_mask, 'MD/C'] = MDC
        df.loc[valid_mask, 'm/C'] = mC
        df.loc[valid_mask, 'r_CF2'] = r_CF2
    
    df['MD'] = df['mz'] - np.round(df['mz'], decimals=0)

    #%%
    # KMD analysis
    # ==============================================================================================

    df_kmd = kmd_analysis(
                mz_vec = df['mz'], 
                RT_vec = df['rt'],
                diffs = diff_kmd, 
                hs_tol = mass_tolerance,
                n_min = n_homologues)

    df = pd.concat([df, df_kmd], axis = 1) 

    print(f'{len(np.unique(df[df["min_homologues"] == True]["hs_number"]))} HS detected') # move inside KMD function

    # Suspect screening
    # ==================================================================================

    df_suspect_screening = suspect_screening(df,
                                             df_suspects,
                                             tol_suspect = mass_tolerance)

    df_suspect_screening.set_index(df.index, inplace = True) # NOTE: check if this always works!
    df = pd.concat([df, df_suspect_screening], axis = 1)


    # ===================================================================================
    # Annotate MS2 spectra with subformulas
    # NOTE: bring this into a function together with match_subformulas_to_peaks

    df_annot = df[(df['mzs_ms2'].notna()) & (df['formulas'].notna())]

    results = []
    for idx in tqdm(df_annot.index, total=len(df_annot), desc="Annotating MS2 spectra with subformulas"):
        mzs_top_n, top_n_idx = get_top_n_peaks(df_annot['mzs_ms2'][idx], df_annot['ints_ms2'][idx], 10)

        matches_col = []
        for n in range(len(df_annot['formulas'][idx])):
            matches = match_subformulas_to_peaks(
                df_annot['formulas'][idx][n],
                mzs_top_n,
                mass_tolerance=0.005,
                max_combinations=100_000
            )
            matches_col.append(matches)

        # Append all results for this row
        results.append({
            'index': idx,
            'top_n_idx_ms2': list(top_n_idx),
            'annot': matches_col
        })

    # Create a new DataFrame from results
    if len(results) > 0:
        annot_df = pd.DataFrame(results).set_index('index')

        annot_df['annot'] = annot_df['annot'].astype('object')
        annot_df['top_n_idx_ms2'] = annot_df['top_n_idx_ms2'].astype('object')

        # add the explained percentage of top_n peaks
        perc_ints_explained_col = []
        for idx in annot_df.index:
            ints = np.array(df_annot['ints_ms2'][idx])
            top_n_idxs = np.array(annot_df['top_n_idx_ms2'][idx])
            if len(top_n_idxs) > 0:
                # NOTE: Check if this is correct! (second else also needs to be checked!)
                ints_total_top_n = np.sum(ints[top_n_idxs])
                perc_ints_explained = []
                for n in range(len(annot_df['annot'][idx])):
                    idx_explained = np.unique(np.array(list(annot_df['annot'][idx][n].keys())))
                    if len(idx_explained) > 0:

                        ints_explained = np.sum(ints[idx_explained])
                        perc_ints_explained.append(ints_explained/ints_total_top_n)
                else:
                    perc_ints_explained.append(0)
            else:
                perc_ints_explained = [0]
            perc_ints_explained_col.append(perc_ints_explained)

        annot_df['perc_explained'] = perc_ints_explained_col

        df = df.merge(annot_df, how='left', left_index=True, right_index=True)

    # =======================================================================================
    # Finalize DataFrame

    df['rt_min'] = df['rt'] / 60
    df = df.sort_values(by = ['m/C'], ascending = False)
    df = df.round(5)

    
    # Define priority columns to show first (left to right)
    priority_columns = ['mz','rt_min', 'adduct', 'C', 'm/C', 'MD/C', 'MD', 'n_diffs', 'n_dias', 'min_homologues', 'unique_homologues', 'intens_mean']
    
    # Get priority columns that exist in the dataframe
    existing_priority_cols = [col for col in priority_columns if col in df.columns]
    
    # Get remaining columns (not in priority list)
    remaining_cols = [col for col in df.columns if col not in priority_columns]
    
    # Reorder columns: priority first, then remaining
    df = df[existing_priority_cols + remaining_cols]
   

    # Generate HTML
    # ========================================================================================

    #plt.switch_backend('Agg')
    #generate_full_html(df_html, sample_names,
    #                os.path.join(output_folder, f'{output_name}.html'), 
    #                f'{output_name}')
    #plt.switch_backend('TkAgg')

    
    # ========================================================================================
    # Report generate files for SIRIUS

    #write_extended_mgf(df[df['mzs_ms2'].notna()], 
    #                   os.path.join(output_folder, f'{output_name}.mgf'),
    #                   polarity = polarity)

    # ========================================================================================
    # Plotting routines

    #plotting.mz_RT(df, output_folder, output_name)
    #plotting.MDC_mC_plot(df, output_folder, output_name)
    #plotting.mC_histogram(df, output_folder, output_name)
    #plotting.KMD_plot(df, output_folder, output_name, mC_limit = 0)
    
    #if save_msms_spectra == True:
    #    for idx in df[df['mz_msms'].notna()].index:
    #        plotting.MS2_spectra_plotter(df, idx, diffs, output_folder, output_name, font_size = 20)
    
    print('PFAScreen evaluation successfully finished!')
    print(f'Took: {(time.time() - start)/60:.2f} minutes!')

    return df