import pandas as pd
import numpy as np
import pyopenms as oms

def kmd_analysis(
        mz_vec, 
        RT_vec, 
        diffs = 'CF2', 
        hs_tol = 0.005, 
        n_min = 3
        ):
    '''
    Optimized function to perform Kendrick mass defect (KMD) analysis
    Parameters: mz_vec: array of m/z values of a MS1 feature list, RT_vec: corresponding retention time array
    diffs: string of repeating unit (e.g., CF2), hs_tol: absolute mass tolerance, n_min: minimum number of homologues
    to be clustered into a homologues series
    '''

    # Convert formula to masses
    if isinstance(diffs, str):
        rep_unit = oms.EmpiricalFormula(diffs).getMonoWeight()
    else:
        rep_unit = diffs
            
    # =========================================================================
    # Create homologous series (HS) based on the user input
    # =========================================================================
    
    # Calculate modulo: Compounds from the same homologous series bear an identical modulo
    modulo = mz_vec % rep_unit

    # Calculation of Kendrick masses
    KM = mz_vec * round(rep_unit) / rep_unit
    KM_round = np.round(KM, decimals=0)
    KMD = KM - KM_round
    
    # Create DataFrame
    Mod_HS_Dataframe = pd.DataFrame({
        'mz': mz_vec,
        'rt': RT_vec,
        'mod': modulo,
        'KMD': KMD
    })
    Mod_sorted_Df = Mod_HS_Dataframe.sort_values('mod').reset_index(drop=False)
    
    # =========================================================================
    # OPTIMIZED: Assign HS numbers - preserves original logic
    # =========================================================================
    
    HS_num = np.zeros(len(mz_vec))
    mod_values = Mod_sorted_Df['mod'].values
    
    # Vectorized approach: find groups where consecutive differences < tolerance
    for n in range(len(mz_vec) - 1):
        if HS_num[n] == 0 and n != 0:  # Keep original condition
            HS_num[n] = n
            # Find consecutive entries within tolerance
            i = n
            while i < len(mz_vec) - 1 and mod_values[i + 1] - mod_values[n] < hs_tol:
                HS_num[i + 1] = n
                i += 1
    
    if HS_num[-1] == 0:
        HS_num[-1] = len(mz_vec) - 1
    
    Mod_sorted_Df['hs_number'] = HS_num
    
    # =========================================================================
    # OPTIMIZED: Calculate homologues count using groupby
    # =========================================================================
    
    # Count members per HS
    hs_counts = Mod_sorted_Df.groupby('hs_number').size()
    Mod_sorted_Df['homologues'] = Mod_sorted_Df['hs_number'].map(hs_counts)
    
    # =========================================================================
    # OPTIMIZED: Calculate unique homologues using groupby
    # =========================================================================
    
    unique_counts = Mod_sorted_Df.groupby('hs_number')['mz'].transform(
        lambda x: len(np.unique(np.round(x, decimals=5)))
    )
    Mod_sorted_Df['unique_homologues'] = unique_counts
    
    # =========================================================================
    # OPTIMIZED: Vectorized min_homologues check
    # =========================================================================
    
    Mod_sorted_Df['min_homologues'] = Mod_sorted_Df['unique_homologues'] >= n_min
    
    # Sort back to original order
    Mod_sorted_Df = Mod_sorted_Df.sort_index().set_index('index')
    
    return Mod_sorted_Df[['KMD', 'hs_number', 'unique_homologues', 'min_homologues']]