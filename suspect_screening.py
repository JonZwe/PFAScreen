"""
Basic suspect screening script based on accurate mass and perform isotope pattern match
"""
import pandas as pd
import numpy as np
from utils import ismembertol
from isotope_matching import isotope_matching

def suspect_screening(df,
                      df_suspects,
                      tol_suspect):
    
    # make cleanup in suspect list (remove charges) before giving it to the functions!

    # NOTE: added as preliminary solution to folder issue
    exact_mass = df_suspects['exact_mass']

    measured_mz_neural = df['mz'].values - df['adduct_mass_diff'].values

    # NOTE: This scales really badly (N2), needs to be fixed.
    # search for matching masses with given tolerance
    _, _, idx_hit_df, idx_hit_df_susp_list = ismembertol(exact_mass, 
                                                         measured_mz_neural, 
                                                         tol_suspect)
    # consider multiple hits
    uniques = np.unique(idx_hit_df)
    names = [np.nan]*len(df)
    smiles = [np.nan]*len(df)
    masses = [np.nan]*len(df)
    formulas = [np.nan]*len(df)
    # collect information for all hits
    for unique in uniques:
        names[unique] = list(df_suspects['name'][idx_hit_df_susp_list[unique == idx_hit_df]])
        smiles[unique] = list(df_suspects['SMILES'][idx_hit_df_susp_list[unique == idx_hit_df]])
        masses[unique] = list(df_suspects['exact_mass'][idx_hit_df_susp_list[unique == idx_hit_df]])
        formulas[unique] = list(df_suspects['formula'][idx_hit_df_susp_list[unique == idx_hit_df]])

    df_suspect_hits = pd.DataFrame(data = {'compound_names':names, 
                                           'SMILES': smiles, 
                                           'exact_masses': masses, 
                                           'formulas': formulas})
    
    # perform matching between measured and theoretical isotope patters via isotope_matching function
    scores_list, mzs_theor_list, ints_theor_list, susp_idx = isotope_matching(df, 
                                                                              df_suspect_hits)

    df_isotope_score = pd.DataFrame(data = {'isotope_scores':scores_list,
                                            'mzs_isotopes_theor':mzs_theor_list,
                                            'ints_isotopes_theor':ints_theor_list})
    df_isotope_score = df_isotope_score.set_index(susp_idx)

    df_suspect_hits = pd.merge(df_suspect_hits, df_isotope_score, how = 'left', left_index=True, right_index=True)
 
    print(f'{len(uniques)} suspect hits found!')
    
    return df_suspect_hits