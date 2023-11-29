"""
Basic suspect screening script based on accurate mass
"""
import pandas as pd
import numpy as np
from ismembertol import ismembertol
from isotope_matching import isotope_matching

def suspect_screening(
        tol_suspect,
        adducts,
        measured_mz,
        Df_FeatureData
        ):
    
    Df_susp_list = pd.read_csv('suspect_list.csv') # read in .csv suspect list
    exact_mass = Df_susp_list['exact_mass']

    # search for matching masses with given tolerance
    if adducts == 1: # 1: M-H, 2: M+H, 3: M+
        _, _, Mbool = ismembertol(exact_mass, measured_mz + 1.0072, tol_suspect)
    elif adducts == 2:
        _, _, Mbool = ismembertol(exact_mass, measured_mz - 1.0072, tol_suspect)
    elif adducts == 3:
        _, _, Mbool = ismembertol(exact_mass, measured_mz, tol_suspect)

    idx_hit_Df_FeatureData = np.where(Mbool)[0]
    idx_hit_Df_susp_list = np.where(Mbool)[1]

    # consider multiple hits
    uniques = np.unique(idx_hit_Df_FeatureData)
    Names = [np.nan]*len(Df_FeatureData)
    Smiles = [np.nan]*len(Df_FeatureData)
    Masses = [np.nan]*len(Df_FeatureData)
    Formulas = [np.nan]*len(Df_FeatureData)
    # collect information for all hits
    for unique in uniques:
        Names[unique] = list(Df_susp_list['name'][idx_hit_Df_susp_list[unique == idx_hit_Df_FeatureData]])
        Smiles[unique] = list(Df_susp_list['SMILES'][idx_hit_Df_susp_list[unique == idx_hit_Df_FeatureData]])
        Masses[unique] = list(Df_susp_list['exact_mass'][idx_hit_Df_susp_list[unique == idx_hit_Df_FeatureData]])
        Formulas[unique] = list(Df_susp_list['formula'][idx_hit_Df_susp_list[unique == idx_hit_Df_FeatureData]])

    Df_suspect_hits = pd.DataFrame(data = {'compound_names':Names, 
                                           'SMILES': Smiles, 
                                           'exact_masses': Masses, 
                                           'formulas': Formulas})
    
    # perform matching between measured and theoretical isotope patters via isotope_matching function
    scores_list, susp_idx = isotope_matching(Df_FeatureData, Df_suspect_hits)

    # write isotope matching data at respective indices of Df_suspect_hits
    Df_suspect_hits['isotope_scores'] = np.nan
    Df_suspect_hits['relative_isotope_intensity_deviation'] = np.nan

    scores_list_object = np.array(scores_list, dtype="object")
    # relative_deviations_list_object = np.array(relative_deviations_list, dtype="object")

    Df_suspect_hits.loc[susp_idx, 'isotope_scores'] = scores_list_object
    # Df_suspect_hits.loc[susp_idx, 'relative_isotope_intensity_deviation'] = relative_deviations_list_object
    # NOTE: Potential BUG in relative_deviations_list_object: if all nested lists have the same dimensions a 3D array is generated! 
    # Needs to be fixed if relative_deviations are used in the future!
 
    print(f'{len(uniques)} suspect hits found')
    
    return Df_suspect_hits