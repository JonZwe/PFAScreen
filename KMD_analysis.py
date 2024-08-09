import pandas as pd
import numpy as np
from pyteomics import mass
import re

def KMD_analysis(
        mz_vec, 
        RT_vec, 
        diffs = 'CF2', 
        hs_tol = 0.005, 
        n_min = 3
        ):
    
    '''
    Function to perform Kendrick mass defect (KMD) analysis based on accurate mass and a predefined repeating unit
    Parameters: mz_vec: array of m/z values of a MS1 feature list, RT_vec: corresponding retention time array
    diffs: string of repeating unit (e.g., CF2), hs_tol: absolute mass tolerance, n_min: minimum number of homologues
    to be clustered into a homologues series
    '''

    # Convert formula to masses
    m_diffs = np.array(np.zeros((len(diffs))))
    for n in range(len(diffs)):
        if re.search('[a-zA-Z]', diffs[n]) == None:
            rep_unit = float(diffs[n])
        else:
            rep_unit = mass.calculate_mass(formula = diffs[n])
            
    # =========================================================================
    # Create homologous series (HS) based on the user input
    # =========================================================================
    
    # Calculate modulo: Compounds from the same homologous series bear an identical modulo
    modulo = mz_vec % rep_unit

    # Calculation of Kendrick masses and features that are in the same homologous series
    KM = mz_vec*round(rep_unit)/rep_unit
    KM_round = np.round_(KM)
    KMD = KM - KM_round
    Mod_HS = {'mz':mz_vec,'rt':RT_vec,'mod':modulo,'KMD':KMD}
    Mod_HS_Dataframe = pd.DataFrame(data=Mod_HS)
    Mod_sorted_Df_temp=Mod_HS_Dataframe.sort_values('mod') # Sort according to ascending modulo
    Mod_sorted_Df = Mod_sorted_Df_temp.reset_index(drop=False)

    # Assign identical HS number to compounds with identical modulo
    HS_num = np.zeros(len(mz_vec))
    exit_loop = 0
    for n in range(len(mz_vec)-1):
        i = n
        if HS_num[n] == 0 and n != 0:
            HS_num[n] = n
            while Mod_sorted_Df['mod'][i+1]-Mod_sorted_Df['mod'][n] < hs_tol and exit_loop == 0:
                HS_num[i+1] = n
                HS_num[n] = n
                if i < (len(mz_vec)-2):
                    i = i+1
                else:
                    exit_loop = 1
    if HS_num[-1] == 0:
        HS_num[-1] = n+1

    Mod_sorted_Df['hs_number'] = HS_num

    # Calculate number of members in HS
    HS_num_temp = np.unique(HS_num, return_counts = True)
    hsnumber = []
    for n in range(len(HS_num_temp[1])):
        for s in range(HS_num_temp[1][n]):
            hsnumber.append(HS_num_temp[1][n])

    Mod_sorted_Df['homologues'] = hsnumber
    
    # Calculate corrected number of members in HS (same masses are neglected)
    HS_num2 = np.zeros(len(mz_vec))
    for n in range(len(HS_num2)):
        HS_num2[n] = len(np.unique(round(Mod_sorted_Df['mz'][Mod_sorted_Df['hs_number'] == Mod_sorted_Df['hs_number'][n]])))
        
    Mod_sorted_Df['unique_homologues'] = HS_num2

    # Check if number of homologues is greater than specified minimum value
    HS_min_check = []
    for n in range(len(hsnumber)):
        if Mod_sorted_Df['unique_homologues'][n] >= n_min:
            HS_min_check.append(True)
        else:
            HS_min_check.append(False)

    Mod_sorted_Df['min_homologues'] = HS_min_check

    # Check the features that were also detected by their fragments
    Mod_sorted_Df.sort_index(inplace = True)

    # Sort according to ascending modulo
    Mod_sorted_Df.sort_values('min_homologues', inplace=True)

    # Sort according to index
    Mod_sorted_Df.set_index('index', inplace=True)

    return Mod_sorted_Df[['KMD', 'hs_number', 'unique_homologues', 'min_homologues']]
