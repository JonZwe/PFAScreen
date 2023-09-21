# %%
# Basic function to perform a blank correction:
# Input data: np.array of mass, RT, intensity for both sample and blank
# m_tol: absolute mass tolerance, RT_tol: retention time tolerance, fold_change: desired fold change
# Output data: boolean indices of overlapping features (not in blank and inverted: in blank)
import numpy as np

def blank_correc(mass_vec, 
                 RT_vec, 
                 intens_vec, 
                 mass_vec_B, 
                 RT_vec_B, 
                 intens_vec_B, 
                 m_tol, RT_tol, 
                 fold_change
                 ):
    
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