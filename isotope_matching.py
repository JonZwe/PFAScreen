# Function to perform matching of theoretical with measured isotope patterns
from pyopenms import EmpiricalFormula, CoarseIsotopePatternGenerator
import numpy as np
import pandas as pd

def isotope_score(intens_array_theoretical, 
                  intens_array_measured):

    # Isotope score similar to Pluskal et al. 2012 (10.1021/ac3000418)
    # which is/was used in MZMine
    intens_difference = abs(intens_array_theoretical - intens_array_measured)
    score = 1
    for n in range(len(intens_difference)):
        score *= (1 - intens_difference[n])
        
    return score

def isotope_matching(df, 
                     df_suspect_hits):

    # NOTE: Add relative/absolute error, error propagation, ppm mass error, n_isotopes in the future?
    scores_list = []
    mzs_theor_list = []
    ints_theor_list = []
    susp_idx = np.where(df_suspect_hits['formulas'].notna())[0]

    susp_idx = np.where(
        df_suspect_hits['formulas'].notna() &
        df['ints_isotopes'].notna()
    )[0]

    for idx in susp_idx:

        feature = df.iloc[idx]

        ints_arr_exp = np.array(feature['ints_isotopes'])
        ints_arr_exp_norm = ints_arr_exp/np.max(ints_arr_exp)

        scores = []
        mzs_theor = []
        ints_theor = []
        for formula_str in df_suspect_hits['formulas'][idx]:

            if isinstance(formula_str, str):
                
                # remove charges
                # NOTE: BUG: Does not work in case of muliple charges CH3+3 -> CH33!!! FIX! 
                formula_str = formula_str.replace("+", "").replace("-", "")


                formula_oms = EmpiricalFormula(formula_str)
                isotopes = formula_oms.getIsotopeDistribution(CoarseIsotopePatternGenerator(len(ints_arr_exp_norm)))
                mzs_arr_theor = np.array([iso.getMZ() for iso in isotopes.getContainer()])
                ints_arr_theor = np.array([iso.getIntensity() for iso in isotopes.getContainer()])
                ints_arr_theor_norm = ints_arr_theor/np.max(ints_arr_theor)

                # NOTE: edgecase: if there are not enough isotopes in the given formula, append zeros
                if len(ints_arr_theor_norm) < len(ints_arr_exp_norm):
                    diff = len(ints_arr_exp_norm) - len(ints_arr_theor_norm)
                    ints_arr_theor_norm = np.pad(ints_arr_theor_norm, (0, diff), 'constant')

                score = isotope_score(ints_arr_theor_norm, ints_arr_exp_norm)

                # relative_deviation = abs( ( (intens_array_measured_normalized - intens_theoretical_normalized) / intens_theoretical_normalized ) )
                scores.append(score)
                mzs_theor.append(mzs_arr_theor)
                ints_theor.append(ints_arr_theor_norm)

            else:
                scores.append(np.nan)
                mzs_theor.append(np.nan)
                ints_theor.append(np.nan)

        scores_list.append(scores)
        mzs_theor_list.append(mzs_theor)
        ints_theor_list.append(ints_theor)

    return scores_list, mzs_theor_list, ints_theor_list, susp_idx

'''
import itertools
import matplotlib.pyplot as plt
from pylab import figure
Scor = list(itertools.chain.from_iterable(scores_list))
figure()
plt.hist(Scor, bins = 100)
plt.xlabel('Score')
plt.ylabel('Counts')
S = np.array(Scor)
S = np.unique(S)
print(np.sum(S<=0.95)/len(S))
'''