# Function to perform matching of theoretical with measured isotope patterns
# Currently only M, M+1 and M+2 isotopes are considered
from pyopenms import EmpiricalFormula, CoarseIsotopePatternGenerator
import numpy as np

def isotope_matching(Df_FeatureData, 
                     Df_suspect_hits):

    # score similar to the isotope score in MZMine https://doi.org/10.1021/ac3000418
    def isotope_score(intens_array_theoretical, intens_array_measured):

        intens_difference = abs(intens_array_theoretical - intens_array_measured)
        score = 1
        for n in range(len(intens_difference)):
            score *= (1 - intens_difference[n])
            
        return score

    # NOTE: Add relative/absolute error, error propagation, ppm mass error, n_isotopes in the future?
    scores_list = []
    # relative_deviations_list = []
    susp_idx = np.where(Df_suspect_hits['formulas'].notna())[0]
    for idx in susp_idx:

        feature = Df_FeatureData.iloc[idx]

        mz_measured = feature['mz']
        intens_array_measured = np.array([feature['mz_area'], feature['mz+1_area'], feature['mz+2_area']])
        intens_array_measured = intens_array_measured[~np.isnan(intens_array_measured)] # remove isotopes if not detected
        intens_array_measured_normalized = intens_array_measured/np.max(intens_array_measured)

        scores = []
        relative_deviations = []
        for formula_str in Df_suspect_hits['formulas'][idx]:

            if isinstance(formula_str, str):
                formula_oms = EmpiricalFormula(formula_str)
                isotopes = formula_oms.getIsotopeDistribution(CoarseIsotopePatternGenerator(len(intens_array_measured_normalized)))
                mzs = np.array([iso.getMZ() for iso in isotopes.getContainer()])
                intens_theoretical = np.array([iso.getIntensity() for iso in isotopes.getContainer()])
                intens_theoretical_normalized = intens_theoretical/np.max(intens_theoretical)

                score = isotope_score(intens_theoretical_normalized, intens_array_measured_normalized)

                # relative_deviation = abs( ( (intens_array_measured_normalized - intens_theoretical_normalized) / intens_theoretical_normalized ) )

                scores.append(score)
                # relative_deviations.append(relative_deviation)
            else:
                scores.append(np.nan)
                # relative_deviations.append(np.nan)

        scores_list.append(scores)
        # relative_deviations_list.append(relative_deviations)

    return scores_list, susp_idx

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
