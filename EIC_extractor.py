# %%
# Function to build extracted ion chromatograms (EICs) based on OpenMS
import numpy as np
import matplotlib.pyplot as plt
from pylab import figure
from pyteomics import mass

def EIC_extractor(
        experiments, 
        mz_list,
        extraction_width = 0.005,
        rep_unit = 'CF2',
        n_rep_unit = 3
        ):

    mz_list_float = []
    for mz in mz_list:
        mz_list_float.append(float(mz))

    if rep_unit != '':
        n_masses = np.arange(-n_rep_unit, n_rep_unit + 1, 1)
        mz_list_float = mz_list_float[0] + mass.calculate_mass(formula = rep_unit) * n_masses
        print(mz_list_float)

    result = {'EIC':[]}
    samples = ['Sample', 'Blank']
    for n, exp in enumerate(experiments):
        for mz in mz_list_float:
            rts = []
            ints = []
            for spec in exp:
                if spec.getMSLevel() == 1:
                    _, intensities = spec.get_peaks()
                    rts.append(spec.getRT())
                    index_highest_peak_within_window = spec.findHighestInWindow(mz,extraction_width/2,extraction_width/2)
                    if index_highest_peak_within_window > -1:
                        ints.append(spec[index_highest_peak_within_window].getIntensity())
                    else:
                        ints.append(0)
            result['EIC'].append({'sample':samples[n] ,'mz': mz, 'rt': rts, 'i': ints})

    leg = []
    figure()
    for n in range(len(result['EIC'])):
        plt.plot(np.array(result['EIC'][n]['rt'])/60, result['EIC'][n]['i'], alpha = 0.8)
        leg.append(f"{result['EIC'][n]['sample']} - {np.round(result['EIC'][n]['mz'], 4)}")
    plt.legend(leg)
    plt.xlabel('RT (min)')
    plt.ylabel('Counts (-)')
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
    plt.title(f'Extraction width = {extraction_width} Da')
    plt.show()