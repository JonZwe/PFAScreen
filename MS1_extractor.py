# %%
import numpy as np
import matplotlib.pyplot as plt
from pylab import figure
from pyopenms import EmpiricalFormula, CoarseIsotopePatternGenerator

def MS1_extractor(
        experiment,
        RT,
        formula 
        ):
    
    RT_array = np.array([exp.getRT() for exp in experiment])
    idx = np.argmin(np.abs(RT_array - RT*60))

    if experiment[int(idx)].getMSLevel() == 2:
        n = idx
        while experiment[int(n)].getMSLevel() == 2:
            n += 1
        idx = n

    mz_array = experiment[int(idx)].get_peaks()[0]
    intens_array = experiment[int(idx)].get_peaks()[1]

    intens_array_normalized = intens_array/np.max(intens_array)
    idx_label = intens_array_normalized > 0.05
    mz_array_label = mz_array[idx_label]
    intens_array_label = intens_array[idx_label]

    figure()
    plt.stem(mz_array, intens_array, 'Black', markerfmt=" ", basefmt=" ")
    markerline, stemlines, baseline = plt.stem(mz_array, intens_array, 'Black',markerfmt=" ", basefmt=" ")
    plt.setp(stemlines, color = 'Black', linewidth= 0.5)

    for i, txt in enumerate(np.round(mz_array_label, 4)):
        plt.annotate(txt, (mz_array_label[i],intens_array_label[i]), color = 'Black', rotation = 20, fontsize=7)

    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
    plt.title(f'MS1 spectrum @ RT = {RT} min')
    plt.xlabel('m/z')
    plt.ylabel('Counts (-)')
    plt.ylim(ymin = 0)
    plt.show()


    if formula != '':
        # Theoretical isotopes
        f = EmpiricalFormula(formula)
        isotopes = f.getIsotopeDistribution( CoarseIsotopePatternGenerator(4) )
        mzs = np.array([iso.getMZ() for iso in isotopes.getContainer()])
        ints = np.array([iso.getIntensity() for iso in isotopes.getContainer()])
        ints_normalized = ints/np.max(ints)*100

        idx_cut = np.logical_and(mz_array > mzs[0]-0.5, mz_array < mzs[-1]-0.5)
        mz_array_cut = mz_array[idx_cut]
        intens_array_cut = intens_array[idx_cut]/np.max(intens_array[idx_cut])*100

        figure()
        plt.stem(mzs, ints_normalized, 'lightseagreen', markerfmt=" ", basefmt=" ")
        markerline, stemlines, baseline = plt.stem(mzs, ints_normalized, 'lightseagreen', markerfmt=" ", basefmt=" ", label='_nolegend_')
        plt.setp(stemlines, color = 'lightseagreen', linewidth= 10)

        for i, txt in enumerate(np.round(mzs, 4)):
            plt.annotate(txt, (mzs[i],ints_normalized[i]), color = 'lightseagreen', rotation = 20, fontsize=7)
        plt.stem(mz_array_cut, intens_array_cut, 'tomato', markerfmt=" ", basefmt=" ")
        markerline2, stemlines2, baseline2 = plt.stem(mz_array_cut, intens_array_cut, 'lightseagreen', markerfmt=" ", basefmt=" ", label='_nolegend_')
        plt.setp(stemlines2, color = 'tomato', linewidth= 2)

        plt.legend(['Theoretical', 'Experimental'])
        plt.title(f'Theoretical isotope match ({formula})')
        plt.xlabel('m/z')
        plt.ylabel('Relative counts (%)')
        plt.ylim(ymin = 0)
        plt.show()