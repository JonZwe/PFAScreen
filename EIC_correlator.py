import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from pylab import figure
from scipy import stats
from textwrap import wrap
import pyteomics

def EIC_correlator(
        exp,
        Df_FeatureData,
        mz_interest,
        extraction_width,
        RT_interest,
        RT_width,
        R2_threshold
        ):

    Df = Df_FeatureData.copy()
    RT_interest = RT_interest*60

    cpd_range = np.logical_and(Df['RT'] > RT_interest - RT_width, Df['RT'] < RT_interest + RT_width)
    masses = Df['m/z'][cpd_range].values
    RTs = Df['RT'][cpd_range].values
    ints = Df['m/z intens'][cpd_range].values

    idx_mz_interest = np.where(np.logical_and(np.abs(masses - mz_interest) < 0.005, 
                                              np.abs(RTs - RT_interest) < 10))[0]
    if len(idx_mz_interest) > 1:
        idx_mz_interest = idx_mz_interest[np.argmax(ints[idx_mz_interest])]
    else:
        idx_mz_interest = idx_mz_interest[0]

    RT_array = np.array([spec.getRT() for spec in exp])

    def EIC(exp, 
            masses, 
            extraction_window, 
            RT_interest, 
            RT_array, 
            RT_width):
        
        all_RT_arrays = [None]*len(masses)
        all_ints_arrays = [None]*len(masses)
        for n, mass in enumerate(masses):
            rts = []
            ints = []
            for m in np.arange(np.argmin(np.abs(RT_array - (RT_interest- RT_width))), np.argmin(np.abs(RT_array - (RT_interest + RT_width)))).tolist():
                if exp[m].getMSLevel() == 1:
                    _, intensities = exp[m].get_peaks()
                    rts.append(exp[m].getRT())
                    index_highest_peak_within_window = exp[m].findHighestInWindow(mass,extraction_window/2,extraction_window/2)
                    if index_highest_peak_within_window > -1:
                        ints.append(exp[m][index_highest_peak_within_window].getIntensity())
                    else:
                        ints.append(0)
            all_RT_arrays[n] = rts
            all_ints_arrays[n] = ints
        return all_RT_arrays, all_ints_arrays

    RT_arrays, ints_arrays = EIC(exp, 
                                 masses, 
                                 extraction_width, 
                                 RT_interest, 
                                 RT_array, 
                                 RT_width)

    R2 = np.zeros(len(RT_arrays))
    for n in range(len(RT_arrays)):
        R2[n] = stats.linregress(np.array(ints_arrays[idx_mz_interest]), np.array(ints_arrays[n])).rvalue**2

    # Plotting
    # ====================================================================================

    fig, axs = plt.subplots(1, 4, figsize=(15, 5))
    plt.subplot(1, 4, 1)
    plt.plot(RT_arrays[idx_mz_interest], ints_arrays[idx_mz_interest], '-', color = 'red')
    plt.xlabel('RT')
    plt.ylabel('Counts')

    plt.subplot(1, 4, 2)
    for n in range(len(RT_arrays)):
        if R2[n] > R2_threshold:
            plt.plot(np.array(RT_arrays[n]),
                     np.array(ints_arrays[n]), '-',alpha = 0.6)
    plt.xlabel('RT')
    plt.ylabel('Counts')

    plt.subplot(1, 4, 3)
    titl = []
    idx_infrags = []
    for n in range(len(RT_arrays)):
        if R2[n] > R2_threshold:
            plt.scatter(ints_arrays[idx_mz_interest], ints_arrays[n], alpha = 0.6)
            titl.append(masses[n])
            idx_infrags.append(n)   
    plt.xlabel('Counts (m/z of interest)')
    plt.ylabel('Counts (all EICs)')

    plt.subplot(1, 4, 4)
    for n in range(len(RT_arrays)):
        if R2[n] > R2_threshold:
            plt.scatter(np.array(ints_arrays[idx_mz_interest])/np.max(np.array(ints_arrays[idx_mz_interest])),
                        np.array(ints_arrays[n])/np.max(np.array(ints_arrays[n])), alpha = 0.6)
    plt.xlabel('Normalized counts (m/z of interest)')
    plt.ylabel('Normalized counts (all EICs)')

    titl_summary = f'{np.round(masses[idx_mz_interest], 4)}, {np.sum(R2 > 0.9)} EICs of {len(RT_arrays)} correlate with R2 > {R2_threshold} | @ RT = {np.round(RT_interest/60, 2)}, RT_width = {np.round(RT_width/60, 2)}, Extr_width = {extraction_width}'
    titl_str = "\n".join(wrap(str(titl), 60))

    plt.suptitle(titl_summary + '\n' + titl_str, fontsize = 9)
    fig.tight_layout()
    plt.show()
    print(titl)

    # == FindPFAS ================================================================
    # match specta

    mz_tol = 0.005
    mz_spec = masses[idx_infrags]
    intens_spec = ints[idx_infrags]

    diffs =  ['CH3COO', 'CF2', 'C2F4', 'C3F6', 'HF', 'H2O', 'CF2O', 'C6H3F9', 'C8H3F13', 'C10H3F17', 'C12H3F21', 'BrH', 'ClH', 'H3F3']

    m_diffs = np.array(np.zeros((len(diffs))))
    for n in range(len(diffs)):
        m_diffs[n] = pyteomics.mass.calculate_mass(formula = diffs[n])

    XX, YY = np.meshgrid(mz_spec, mz_spec)
    diff_matrix = abs(XX-YY)

    l_bound = m_diffs - mz_tol
    u_bound = m_diffs + mz_tol

    # create empty list for logicals
    dist_logical = [None] *len(m_diffs)
    idx_differences = [None] *len(m_diffs)
    for n in range(len(m_diffs)):
        # find desired differences in difference matrices 
        dist_logical[n] = np.logical_and(diff_matrix >= l_bound[n], diff_matrix <= u_bound[n]).astype(int)
        # create empty index arrays and fill with indices 
        idx_differences[n] = np.array(np.zeros((len(mz_spec))))
        idx_differences[n] = np.array(np.sum(dist_logical[n], axis = 0))
        idx_differences[n] = idx_differences[n] > 0

    spec_mz_diffs = [None]*len(m_diffs)
    spec_intens_diffs = [None]*len(m_diffs)
    for n in range(len(m_diffs)):
        spec_mz_diffs[n] =  np.array(np.zeros((sum(idx_differences[n]))))
        spec_intens_diffs[n] = np.array(np.zeros((sum(idx_differences[n]))))
        spec_mz_diffs[n] = mz_spec[idx_differences[n]]
        spec_intens_diffs[n] = intens_spec[idx_differences[n]]

    # ====================================================================================

    figure()
    _, sl, _ = plt.stem(masses[idx_infrags], ints[idx_infrags], markerfmt=" ", basefmt=" ", label='_nolegend_')
    plt.setp(sl, color = 'Darkgrey', linewidth = 5)
    for i, txt in enumerate(np.round(masses[idx_infrags], 4)):
        plt.annotate(txt, (masses[idx_infrags][i], ints[idx_infrags][i]),
                            rotation = 85, fontsize = 9)
    leg_diff = []
    col = plt.cm.rainbow(np.linspace(0, 1, len(m_diffs)))
    linew = 4
    for n in range(len(m_diffs)):
        if len(spec_mz_diffs[n]) > 0:
            _, sl, _ = plt.stem(spec_mz_diffs[n], spec_intens_diffs[n], markerfmt=" ", basefmt=" ")
            linew = linew - 0.7
            plt.setp(sl, color=col[n], linewidth = linew)
            leg_diff.append(diffs[n])

    plt.title(f' In-source fragementation spectrum of m/z = {mz_interest}')
    plt.xlabel('m/z (Insource)')
    plt.ylabel('Peak area (MS1)')
    plt.ylim(ymin = 0)
    plt.legend(leg_diff)
    plt.show()