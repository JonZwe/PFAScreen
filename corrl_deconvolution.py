from scipy import stats
from scipy.ndimage import gaussian_filter1d
import numpy as np
import matplotlib.pyplot as plt

def corrl_deconvolution(exp, mz, rt, 
                            rt_width, 
                            int_thresh, 
                            corr_thresh,
                            smoothing = False, 
                            md_range = 1, 
                            extraction_window = 0.005,
                            ms_level = 2,
                            mass_range = None,
                            full_spectrum_range = False, 
                            plotting = False):

    '''
    Basic deconvolution of MS2 DIA (all-ions), or MS1 data based on EIC correlation, highly dependend on R2 threshold -> needs validation
    ToDo:
    - refine MD filter
    - combine multiple MS2 scans instead of only one
    - get combined intensities
    - consider changing interpolation
    - make correlation dependent on intensity
    - add spectra subtraction (before and after peak, e.g. if MS1 peak < 5% of apex intensity)
    - Optional: Remove persistent background masses from all scans
    '''

    RT_array = np.array([spec.getRT() for spec in exp])

    def EIC(exp, RT_array, mass, extraction_window, MS_level, RT_interest, RT_width):
        
        rts = []
        ints = []
        scan_idx = []
        spec_range = np.arange(np.argmin(np.abs(RT_array - (RT_interest- RT_width))), np.argmin(np.abs(RT_array - (RT_interest + RT_width)))).tolist()

        for m in spec_range:
            if exp[m].getMSLevel() == MS_level:
                _, intensities = exp[m].get_peaks()
                rts.append(exp[m].getRT())
                scan_idx.append(m)
                index_highest_peak_within_window = exp[m].findHighestInWindow(mass,extraction_window/2,extraction_window/2)
                if index_highest_peak_within_window > -1:
                    ints.append(exp[m][index_highest_peak_within_window].getIntensity())
                else:
                    ints.append(0)
        return rts, ints, scan_idx

    # find apex scan
    idx_spec_apex = int(np.argmin(np.abs(RT_array - rt)))


    if exp[idx_spec_apex].getMSLevel() == ms_level:
        idx_MS2_spec_apex = idx_spec_apex
    else: 
        if ms_level == 1:
            other_ms_level = 2
        elif ms_level == 2:
            other_ms_level = 1

        n = idx_spec_apex
        while exp[int(n)].getMSLevel() == other_ms_level:
            n += 1
        idx_MS2_spec_apex = n

    ms2_spec_apex_mz_all = exp[idx_MS2_spec_apex].get_peaks()[0]
    ms2_spec_apex_int_all = exp[idx_MS2_spec_apex].get_peaks()[1]

    ms2_spec_apex_mz_thr = ms2_spec_apex_mz_all[ms2_spec_apex_int_all > int_thresh]
    ms2_spec_apex_int_thr = ms2_spec_apex_int_all[ms2_spec_apex_int_all > int_thresh]

    if full_spectrum_range == True:
        ms2_spec_apex_mz = ms2_spec_apex_mz_thr
        ms2_spec_apex_int = ms2_spec_apex_int_thr
    else:
        ms2_spec_apex_mz = ms2_spec_apex_mz_thr[ms2_spec_apex_mz_thr < (mz + 1)]
        ms2_spec_apex_int = ms2_spec_apex_int_thr[ms2_spec_apex_mz_thr < (mz + 1)]

    if mass_range:
        idx_mass_range = (ms2_spec_apex_mz > mass_range[0]) & (ms2_spec_apex_mz < mass_range[1])
        ms2_spec_apex_mz = ms2_spec_apex_mz[idx_mass_range]
        ms2_spec_apex_int = ms2_spec_apex_int[idx_mass_range]


    # calculate the mass defect of precursor mz
    md_prec = mz - np.round(mz)
    # calculatre the mass defect of MS2 apex mzs
    ms2_spec_apex_mz_MD = ms2_spec_apex_mz - np.round(ms2_spec_apex_mz)

    indices_md_range = (ms2_spec_apex_mz_MD < (md_prec + md_range/2)) & (ms2_spec_apex_mz_MD > (md_prec - md_range/2))
    print(f'{np.sum(indices_md_range)/len(ms2_spec_apex_mz)*100:.2f} % of MS2 apex mzs remained after mass defect filtering')

    ms2_spec_apex_mz = ms2_spec_apex_mz[indices_md_range]
    ms2_spec_apex_int = ms2_spec_apex_int[indices_md_range]

    corr_coeffs = np.zeros(len(ms2_spec_apex_mz))
    rts_feature, ints_feature, scan_idx_feature = EIC(exp, RT_array, mz, extraction_window, 1, rt, rt_width)
    correl_rts = []
    correl_ints = []

    for n, frag_mz in enumerate(ms2_spec_apex_mz):

        rts, ints, scan_idx = EIC(exp, RT_array, frag_mz, extraction_window, ms_level, rt, rt_width)

        # current solution: interpolate EIC of MS1 data so that it matches the scan indices of the MS2 EIC
        # NOTE: Check if there is a better option than interpolation!
        if ms_level == 2:
            ints_feature_interpolated = np.interp(scan_idx, scan_idx_feature, ints_feature)
        elif ms_level == 1:
            ints_feature_interpolated = ints_feature	

        # if smoothing == True:
            # idea: smooth all that have their apex either in the same or direct neighbour scan
            # dist_peak_max = abs(np.argmax(ints) - np.argmax(ints_feature_interpolated))
            #   if dist_peak_max <= 1:
        #     if np.argmax(ints) == np.argmax(ints_feature_interpolated):
        #         ints = gaussian_filter1d(ints, sigma=1)
                # Considere also peaks that have also an adjacent peak apex

        # consideration: smooth EIC of MS2 prior to checking if apex it the same
        if smoothing == True:
            if max(ints)/max(ints_feature_interpolated) < 0.25:
                ints_smoothed = gaussian_filter1d(ints, sigma=1)
                if abs(np.argmax(ints_smoothed) - np.argmax(ints_feature_interpolated)):
                    ints = ints_smoothed
        
        ints = np.array(ints)
        ints_feature_interpolated = np.array(ints_feature_interpolated)

        if np.all(ints == ints[0]) or np.all(ints_feature_interpolated == ints_feature_interpolated[0]):
            # if all intensities are the same, the correlation coefficient is 0
            R2 = 0
            corr_coeffs[n] = R2
        else:
            R2 = stats.linregress(ints_feature_interpolated, ints).rvalue**2
            corr_coeffs[n] = R2
        
        if R2 > corr_thresh:
            correl_rts.append(rts)    
            correl_ints.append(ints)

    idx_corr = corr_coeffs > corr_thresh
    if np.sum(idx_corr) > 0:
        ms2_spec_apex_mz_deconv = ms2_spec_apex_mz[idx_corr]
        ms2_spec_apex_int_deconv = ms2_spec_apex_int[idx_corr]

        orig_spec = (ms2_spec_apex_mz, ms2_spec_apex_int)
        deconv_spec = (ms2_spec_apex_mz_deconv, ms2_spec_apex_int_deconv)
    else:
        print('No correlation fragment detected')
        orig_spec = (ms2_spec_apex_mz, ms2_spec_apex_int)
        deconv_spec = None

    if (plotting == True) and (orig_spec != None) and (len(ms2_spec_apex_mz_deconv) > 1):
        # =====================================
        # Plotting
        fig = plt.figure(figsize=(15,5))
        plt.subplot(1, 3, 1)
        plt.plot(rts_feature, ints_feature, linewidth = 2, color = 'blue')
        #plt.plot(rts, ints_feature_interpolated, linewidth = 1, color = 'green')
        for n in range(len(correl_rts)):
            plt.plot(correl_rts[n], correl_ints[n], color = 'red', alpha = 0.6)
        plt.xlabel('RT (s)')
        plt.ylabel('Counts ()')

        plt.subplot(1, 3, 2)
        plt.stem(ms2_spec_apex_mz, -ms2_spec_apex_int, 'Black', markerfmt=" ", basefmt=" ")
        _, stemlines, _ = plt.stem(ms2_spec_apex_mz, -ms2_spec_apex_int, 'Black',markerfmt=" ", basefmt=" ")

        if np.sum(idx_corr) > 0:
            plt.stem(ms2_spec_apex_mz_deconv, ms2_spec_apex_int_deconv, 'Red', markerfmt=" ", basefmt=" ")
            _, stemlines, _ = plt.stem(ms2_spec_apex_mz_deconv, ms2_spec_apex_int_deconv, 'Red',markerfmt=" ", basefmt=" ")
            plt.setp(stemlines, color = 'Red', linewidth = 1)

            for i, txt in enumerate(np.round(ms2_spec_apex_mz_deconv, 4)):
                plt.annotate(txt, (ms2_spec_apex_mz_deconv[i], ms2_spec_apex_int_deconv[i]), fontsize = 7, color = 'red')
            plt.xlabel('m/z')
            plt.ylabel('Counts ()')

            fig.suptitle(f'm/z = {mz} | n_peaks = {len(ms2_spec_apex_int)} -> {len(ms2_spec_apex_mz_deconv)} peaks after deconvol @ R$^2$ = {corr_thresh}')

        else:
            plt.xlabel('m/z')
            plt.ylabel('Counts ()')
            fig.suptitle('No fragments found')

        plt.subplot(1, 3, 3)
        md_ms2_spec_apex_mz = ms2_spec_apex_mz - np.round(ms2_spec_apex_mz) 
        plt.scatter(ms2_spec_apex_mz, md_ms2_spec_apex_mz, color = 'grey', alpha=0.5)
        plt.scatter(mz, (mz - np.round(mz)), 70, color = 'blue')
        if np.sum(idx_corr) > 0:
            md_ms2_spec_apex_mz_deconv = ms2_spec_apex_mz_deconv - np.round(ms2_spec_apex_mz_deconv)
            plt.scatter(ms2_spec_apex_mz_deconv, md_ms2_spec_apex_mz_deconv, color = 'red', alpha=0.8)
        plt.ylim([-0.5, 0.5])
        plt.xlabel('m/z')
        plt.ylabel('MD')
        fig.tight_layout()
        plt.show()

    return orig_spec, deconv_spec, idx_corr, corr_coeffs