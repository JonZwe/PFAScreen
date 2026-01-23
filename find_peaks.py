import scipy
import numpy as np
import scipy.signal
from chromatogram_utils import get_full_eic

def find_peaks(exp, 
               mz, 
               extraction_window=0.005, 
               smoothing=False, 
               prominence=1000, 
               height=None, 
               baseline_percent=90, 
               ms_level=1,
               print_warnings=True):

    """
    Very basic peak detection (scipy.signal.find_peaks) based on local maxima.
    If highly noisy chromatograms, too much peaks are detected
    NOTE: Consider improving robustness
    """
    
    rts, ints = get_full_eic(exp, mz, extraction_window = extraction_window, ms_level = ms_level)
    ints_orig = ints

    if smoothing == True:
        ints = scipy.ndimage.gaussian_filter1d(ints, sigma=1)
        ints_smoothed = ints
    else:
        ints_smoothed = None

    peak_indices, _ = scipy.signal.find_peaks(ints, 
                                              prominence=prominence, 
                                              height = height) 
    peak_width = scipy.signal.peak_widths(ints, peak_indices, rel_height=1)


    def estimate_baseline_from_histogram(intensities, lower_percent=90):
        """
        Estimate baseline intensity by taking the average of the lower `lower_percent` 
        of intensity values based on a histogram approach.
        """
        # sort them
        sorted_intensities = np.sort(intensities)
        
        # Determine the number of values to include in the baseline calculation
        cutoff_index = int(len(sorted_intensities) * (lower_percent / 100))
        
        # Take the average of the lower `lower_percent` intensity values
        baseline_intensity = np.mean(sorted_intensities[:cutoff_index])
        #baseline_std = np.std(sorted_intensities[:cutoff_index])
        
        return baseline_intensity

    baseline = estimate_baseline_from_histogram(ints, lower_percent=baseline_percent)

    if print_warnings == True:
        if len(peak_indices) > 1:
            print(f'Warning: {len(peak_indices)} peaks in EIC detected!')
        elif peak_indices.size == 0:
            print('Warning!: No peaks in EIC detected!')

        # NOTE: Could be removed in the future: for testing purposes
        #plt.figure()
        #plt.plot(rts, ints, color = 'green', linewidth = 2)
        # NOTE: plot_dir is not a variable in this function
        #plt.savefig(os.path.join(plot_dir, f'{np.round(mz,4)}_not_detected.png'), dpi = 50)
        #plt.close()

    return rts, ints_orig, ints_smoothed, peak_indices, peak_width, baseline