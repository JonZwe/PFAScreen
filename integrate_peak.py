
import numpy as np
import scipy

def integrate_peak(x, 
                   y, 
                   smoothing=True, 
                   percent_of_max=5):
    """
    Get FWHM and perform a simplified peak integration at X% of the peak maximum
    NOTE: Uncommon peak shapes lead to wrong integration, consider different integration method
    """
    # 1) Determine FWHM
    # Find the index of the peak (maximum value)
    peak_idx = np.argmax(y)
    # Get the maximum and the half maximum
    peak_max = y[peak_idx]
    peak_height = peak_max
    half_max = peak_max / 2.0
    # Find where the signal crosses the half maximum on both sides of the peak
    # Check all edge cases:

    if (len(np.where(y[peak_idx:] <= half_max)[0]) == 0) and (len(np.where(y[:peak_idx] <= half_max)[0]) > 0):
        left_idx = 0
        # right_idx = np.where(y[peak_idx:] <= half_max)[0][0] + peak_idx
        right_idx = np.argmin(np.abs(y[peak_idx:] - half_max)) + peak_idx
    elif (len(np.where(y[:peak_idx] <= half_max)[0]) == 0) and (len(np.where(y[peak_idx:] <= half_max)[0]) > 0):
        right_idx = len(y) - 1
        # left_idx = np.where(y[:peak_idx] <= half_max)[0][-1]
        left_idx = np.argmin(np.abs(y[:peak_idx] - half_max))
    elif (len(np.where(y[:peak_idx] <= half_max)[0]) == 0) and (len(np.where(y[peak_idx:] <= half_max)[0]) == 0):
        left_idx = 0
        right_idx = len(y) - 1
    else:
        # Normal case where peak is in the middle of the EIC!
        # left_idx = np.where(y[:peak_idx] <= half_max)[0][-1]
        # right_idx = np.where(y[peak_idx:] <= half_max)[0][0] + peak_idx

        left_idx = np.argmin(np.abs(y[:peak_idx] - half_max))
        right_idx = np.argmin(np.abs(y[peak_idx:] - half_max)) + peak_idx
        # FWHM is the distance between the two crossing points
    fwhm_value = x[right_idx] - x[left_idx]

    # ====================================================
    # 2) Integrate peak
    percent_of_max_int = peak_max*percent_of_max/100

    if (len(np.where(y[peak_idx:] <= percent_of_max_int)[0]) == 0) and (len(np.where(y[:peak_idx] <= percent_of_max_int)[0]) > 0):
        left_idx = 0
        right_idx = np.argmin(np.abs(y[peak_idx:] - percent_of_max_int)) + peak_idx
    elif (len(np.where(y[:peak_idx] <= percent_of_max_int)[0]) == 0) and (len(np.where(y[peak_idx:] <= percent_of_max_int)[0]) > 0):
        right_idx = len(y) - 1
        left_idx = np.argmin(np.abs(y[:peak_idx] - percent_of_max_int))
    elif (len(np.where(y[:peak_idx] <= percent_of_max_int)[0]) == 0) and (len(np.where(y[peak_idx:] <= percent_of_max_int)[0]) == 0):
        left_idx = 0
        right_idx = len(y) - 1
    else:
        # Normal case where peak is in the middle of the EIC!
        left_idx = np.argmin(np.abs(y[:peak_idx] - percent_of_max_int))
        right_idx = np.argmin(np.abs(y[peak_idx:] - percent_of_max_int)) + peak_idx
    
    y_integration = y[left_idx:right_idx]
    if smoothing == True:
        y_integration = scipy.ndimage.gaussian_filter1d(y_integration, sigma=1)

    peak_area = np.trapz(y_integration)

    return peak_area, peak_height, left_idx, right_idx, fwhm_value