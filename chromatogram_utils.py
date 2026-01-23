import numpy as np
import matplotlib.pyplot as plt

"""
# Collection of function associated with building extracted ion chromatograms (EICs) (for PFAScreen).
# Currently inludes:
# get_full_eic, get_eic, get_eic_start_stop, plot_full_eic, plot_eics_from_sample_info, get_all_eics, get_all_eics_2
# Jonathan Zweigle, 06/2025
"""

def get_full_tic(exp, 
                 md_range=(-0.5, 0.5), 
                 ms_level=1) -> tuple:
    """
    Function to extract a full TIC from oms.MSExperiment with option to filter by mass defect (MD).
    """
    rts = []
    ints = []
    for m in range(exp.getNrSpectra()):
        if exp[m].getMSLevel() == ms_level:
            mzs, intensities = exp[m].get_peaks()
            rts.append(exp[m].getRT())
            md = mzs - np.round(mzs)
            intensities = intensities[(md >= md_range[0]) & (md <= md_range[1])]
            ints.append(np.sum(intensities))
    
    return np.array(rts), np.array(ints)


def get_full_kmd_eic(exp, 
                     kmd_range=(-0.5, 0.5), 
                     rep_unit=49.99680, 
                     ms_level=1) -> tuple:
    """
    Function to extract a full KMD EIC from oms.MSExperiment with option to filter by KMD.
    """
    rts = []
    ints = []
    for m in range(exp.getNrSpectra()):
        if exp[m].getMSLevel() == ms_level:
            mzs, intensities = exp[m].get_peaks()
            rts.append(exp[m].getRT())
            km = mzs * np.round(rep_unit)/rep_unit
            kmd = km - np.round(km)
            intensities = intensities[(kmd >= kmd_range[0]) & (kmd <= kmd_range[1])]
            ints.append(np.sum(intensities))
    
    return np.array(rts), np.array(ints)


def get_full_eic(exp, 
                 mass, 
                 extraction_window=0.005, 
                 ms_level=1) -> tuple:
    """
    Function to extract a full EIC for a given mass and extraction window from oms.MSExperiment.
    """
    rts = []
    ints = []
    for m in range(exp.getNrSpectra()):
        if exp[m].getMSLevel() == ms_level:
            _, intensities = exp[m].get_peaks()
            rts.append(exp[m].getRT())
            index_highest_peak_within_window = exp[m].findHighestInWindow(mass,extraction_window/2,extraction_window/2)
            if index_highest_peak_within_window > -1:
                ints.append(exp[m][index_highest_peak_within_window].getIntensity())
            else:
                ints.append(0)
    
    return np.array(rts), np.array(ints)


def get_eic(exp, 
            mass, 
            rt, 
            rt_width, 
            extraction_window=0.005,
            ms_level=1) -> tuple:

    """
    Function to extract a an EIC for a given mass, rt_width and extraction window from oms.MSExperiment.
    """
    rt_array = np.array([spec.getRT() for spec in exp])
    
    rts = []
    ints = []
    scan_idx = []
    spec_range = np.arange(np.argmin(np.abs(rt_array - (rt - rt_width))), np.argmin(np.abs(rt_array - (rt + rt_width)))).tolist()

    for m in spec_range:
        if exp[m].getMSLevel() == ms_level:
            _, intensities = exp[m].get_peaks()
            rts.append(exp[m].getRT())
            scan_idx.append(m)
            index_highest_peak_within_window = exp[m].findHighestInWindow(mass,extraction_window/2,extraction_window/2)
            if index_highest_peak_within_window > -1:
                ints.append(exp[m][index_highest_peak_within_window].getIntensity())
            else:
                ints.append(0)

    rts = np.array(rts)
    ints = np.array(ints)
    scan_idx = np.array(scan_idx)

    return rts, ints, scan_idx


def get_eic_start_stop(exp, 
                       mass, 
                       rt_start, 
                       rt_stop, 
                       extraction_window=0.005, 
                       ms_level=1):

    """
    Function same as get_eic, but with start and stop RT.
    """
    rt_array = np.array([spec.getRT() for spec in exp])
    
    rts = []
    ints = []
    scan_idx = []
    spec_range = np.arange(np.argmin(np.abs(rt_array - (rt_start))), np.argmin(np.abs(rt_array - (rt_stop)))).tolist()

    for m in spec_range:
        if exp[m].getMSLevel() == ms_level:
            _, intensities = exp[m].get_peaks()
            rts.append(exp[m].getRT())
            scan_idx.append(m)
            index_highest_peak_within_window = exp[m].findHighestInWindow(mass,extraction_window/2,extraction_window/2)
            if index_highest_peak_within_window > -1:
                ints.append(exp[m][index_highest_peak_within_window].getIntensity())
            else:
                ints.append(0)

    rts = np.array(rts)
    ints = np.array(ints)
    scan_idx = np.array(scan_idx)

    return rts, ints, scan_idx


def plot_eics_from_sample_info(exps, 
                               sample_names, 
                               mass, 
                               rt, 
                               rt_width=30, 
                               extraction_window=0.005, 
                               ms_level=1) -> None:

    """
    Function to plot EICs from multiple experiments 
    based on sample names and mass.
    """   
    def get_set1_colors(n):
        cmap = plt.get_cmap('Dark2')
        return [cmap(i) for i in range(n)]

    rts_col = []
    ints_col = []
    for exp in exps:

        rt_array = np.array([spec.getRT() for spec in exp])
        
        rts = []
        ints = []
        scan_idx = []
        spec_range = np.arange(np.argmin(np.abs(rt_array - (rt - rt_width))), np.argmin(np.abs(rt_array - (rt + rt_width)))).tolist()

        for m in spec_range:
            if exp[m].getMSLevel() == ms_level:
                _, intensities = exp[m].get_peaks()
                rts.append(exp[m].getRT())
                scan_idx.append(m)
                index_highest_peak_within_window = exp[m].findHighestInWindow(mass,extraction_window/2,extraction_window/2)
                if index_highest_peak_within_window > -1:
                    ints.append(exp[m][index_highest_peak_within_window].getIntensity())
                else:
                    ints.append(0)

        rts = np.array(rts)
        ints = np.array(ints)

        rts_col.append(rts)
        ints_col.append(ints)

    cols = get_set1_colors(len(sample_names))
    #cols = ['red', 'red', 'red', 'blue', 'blue', 'blue']
    plt.figure(figsize=(6,4))
    for n in range(len(sample_names)):
        plt.plot(rts_col[n], ints_col[n], color = cols[n], alpha = 0.7)
    plt.legend(sample_names)
    plt.xlabel('RT (s)')
    plt.ylabel('Counts (-)')
    plt.title(f'm/z = {mass} @ {extraction_window} Da')


def get_all_eics(exp, 
                 fm) -> tuple:

    """
    Function to extract all EICs from a given experiment and feature map.
    Currenly df_dict is the major bottleneck for performance (still best solution so far).
    It also returns the maximum m/z deviation for each EIC, which can be useful for quality control
    and to find features which are likely within detector saturation.
    Could be improved if necessary.
    """

    df_exp = exp.get_df(long='True')
    df_dict = {(rt, mz): inty for rt, mz, inty in zip(df_exp['RT'], df_exp['mz'], df_exp['inty'])}

    rts, ints = [], []
    max_mz_dev = []
    for f in fm:

        arr = f.getConvexHulls()[0].getHullPoints()
        keys = [tuple(row) for row in arr]
        inty = [df_dict.get(key, np.nan) for key in keys]

        rts.append(arr[:, 0])
        ints.append(inty)

        max_mz_dev.append( np.max(arr[:, 1]) - np.min(arr[:, 1]) )

    return rts, ints, max_mz_dev

"""
Example usage:
rts, ints = get_all_eics(exp, fm)
plt.figure()
for rt, inty in zip(rts, ints):
    plt.plot(rt, inty, alpha=0.5)
plt.show()
"""

def get_all_eics_2(exp, 
                   fm, 
                   extraction_window=0.005) -> tuple:
    """
    Different implementation
    NOTE: Much lower performance than get_all_eics, 
    cannot capture EICs with large mass deviations (e.g., detector saturation).
    """
    df_fm = fm.get_df()

    rts, ints = [], []
    for n in range(len(df_fm)):
        rt, i, _ = get_eic_start_stop(exp,
                                   mass = df_fm['mz'].iloc[n],
                                   rt_start = df_fm['RTstart'].iloc[n],
                                   rt_stop = df_fm['RTend'].iloc[n],
                                   extraction_window = extraction_window,
                                   ms_level = 1)
        rts.append(rt)
        ints.append(i)

    return rts, ints


def get_neutral_loss_eic(exp, 
                         neutral_loss, 
                         extraction_window=0.01, 
                         ms_level=1, 
                         mass_defect_range=(-0.5, 0.5), 
                         intensity_threshold=100):
    """
    Extract an extracted ion chromatogram (EIC) for a neutral loss detected across the entire spectrum.
    This function sums the intensities of peaks corresponding to the same neutral loss event.
    It also applies a mass defect filter to remove irrelevant peaks.
    
    Parameters:
    - exp: MSExperiment object (pyopenms)
    - neutral_loss: The neutral loss mass (e.g., for H2O, NH3, etc.) in Da
    - extraction_window: The width of the extraction window around the potential neutral loss peaks (in Da)
    - MS_level: The MS level to extract data from (1 for MS1, 2 for MS2)
    - mass_defect_range: Tuple of (min_defect, max_defect) to filter peaks based on mass defect (default is (-0.2, 0.05))
    - intensity_threshold: Intensity threshold for peak inclusion (default is 0)
    
    Returns:
    - rts: Array of retention times corresponding to the detected neutral loss peaks
    - ints: Array of summed intensities for each neutral loss event
    """
    rts = []
    ints = []

    # Loop through each spectrum in the MS experiment
    for m in range(exp.getNrSpectra()):
        if exp[m].getMSLevel() == ms_level:  # Only process spectra of the specified MS level
            mz_vals, intens_vals = exp[m].get_peaks()  # Get m/z and intensities for the spectrum
            rt = exp[m].getRT()  # Retention time
            
            # Step 1: Apply intensity threshold to filter low-intensity peaks
            valid_peaks_mask = intens_vals >= intensity_threshold
            mz_vals = mz_vals[valid_peaks_mask]
            intens_vals = intens_vals[valid_peaks_mask]

            # Step 2: Calculate mass defect for all peaks
            mass_defects = mz_vals - np.round(mz_vals)

            # Apply mass defect filter (exclude peaks outside the specified range)
            valid_peaks_mask = (mass_defects >= mass_defect_range[0]) & (mass_defects <= mass_defect_range[1])
            mz_vals = mz_vals[valid_peaks_mask]
            intens_vals = intens_vals[valid_peaks_mask]

            # Step 3: Find peaks corresponding to the neutral loss (i.e., mz1 - mz2 = neutral_loss)
            summed_intensity = 0
            found_neutral_loss = False  # Flag to check if we found a neutral loss peak
            
            # Vectorized approach: Calculate all pairwise differences between mz_vals
            mz_diff_matrix = mz_vals[:, None] - mz_vals  # Broadcasted difference calculation
            intensity_matrix = intens_vals[:, None] + intens_vals  # Broadcasted intensity sum

            # Find all mass differences within the neutral loss window
            neutral_loss_mask = np.abs(mz_diff_matrix - neutral_loss) <= extraction_window
            neutral_loss_mask[np.triu(neutral_loss_mask)] = False  # Ignore duplicates (upper triangular part of the matrix)

            # Sum intensities of pairs that match the neutral loss (across the entire spectrum)
            for i in range(len(mz_vals)):
                if np.any(neutral_loss_mask[i]):
                    found_neutral_loss = True
                summed_intensity += np.sum(intensity_matrix[neutral_loss_mask[i]])

            # If no neutral loss was found, set the intensity to 0
            if not found_neutral_loss:
                summed_intensity = 0

            # Append retention time and summed intensity (or 0 if no neutral loss found)
            rts.append(rt)
            ints.append(summed_intensity)

    # Convert lists to numpy arrays for easy manipulation
    rts = np.array(rts)
    ints = np.array(ints)

    return rts, ints

#rts, ints = get_neutral_loss_eic(exp, neutral_loss = 'CF2', mass_defect_range=(-0.2, 0.05), intensity_threshold=100)


def manual_eic_integration(rt, 
                           ints):

    """
    Function to manually integrate an EIC by clicking on the plot.
    NOTE: Very preliminary, not tested yet, does not work currently.
    """
    num_clicks = 0
    click1 = None
    def on_click(event):
        global num_clicks, click1
        if event.button == 3:  # Check for right click event
            if event.inaxes:
                if num_clicks == 0:
                    click1 = (event.xdata, event.ydata)
                num_clicks += 1
                if num_clicks == 2:
                    click2 = (event.xdata, event.ydata)
                    show_data(click1, click2, rt, ints)
                    num_clicks = 0

    def show_data(click1, click2, rt, intens):
        
        x1, x2 = click1[0], click2[0]
        y1, y2 = click1[1], click2[1]

        m = (y2 - y1)/(x2 - x1)
        b = y1 - m * x1
        line = m*rt + b
        idx = np.argwhere(np.diff(np.sign(intens - line))).flatten()
        plt.plot(rt[idx], intens[idx], 'ro')
        plt.plot([x1, x2], [y1, y2], color='tomato')
        plt.fill_between(rt[idx[0]:idx[1]], intens[idx[0]:idx[1]], line[idx[0]:idx[1]], 
                         where=intens[idx[0]:idx[1]] > line[idx[0]:idx[1]], color='blue', alpha=0.3)
        area = np.trapz(intens[idx[0]:idx[1]], x=rt[idx[0]:idx[1]])
        plt.title(f'Area = {np.round(area,2)}')

    plt.figure(figsize=(8, 8))
    plt.plot(rt, ints, color = 'black')
    plt.gcf().canvas.mpl_connect('button_press_event', on_click)
    plt.show()
