import numpy as np
import pyopenms as oms
from itertools import groupby
import matplotlib.pyplot as plt


def spectra_2D_to_matrix(exp, modulation_time = 6):

    """
    Bring spectra from GCxGC run into matrix format by specifying a constant modulation time.
    """

    # create 2D spectra matrix
    rt_arr = np.array([spec.getRT() for spec in exp])

    # Compute bin indices (integer division groups numbers into bins)
    bin_indices = ((rt_arr - rt_arr[0]) // modulation_time).astype(int)

    # Get unique bin indices and their counts
    uni, counts = np.unique(bin_indices, return_counts = True)

    n_dim = counts.max()
    m_dim = uni.max()

    # Define a structured NumPy dtype with two array fields
    dt = np.dtype([('mz', 'O'), ("intensity", 'O'), ('spec_index', 'O')])

    # NOTE: here are still bugs! Spectra must be correctly placed in the matrix
    # There are still empty pixels (can be seen in if chrom_type == TIC)
    spectra_matrix = np.empty((n_dim, m_dim + 1), dtype=dt)

    row_indices = np.concatenate([np.arange(len(list(group))) for _, group in groupby(bin_indices)])

    for i in range(exp.size()):
        
        n = row_indices[i]
        m = bin_indices[i]

        spectra_matrix[n, m]["mz"] = exp[int(i)].get_peaks()[0]
        spectra_matrix[n, m]["intensity"] = exp[int(i)].get_peaks()[1]
        spectra_matrix[n, m]["spec_index"] = int(i)

    return spectra_matrix


def get_2D_chromatogram(spectra_matrix, 
                        chrom_type ='EIC', 
                        mz_interest=None, 
                        formula=None, 
                        md_range=(-0.5, 0.5),
                        rep_unit=49.99680, 
                        kmd_range=(-0.05, 0.05),
                        mz_tol=0.01, 
                        cmap='coolwarm', 
                        type='img', 
                        plotting=True):
    
    """
    Interactive 2D chromatogram (EIC or TIC), which 
    allows spectra visualization by clicking on pixels
    """

    global clicked_points
    clicked_points = []

    if formula:
        mz_interest = oms.EmpiricalFormula(formula).getMonoWeight()

    chromatogram_2D = np.zeros((spectra_matrix.shape[0], spectra_matrix.shape[1]))

    for n in range(spectra_matrix.shape[0]):
        for m in range(spectra_matrix.shape[1]):
            if chrom_type == 'TIC':
                chromatogram_2D[n, m] = np.sum(spectra_matrix[n, m]["intensity"])
            elif chrom_type == 'EIC':
                if spectra_matrix[n, m]["mz"] is None:
                    chromatogram_2D[n, m] = 0
                else:
                    idx = (spectra_matrix[n, m]["mz"] > mz_interest - mz_tol) & (spectra_matrix[n, m]["mz"] < mz_interest + mz_tol)
                    chromatogram_2D[n, m] = np.sum(spectra_matrix[n, m]["intensity"][idx]) if np.sum(idx) > 0 else 0
            elif chrom_type == 'MD':
                if spectra_matrix[n, m]["mz"] is None:
                    chromatogram_2D[n, m] = 0
                else:
                    md = spectra_matrix[n, m]["mz"] - np.round(spectra_matrix[n, m]["mz"])
                    idx = (md > md_range[0]) & (md < md_range[1])
                    chromatogram_2D[n, m] = np.sum(spectra_matrix[n, m]["intensity"][idx]) if np.sum(idx) > 0 else 0
            elif chrom_type == 'KMD':
                if spectra_matrix[n, m]["mz"] is None:
                    chromatogram_2D[n, m] = 0
                else:
                    km = spectra_matrix[n, m]["mz"] * rep_unit/np.round(rep_unit)
                    kmd = km - np.round(km)
                    idx = (kmd > kmd_range[0]) & (kmd < kmd_range[1])
                    chromatogram_2D[n, m] = np.sum(spectra_matrix[n, m]["intensity"][idx]) if np.sum(idx) > 0 else 0


    if plotting:
        #plt.ion()
        fig, ax = plt.subplots(figsize=(10, 6))

        if type == 'img':
            im = ax.imshow(chromatogram_2D, aspect='auto', cmap=cmap)
            plt.xlabel(r'$^1$D')
            plt.ylabel(r'$^2$D')
            plt.gca().invert_yaxis()
            plt.colorbar(im)
            
            def onclick(event):
                if event.xdata is not None and event.ydata is not None:
                    x, y = int(event.xdata), int(event.ydata)
                    clicked_points.append((x, y))
                    print(f"Clicked: ({x}, {y}), Value: {chromatogram_2D[y, x]}")  

                    # Plot the spectrum at the clicked point
                    mz_values = spectra_matrix[y, x]["mz"]
                    intensity_values = spectra_matrix[y, x]["intensity"]

                    if mz_values is not None and intensity_values is not None:
                        plt.figure(figsize=(8, 5))
                        _, sl1, _ = plt.stem(mz_values, intensity_values, 'Black', markerfmt=" ", basefmt=" ")
                        plt.setp(sl1, color='Black', linewidth=1)
                        
                        # Label peaks above 5% intensity
                        mz_values_label = mz_values[intensity_values / np.max(intensity_values) > 0.05]
                        intensity_values_label = intensity_values[intensity_values / np.max(intensity_values) > 0.05]
                        
                        for i, txt in enumerate(np.round(mz_values_label, 4)):
                            plt.annotate(txt, (mz_values_label[i], intensity_values_label[i]), 
                                        color='Black', rotation=0, fontsize=9)
                        plt.show(block=False)

            fig.canvas.mpl_connect('button_press_event', onclick)

        elif type == 'contour':
            x = np.arange(chromatogram_2D.shape[1])  
            y = np.arange(chromatogram_2D.shape[0])  
            X, Y = np.meshgrid(x, y)
            contour = ax.contour(X, Y, chromatogram_2D, cmap=cmap)
            plt.clabel(contour, inline=True, fontsize=8)
            plt.xlabel(r'$^1$D')
            plt.ylabel(r'$^2$D')
            plt.gca().invert_yaxis()

        elif type == '3D':
            from mpl_toolkits.mplot3d import Axes3D
            x = np.arange(chromatogram_2D.shape[1])  
            y = np.arange(chromatogram_2D.shape[0])  
            X, Y = np.meshgrid(x, y)
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X, Y, chromatogram_2D, cmap=cmap)
            plt.gca().invert_yaxis()

        plt.show()

    return chromatogram_2D
"""

while True:
    print("\nEnter new parameters for get_2D_chromatogram or type 'exit' to quit.")
    
    chrom_type = input("Chromatogram Type (EIC/TIC) [default: EIC]: ") or 'EIC'
    formula = input("Enter molecular formula (or press Enter to skip): ") or None
    mz_interest = None
    if not formula:
        mz_interest = input("Enter mz interest (or press Enter to skip): ")
        mz_interest = float(mz_interest) if mz_interest else None

    mz_tol = input("Enter mz tolerance [default: 0.01]: ") or 0.01
    mz_tol = float(mz_tol)
    
    cmap = input("Enter colormap (default: viridis): ") or 'coolwarm'
    plot_type = input("Enter plot type (img/contour/3D) [default: img]: ") or 'img'

    if chrom_type.lower() == 'exit' or cmap.lower() == 'exit' or plot_type.lower() == 'exit':
        print("Exiting...")
        break

    get_2D_chromatogram(spectra_matrix, chrom_type=chrom_type, mz_interest=mz_interest, formula=formula, mz_tol=mz_tol, cmap=cmap, type=plot_type)

"""