# Functions with plotting routines (e.g., spectra, isotope mirror, KMD, MD/C-m/C, m/z vs. RT) for PFAScreen
import os
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import altair as alt
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
#from IPython.display import display
import pandas as pd
from kmd_analysis import kmd_analysis
import mplcursors
import pyperclip
from find_peaks import find_peaks
from chromatogram_utils import get_full_eic

"""
Collection of various visualization functions (for PFAScreen)
Currently including:
spectrum_plot, smiles_subplot, save_smiles_to_svg, isotope_mirror_plot, ms2_spectrum_plot, isotope_patters_raw_plot
mz_rt_plot_interactive, mdc_mc_plot_interactive, mc_histogram_plot_interactive, kmd_plot_interactive
kmd_plot_interactive_BACKUP, MS2_spectra_plotter, mz_RT_MSMS, mass_spec_3d_html, feature_map_plotter
show_data_reduction
Jonathan Zweigle, 06/2025
"""

def spectrum_plot(mzs_arr, 
                  ints_arr, 
                  top_n_label = 10, 
                  matches = None, 
                  top_n_idx = None, 
                  title=None):

    mzs_arr = np.array(mzs_arr)
    ints_arr = np.array(ints_arr)

    fig = plt.figure(figsize=(6,4))

    plt.stem(mzs_arr, ints_arr, 'Black', markerfmt=" ", basefmt=" ")
    _, stemlines, _ = plt.stem(mzs_arr, ints_arr, 'Black',markerfmt=" ", basefmt=" ")
    plt.setp(stemlines, color = 'Black', linewidth= 0.5)

    if len(mzs_arr) < top_n_label:
        top_n_idx = np.arange(len(mzs_arr))
    else:
        top_n_idx = np.argpartition(ints_arr, -top_n_label)[-top_n_label:]
        # Sort these indices by intensity for nicer annotation order (optional)
        top_n_idx = top_n_idx[np.argsort(-ints_arr[top_n_idx])]
        
    #idx_percent = ints_arr/np.max(ints_arr) > percent_label
    for i, txt in enumerate(np.round(mzs_arr[top_n_idx], 4)):
        plt.annotate(txt, (mzs_arr[top_n_idx][i], ints_arr[top_n_idx][i]), color = 'Black', rotation = 20, fontsize=7)

    if isinstance(matches, dict):

        for k, formulas in matches.items():

            best = sorted(formulas, key=lambda x: x['delta'])[0]
            orig_idx = top_n_idx[k]
            mz = mzs_arr[orig_idx]
            intensity = ints_arr[orig_idx]
            
            offset_formula = np.max(ints_arr) * 0.05
            offset_delta = np.max(ints_arr) * 0.10

            plt.annotate(best['formula'],
                        (mz, intensity + offset_formula),
                        color='blue', rotation=20, fontsize=7)

            plt.annotate(f"Δ={best['delta']:.4f}",
                        (mz, intensity + offset_delta),
                        color='blue', rotation=20, fontsize=7)

    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)

    if title:
        plt.title(title)
    plt.xlabel('m/z')
    plt.ylabel('Counts (-)')
    plt.ylim(bottom=0)
    fig.tight_layout()
    plt.show()

    return fig

def feature_overview_plot(df, 
                          feature_index, 
                          sample_names=None):
    """
    IN DEVELOPMENT!
    """
    feature = df.loc[feature_index]

    plt.figure(figsize=(10,8))

    plt.subplot(2,3,1)
    if isinstance(sample_names, pd.core.series.Series) or isinstance(sample_names, list):
        # Simple color cycling through matplotlib's default color cycle
        colors = plt.cm.tab10(np.arange(len(sample_names)) % 10)
        
        values = feature[sample_names].values
        plt.bar(range(len(sample_names)), values, color=colors)
        plt.xticks(range(len(sample_names)), sample_names, rotation=90, fontsize=7)
        plt.ylabel('Intensity', fontsize=12)
        plt.tight_layout()
    else:
        plt.text(0.5, 0.5, 'No sample names provided', fontsize=14, ha='center', va='center')
        plt.axis('off')

    plt.subplot(2,3,2)
    top_n_label = 20
    if isinstance(feature['mzs_ms2'], list):
        mzs_arr = np.array(feature['mzs_ms2'])
        ints_arr = np.array(feature['ints_ms2'])
        plt.stem(mzs_arr, ints_arr, 'Black', markerfmt=" ", basefmt=" ")
        _, stemlines, _ = plt.stem(mzs_arr, ints_arr, 'Black',markerfmt=" ", basefmt=" ")
        plt.setp(stemlines, color = 'Black', linewidth= 0.5)

        if len(mzs_arr) < top_n_label:
            top_n_idx = np.arange(len(mzs_arr))
        else:
            top_n_idx = np.argpartition(ints_arr, -top_n_label)[-top_n_label:]
            # Sort these indices by intensity for nicer annotation order (optional)
            top_n_idx = top_n_idx[np.argsort(-ints_arr[top_n_idx])]
            
        #idx_percent = ints_arr/np.max(ints_arr) > percent_label
        for i, txt in enumerate(np.round(mzs_arr[top_n_idx], 4)):
            plt.annotate(txt, (mzs_arr[top_n_idx][i], ints_arr[top_n_idx][i]), color = 'Black', rotation = 20, fontsize=7)

        plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
        plt.xlabel('m/z')
        plt.ylabel('Counts (-)')
        plt.ylim(bottom=0)
    else:
        plt.text(0.5, 0.5, 'No MS2 data available', fontsize=14, ha='center', va='center')
        plt.axis('off')

    plt.subplot(2,3,3)
    if isinstance(feature['SMILES'], list):
        img = Draw.MolsToGridImage([Chem.MolFromSmiles(x) for x in feature['SMILES']],
                                    legends=[x for x in feature['compound_names']])
        plt.imshow(img)
        plt.axis('off')
    else:
        plt.text(0.5, 0.5, 'No SMILES data available', fontsize=14, ha='center', va='center')
        plt.axis('off')

    plt.subplot(2,3,4)
    if isinstance(sample_names, pd.core.series.Series) or isinstance(sample_names, list):
        # Simple color cycling through matplotlib's default color cycle (consistent with bar plot)
        colors = plt.cm.tab10(np.arange(len(sample_names)) % 10)
        
        has_eic_data = False
        
        for i, sample_name in enumerate(sample_names):
            rt_col = f"{sample_name}_EIC_rt"
            int_col = f"{sample_name}_EIC_intensity"
            
            if rt_col in df.columns and feature[rt_col] is not None:
                rt_array = feature[rt_col]
                intensity_array = feature[int_col]
                if rt_array is not None and intensity_array is not None:
                    plt.plot(rt_array/60, intensity_array, label=sample_name, alpha=0.8, color=colors[i])
                    has_eic_data = True
        
        if has_eic_data:
            plt.xlabel('RT (min)', fontsize=12)
            plt.ylabel('Intensity', fontsize=12)
            plt.legend(fontsize=8)
            plt.title('EICs', fontsize=12)
        else:
            plt.text(0.5, 0.5, 'No EIC data available', fontsize=14, ha='center', va='center')
            plt.axis('off')
    else:
        plt.text(0.5, 0.5, 'No sample names provided', fontsize=14, ha='center', va='center')
        plt.axis('off')

    plt.subplot(2,3,5)
    if isinstance(feature['mzs_isotopes'], list):
        mzs_arr = np.array(feature['mzs_isotopes'])
        ints_arr = np.array(feature['ints_isotopes'])
        plt.stem(mzs_arr, ints_arr, 'blue', markerfmt=" ", basefmt=" ")
        _, stemlines, _ = plt.stem(mzs_arr, ints_arr, 'blue',markerfmt=" ", basefmt=" ")
        plt.setp(stemlines, color = 'blue', linewidth= 3)

        if len(mzs_arr) < top_n_label:
            top_n_idx = np.arange(len(mzs_arr))
        else:
            top_n_idx = np.argpartition(ints_arr, -top_n_label)[-top_n_label:]
            # Sort these indices by intensity for nicer annotation order (optional)
            top_n_idx = top_n_idx[np.argsort(-ints_arr[top_n_idx])]
            
        #idx_percent = ints_arr/np.max(ints_arr) > percent_label
        for i, txt in enumerate(np.round(mzs_arr[top_n_idx], 4)):
            plt.annotate(txt, (mzs_arr[top_n_idx][i], ints_arr[top_n_idx][i]), color = 'Black', rotation = 20, fontsize=7)

        plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
        plt.xlabel('m/z')
        plt.ylabel('Counts (-)')
        plt.ylim(bottom=0)

    plt.suptitle(f'Feature {feature_index} | m/z: {feature['mz']:.4f} | RT: {feature['rt']/60:.2f} min', fontsize=14)
    plt.tight_layout()
    plt.show(block=False)
    #plt.pause(0.1)
    #plt.draw()
    #plt.pause(0.1)

def schedule_plot_overview(df, feature_idx):
    # Use matplotlib’s GUI-safe event loop
    def run():
        feature_overview_plot(df, feature_idx)
    plt.gcf().canvas.manager.window.after(10, run)

def scatter_plot(df, x, y, col=None, highlight=None, size=100, alpha=0.7, cmap_name='coolwarm'):
    
    df = df.copy()
    df['feature_index'] = df.index
   
    _, ax = plt.subplots(figsize=(10, 8))

    sc = ax.scatter(df[x], df[y], s=size, c=df[col] if col else 'darkcyan', cmap=cmap_name, alpha=alpha)
    if isinstance(highlight, pd.core.series.Series):
        ax.scatter(df[x][highlight], df[y][highlight], facecolors='none', edgecolors='red')

    cur = mplcursors.cursor(sc, hover=False)

    @cur.connect("add")
    def on_add(sel):
        idx = sel.index
        feature_idx = df.iloc[idx]['feature_index']
        feature_overview_plot(df, feature_idx)

    ax.set_xlabel(x, fontsize=16)
    ax.set_ylabel(y, fontsize=16)
    plt.tight_layout()
    plt.show()


def kmd_plot(df, 
             diffs='CF2', 
             hs_tol=0.005, 
             n_min=3):
    """
    IN DEVELOPMENT
    BUG: Lines are not correct!
    """
    mz_vec = df['mz'].values

    df_kmd = kmd_analysis(mz_vec=mz_vec,
                          RT_vec=np.zeros(len(mz_vec)),
                          diffs=diffs,
                          hs_tol=hs_tol,
                          n_min=n_min)

    idx_hit = np.where(df_kmd['min_homologues'] == True)[0]
    df_kmd['feature_index'] = df.index

    plt.figure(figsize=(10, 8))
    scatter_list = []
    hs_num = np.unique(df_kmd['hs_number'][idx_hit].values)
    markers = ['o', 's', '^', 'v', 'D', '*', 'p', 'h', 'x', '+'] * len(hs_num)

    for n, hs_n in enumerate(hs_num):
        idx = np.where(df_kmd['hs_number'] == hs_n)[0]
        if len(idx) > 0:
            mz_array = np.array(mz_vec[idx])
            kmd_array = np.array(df_kmd['KMD'].iloc[idx])
            feature_indices = df_kmd['feature_index'].iloc[idx]

            sorted_idx = np.argsort(mz_array)
            sc = plt.scatter(mz_array, kmd_array, 50, marker=markers[n], alpha=0.7)

            # Attach feature index to each point
            sc.feature_indices = feature_indices.values 
            scatter_list.append(sc)

            plt.plot(mz_array[sorted_idx], kmd_array[sorted_idx], alpha=0.5)

    plt.xlabel('m/z', fontsize=16)
    plt.ylabel('KMD', fontsize=16)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.tight_layout()

    cur = mplcursors.cursor(scatter_list, hover=False)

    @cur.connect("add")
    def on_add(sel):
        x_val, y_val = sel.target
        idx = sel.index
        feature_idx = int(sel.artist.feature_indices[idx])

        sel.annotation.set_text(
            f"index={feature_idx}\nx={x_val:.4f}\ny={y_val:.4f}")

        try:
            pyperclip.copy(f"{x_val:.4f}")
            print(f"Copied to clipboard: x = {x_val:.4f}")
        except Exception as e:
            print(f"Clipboard copy failed: {e}")

        feature_overview_plot(df, feature_idx)
    plt.show()


def plot_spectrum_with_differences(masses, 
                                   intensities,
                                   ref_mass, 
                                   tol=0.005, 
                                   title=None, 
                                   show_plot=False, 
                                   reference_spec=None):
    """
    Plots a stem plot of a mass spectrum and annotates positive mass differences from a reference peak.
    
    Parameters:
    - masses: Array of m/z values.
    - intensities: Corresponding intensity values.
    - ref_mass: Mass of the reference peak.
    - tol: Tolerance for considering a mass difference as labeled.
    - labels: Dictionary of named mass differences to annotate (e.g., {'Na': 23.003, 'K': 45.342}).
    """
        
    labels = {'-H': -1.0072, 
              'H':1.0072, 
              'Na': 21.9823, 
              'K': 37.9565, 
              'NH4': 17.0265,
              'Cl': 34.96885,
              'CH3COOH':60.021130, 
              'Br': 78.91833, 
              'COOH': 44.99765, 
              'NaCOOH': 67.987425,
              'H2O':18.01056,
              'CO2':43.98983}

    ref_mass_idx = np.argmin(np.abs(masses - ref_mass)) 

    ref_mass = masses[ref_mass_idx]  # Get the reference mass
    mass_diffs = np.abs(masses - ref_mass)  # Compute absolute mass differences 

    fig, ax = plt.subplots(figsize=(10, 6))

    if reference_spec is not None:
        # normalize spec to max

        reference_spec[1] = reference_spec[1]*np.max(intensities)/np.max(reference_spec[1])

        idx_smaller = reference_spec[0] > (np.min(masses) - 5)

        reference_spec[0] = reference_spec[0][idx_smaller]
        reference_spec[1] = reference_spec[1][idx_smaller]

        if len(reference_spec[0] > 0):

            ax.stem(reference_spec[0], -reference_spec[1], 'violet', markerfmt=" ", basefmt=" ")  # Plot the spectrum
            _, stemlines, _ = ax.stem(reference_spec[0], -reference_spec[1], 'violet',markerfmt=" ", basefmt=" ")
            plt.setp(stemlines, color = 'violet', linewidth = 1)

    ax.stem(masses, intensities, 'black', markerfmt=" ", basefmt=" ")  # Plot the spectrum
    _, stemlines, _ = ax.stem(masses, intensities, 'black',markerfmt=" ", basefmt=" ")
    plt.setp(stemlines, color = 'black', linewidth = 1)

    ax.stem(masses[ref_mass_idx], intensities[ref_mass_idx], 'blue', markerfmt=" ", basefmt=" ")  # Highlight the reference peak

    for i, txt in enumerate(np.round(masses, 4)):
        plt.annotate(txt, (masses[i],intensities[i] - intensities[i] * 0.05), color = 'Black', rotation = 0, fontsize=8)

    # Annotate mass differences
    for i, diff in enumerate(mass_diffs):
        if i == ref_mass_idx or diff <= tol:  # Skip the reference peak and small differences
            continue

        # Check if difference matches any given label
        matched_label = None
        if labels:
            for name, value in labels.items():

                if abs(diff - value) <= tol:  
                    matched_label = name
                    break

        # Choose annotation text: Name if matched, else just the numeric value
        annotation_text = matched_label if matched_label else f"{diff:.4f}"

        ax.annotate(annotation_text,
                    xy=(masses[i], intensities[i]),
                    xytext=(masses[i], intensities[i] + max(intensities) * 0.05),
                    arrowprops=dict(arrowstyle="->", color='red', lw=1),
                    fontsize=8, color='red')
        
    # Labels and title
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")
    ax.axhline(0, color='black', linewidth=1)
    #plt.ylim(ymin = 0)
    ax.set_title(title)

    if show_plot == True:
        plt.show()

    return fig, ax


def smiles_subplot(smiles_list, 
                   name_list, 
                   molsPerRow=4):

    img = Draw.MolsToGridImage([Chem.MolFromSmiles(x) for x in smiles_list],
                               molsPerRow=molsPerRow,
                               subImgSize=(200,200),
                               legends=[x for x in name_list])
    plt.figure()
    plt.imshow(img)
    plt.axis('off')
    plt.tight_layout()
    plt.show()
    #display(img)


def save_smiles_to_svg(smiles: str, 
                       output_path: str, 
                       size: int = 300):
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    drawer = Draw.MolDraw2DSVG(size, size)
    drawer.drawOptions().backgroundColour = None  # Request transparent background
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    
    svg = drawer.GetDrawingText()
    # Remove white background manually (the white <rect> tag)
    svg = svg.replace('fill:#FFFFFF', 'fill:none')

    with open(output_path, 'w') as f:
        f.write(svg)


def isotope_mirror_plot(mz_arr_exp, 
                        ints_arr_exp, 
                        mz_arr_theo, 
                        ints_arr_theor, 
                        title=None):
    
    fig = plt.figure(figsize=(5,4))

    col_exp = 'red'
    col_theor = 'blue'

    plt.stem(mz_arr_exp, ints_arr_exp, col_exp, markerfmt=" ", basefmt=" ")
    _, stemlines, _ = plt.stem(mz_arr_exp, ints_arr_exp, col_exp, markerfmt=" ", basefmt=" ", label='_nolegend_')
    plt.setp(stemlines, color=col_exp, linewidth=5)

    for i, txt in enumerate(np.round(mz_arr_exp, 4)):
        plt.annotate(txt, (mz_arr_exp[i], ints_arr_exp[i]), color=col_exp, rotation=20, fontsize=7, alpha = 0.7)
    
    plt.stem(mz_arr_theo, -ints_arr_theor, col_theor, markerfmt=" ", basefmt=" ")
    _, stemlines, _ = plt.stem(mz_arr_theo, -ints_arr_theor, col_theor, markerfmt=" ", basefmt=" ", label='_nolegend_')
    plt.setp(stemlines, color=col_theor, linewidth=5)

    for i, txt in enumerate(np.round(mz_arr_theo, 4)):
        plt.annotate(txt, (mz_arr_theo[i], -ints_arr_theor[i]), color=col_theor, rotation=20, fontsize=7, alpha = 0.7)

    plt.legend(['exp', 'theo'])
    plt.axhline(y=0, color='k', linestyle='-')
    plt.xlabel('m/z')
    plt.ylabel('Intensity')
    if title:
        plt.title(title)
    plt.tight_layout()

    return fig


def eics_overview(exps, 
                  mzs, 
                  num_labeled_peaks=2, 
                  color_by_sample=True, 
                  eic_colors=None, 
                  sample_colors=None, 
                  sample_idx=None, 
                  extraction_window=0.005, 
                  smoothing=True, 
                  prominence=3000, 
                  height=None, 
                  ms_level=1):
    """
    Plot EICs for all samples and m/z values, or for a specific sample (sample_idx).
    If color_by_sample is True, all EICs from the same sample get the same color (sample_colors).
    If False, each EIC gets its own color (eic_colors).
    Set sample_idx to an integer to plot only that sample (index in exps).
    """

    if sample_colors is None:
        sample_colors = ['blue', 'tomato', 'darkcyan', 'crimson', 'darkmagenta', 'dimgray', 'fuchsia']
    if eic_colors is None:
        eic_colors = ['blue', 'tomato', 'darkcyan', 'crimson', 'darkmagenta', 'dimgray', 'fuchsia']*len(mzs)
    plt.figure()
    c = 0
    # Determine which samples to plot
    if sample_idx is not None:
        sample_indices = [sample_idx]
    else:
        sample_indices = range(len(exps))
    for s in sample_indices:
        exp = exps[s]
        for n in range(len(mzs)):
            rts, ints_orig, _, peak_indices, peak_width, _ = find_peaks(
                exp,
                mzs[n],
                extraction_window=extraction_window,
                smoothing=smoothing,
                prominence=prominence,
                height=height,
                ms_level=ms_level
            )
            rts = rts/60
            if color_by_sample:
                color = sample_colors[s % len(sample_colors)]
            else:
                color = eic_colors[c % len(eic_colors)]
            c += 1
            if len(peak_indices) > 0:
                top_idxs = np.argsort(ints_orig[peak_indices])[-num_labeled_peaks:][::-1]
            else:
                top_idxs = []
            for m in range(len(peak_indices)):
                peak_range = np.arange(peak_width[2][m], peak_width[3][m]).astype(int)
                plt.plot(rts[peak_range], ints_orig[peak_range], color=color, alpha=0.8, linewidth=2)
                if m in top_idxs:
                    plt.text(rts[peak_indices[m]], ints_orig[peak_indices[m]], f'{mzs[n]:.4f}', fontsize = 8, color = color)
    plt.xlabel('RT (min)')
    plt.ylabel('Intensity')
    plt.tight_layout()


def isotope_patters_raw_plot(exp, 
                             mz, 
                             rt, 
                             md_deviation=0.02):

    """
    Function plot a cleaned isotope patter from raw data
    """

    # Extraction from measured raw data
    rt_array = np.array([spec.getRT() for spec in exp])
    # closest index to the target retention time
    idx = np.argmin(np.abs(rt_array - rt))
    # array with MS levels
    ms_level_arr = np.array([spec.getMSLevel() for spec in exp])
    # find indices where the MS1 level is 1
    indices_ms1 = np.where(ms_level_arr == 1)[0]
    # compute the distances to the target index
    distances = np.abs(indices_ms1 - idx)   
    # find the closest 1 (ties are automatically resolved)
    idx = indices_ms1[np.argmin(distances)]

    mz_array = exp[int(idx)].get_peaks()[0]
    ints_array = exp[int(idx)].get_peaks()[1]

    # up to the fifth isotopes
    idx_cut = np.logical_and(mz_array > mz - 0.1, mz_array < mz + 5.1)

    mz_arr_exp = mz_array[idx_cut]
    ints_arr_exp = ints_array[idx_cut]/np.max(ints_array[idx_cut])*100

    # get closest isotopes to perform matching
    md_precursor = mz - np.round(mz)
    md_mz_arr_exp = mz_arr_exp - np.round(mz_arr_exp)

    # max MD deviations: 3*13C = 0.01008, 37Cl = -0.0029, 81Br = -0.0020, 34S = -0.0042
    # assume 0.02 Da deviation to be conservative enough
    idx_md_match = np.abs(md_mz_arr_exp - md_precursor) < md_deviation
    mz_arr_exp_md_cleaned = mz_arr_exp[idx_md_match]
    ints_arr_exp_md_cleaned = ints_arr_exp[idx_md_match]

    plt.figure(figsize=(6,4))
    col = 'blue'

    plt.stem(mz_arr_exp, ints_arr_exp, 'grey', markerfmt=" ", basefmt=" ")
    markerline, stemlines, baseline = plt.stem(mz_arr_exp, ints_arr_exp, 'grey', markerfmt=" ", basefmt=" ", label='_nolegend_')
    plt.setp(stemlines, color='grey', linewidth=2)
    plt.setp(markerline, alpha=0.1)     
    plt.setp(stemlines, alpha=0.1)      
    plt.setp(baseline, alpha=0.1) 

    plt.stem(mz_arr_exp_md_cleaned, ints_arr_exp_md_cleaned, col, markerfmt=" ", basefmt=" ")
    _, stemlines, _ = plt.stem(mz_arr_exp_md_cleaned, ints_arr_exp_md_cleaned, col, markerfmt=" ", basefmt=" ", label='_nolegend_')
    plt.setp(stemlines, color=col, linewidth=2)

    for i, txt in enumerate(np.round(mz_arr_exp_md_cleaned, 4)):
        plt.annotate(txt, (mz_arr_exp_md_cleaned[i], ints_arr_exp_md_cleaned[i]), color=col, rotation=20, fontsize=7, alpha = 0.7)

    plt.xlabel('m/z')
    plt.ylabel('Counts')
    plt.ylim(ymin = 0)
    plt.title(f'm/z = {np.round(mz, 4)} @ RT = {np.round(rt, 1)}')
    plt.tight_layout()
    plt.show()


def plot_feature_eics(df_alignment, 
                      sample_names, 
                      feature_idx, 
                      figsize=(8, 6)):
    """
    Plot EICs for a specific feature across all samples
    
    Parameters:
    df_alignment: DataFrame - the alignment table with EIC columns
    sample_names: list - list of sample names
    feature_idx: int - index of the feature in df_alignment
    figsize: tuple - figure size (width, height)
    """
    row = df_alignment.loc[feature_idx]
    
    plt.figure(figsize=figsize)
    
    colors = plt.cm.tab10(range(len(sample_names)))  # Different colors for each sample
    
    for i, sample_name in enumerate(sample_names):
        rt_col = f"{sample_name}_EIC_rt"
        int_col = f"{sample_name}_EIC_intensity"
        
        if rt_col in df_alignment.columns and row[rt_col] is not None:
            rt_array = row[rt_col]
            intensity_array = row[int_col]
            if rt_array is not None and intensity_array is not None:
                plt.plot(rt_array/60, intensity_array, label=sample_name, alpha=0.8, color=colors[i])

    plt.xlabel('RT (min)')
    plt.ylabel('Intensity')
    
    # Use 'mz' and 'RT' columns or fallback values
    mz_val = row.get('mz', row.get('MZ', 'Unknown'))
    rt_val = row.get('rt', row.get('RT', row.get('Rt', 'Unknown')))

    plt.title(f'FeatureID {feature_idx} - m/z: {mz_val:.4f}, RT: {rt_val:.2f}')
    plt.legend()
    plt.tight_layout()
    plt.show()

# Example usage:
# plot_feature_eics(df_alignment, sample_names, 0)  # Plot EICs for first feature
# plot_feature_eics(df_alignment, sample_names, 10, figsize=(15, 8))  # Custom size


def plot_full_eic(exp, 
                  mass, 
                  extraction_window=0.005, 
                  ms_level=1) -> None:

    """
    Function to directly plot a full EIC using matplotlib.
    """
    rts, ints = get_full_eic(exp, mass, extraction_window, ms_level)

    plt.figure()
    plt.plot(rts, ints)
    plt.title(f'm/z = {np.round(mass, 4)} | extraction width = {extraction_window}')


def mz_rt_plot_interactive(df):
    
    fig = px.scatter(df, 
                     x = df['rt']/60, 
                     y = 'mz',
                     color =  np.log10(df['intens_mean']), color_continuous_scale = px.colors.sequential.Viridis,
                     hover_name = df.index, 
                     hover_data=['mz', 'rt', 'unique_homologues', 'formulas'])

    fig.update_traces(marker=dict(size=20,
                      line=dict(width=1, color='DarkSlateGrey')),
                      selector=dict(mode='markers'),
                      opacity=0.5)
    
    fig.update_layout(xaxis_title="RT (min)", yaxis_title="m/z", font=dict(size=20), showlegend=False)

    return fig.to_html(full_html=False, include_plotlyjs=False)


def mdc_mc_plot_interactive(df):
    
    # constants
    m_CF = -8.40596e-05
    m_CHF = -0.0005237
    intercept_CF = 0.0010087
    intercept_CHF = 0.0229902

    x = np.linspace(np.min(df['m/C']), np.max(df['m/C']), 100)
    y_CF = m_CF * x + intercept_CF
    y_CHF = m_CHF * x + intercept_CHF

    fig1 = px.scatter(df, x='m/C', y='MD/C', color = np.log10(df['intens_mean']),#Df_FeatureData['Score_scaled']
                      color_continuous_scale = px.colors.sequential.Turbo,
                      hover_name = df.index, 
                      hover_data=['mz', 'rt', 'intens_mean', 'unique_homologues', 'formulas'])
    
    fig2 = px.line(x=x, y=y_CF, markers=False)
    fig3 = px.line(x=x, y=y_CHF, markers=False)
    fig4 = px.line(x = np.ones(200)*50, y = np.linspace(np.min(df['MD/C']), np.max(df['MD/C']), 200), markers=False)
    
    fig1.update_traces(marker=dict(size=15,
                      line=dict(width=1, color='DarkSlateGrey')),
                      selector=dict(mode='markers'),
                      opacity=0.5)
    
    fig_all = go.Figure(data=fig1.data + fig2.data + fig3.data + fig4.data)
    fig_all.update_layout(xaxis_title="m/C", yaxis_title="MD/C", font=dict(size=22), showlegend=False)
    
    return fig_all.to_html(full_html=False, include_plotlyjs=False)


def mc_histogram_plot_interactive(df):

    fig = px.histogram(df, x="m/C")
    fig.update_layout(xaxis_title="m/C", 
                      yaxis_title="Counts", 
                      font=dict(size=22), 
                      showlegend=False)
    fig.update_xaxes(range=[0, 100])

    return fig.to_html(full_html=False, include_plotlyjs=False)


def kmd_plot_interactive(df, mC_limit=0):

    df = df[['mz', 'rt', 'unique_homologues', 'hs_number',
                                     'intens_mean', 'min_homologues', 'm/C', 'KMD', 'formulas']]

    HS_pos = df[df['min_homologues'] == True]
    HS_neg = df[df['min_homologues'] == False]
    HS_pos = HS_pos[HS_pos['m/C'] > mC_limit]

    if HS_pos.empty:
        return "<p>No data available for KMD plot.</p>"

    # Selection
    selection = alt.selection_point(on='mouseover', fields=['hs_number'])
    color = alt.condition(selection,
                          alt.Color('hs_number:N', legend=None, scale=alt.Scale(scheme="set1")),
                          alt.value('lightgrey'))

    HS_features = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='mz',
        y='KMD',
        color=color,
        tooltip=["mz", "rt", "unique_homologues", "hs_number", "intens_mean", "formulas"]
    ).properties(width=500, height=500).interactive().add_params(selection)

    HS_negative = alt.Chart(HS_neg).mark_circle(size=100, opacity=0.01).encode(
        x='mz',
        y='KMD',
        color=alt.value('lightgray'),
        tooltip=["mz", "rt"]
    )

    selection2 = alt.selection_point(fields=['hs_number'])
    opacity = alt.condition(selection2, alt.value(0), alt.value(1))
    color_click = alt.condition(selection2,
                                alt.Color('set1:N', legend=None, scale=alt.Scale(scheme="set1")),
                                alt.value('lightgray'))

    HS_click = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='mz',
        y='KMD',
        color=color_click,
        opacity=opacity,
        tooltip=["mz", "rt", "unique_homologues", "hs_number", "intens_mean", 'formulas']
    ).properties(width=500, height=500).interactive().add_params(selection2)

    RT_chart1 = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='rt',
        y='mz',
        color=color
    ).properties(width=500, height=500).add_params(selection)

    RT_chart2 = alt.Chart(HS_pos).mark_circle(size=100).encode(
        x='rt',
        y='mz',
        color=color_click,
        opacity=opacity,
        tooltip=["mz", "rt", "unique_homologues", "hs_number", "intens_mean", "formulas"]
    ).properties(width=500, height=500).interactive().add_params(selection2)

    # Combine charts
    RT_chart = RT_chart1 + RT_chart2
    HS_concat = alt.hconcat(HS_features + HS_negative + HS_click, RT_chart).configure_axis(
        labelFontSize=20,
        titleFontSize=20
    )

    return HS_concat.to_html()


# Old functions, not used right now.

# KMD vs. m/z plot linked together with m/z vs. RT plot
def kmd_plot_interactive_BACKUP(Df_FeatureData, mC_limit = 0):

    Df_FeatureData = Df_FeatureData[['mz','rt','unique_homologues','hs_number', 'intens_mean', 'min_homologues', 'm/C', 'KMD', 'formulas']]

    HS_pos = Df_FeatureData[Df_FeatureData['min_homologues'] == True]
    HS_neg = Df_FeatureData[Df_FeatureData['min_homologues'] == False]

    HS_pos = HS_pos[HS_pos['m/C'] > mC_limit]

    if not HS_pos.empty:
        # Select features of same homologous series on mouseover
        selection = alt.selection_point(on='mouseover', fields=['hs_number'])

        # Set color of features in the same homologous series
        color = alt.condition(selection,
                            alt.Color('hs_number:N', legend=None, scale=alt.Scale(scheme="set1")),
                            alt.value('lightgrey'))

        # Create layer 1: Figure contains features in homologous series
        HS_features = alt.Chart(HS_pos).mark_circle(size=100).encode(
            x='mz',
            y='KMD',
            color=color,
            tooltip=["mz","rt","unique_homologues","hs_number", "intens_mean", "formulas"]
        ).properties(
            width=500,
            height=500  
        ).interactive(
        ).add_params(
            selection
        )

        # Create layer 2: Figure contains features not in homologous series
        HS_negative = alt.Chart(HS_neg).mark_circle(size=100, opacity = 0.01).encode(
            x='mz',
            y='KMD',
            color=alt.value('lightgray'),
            tooltip=["mz","rt"],
        )

        # Create layer 3: Activated when a certain homologous series is clicked
        selection2 = alt.selection_point(fields=['hs_number'])
        opacity = alt.condition(selection2, alt.value(0), alt.value(1))


        color_click = alt.condition(selection2,
                            alt.Color('set1:N', legend=None, scale=alt.Scale(scheme="set1")),
                            alt.value('lightgray'))

        HS_click = alt.Chart(HS_pos).mark_circle(size=100).encode(
            x='mz',
            y='KMD',
            color=color_click,
            opacity=opacity,
            tooltip=["mz","rt","unique_homologues","hs_number", "intens_mean", 'formulas']
        ).properties(
            width=500,
            height=500  
        ).interactive(
        ).add_params(
            selection2
        )
            
        # Combine all 3 layers for first Figure HS_comb 
        HS_comb = HS_features + HS_negative + HS_click
            

        # Create second figure (RT vs. mz)
        # Create layer 1: Highlight on mouseover
        RT_chart1 = alt.Chart(HS_pos).mark_circle(size=100).encode(
            x='rt',
            y='mz',
            color=color
        ).properties(
            height=500,
            width=500
        ).add_params(
            selection
        )

        # Create layer 2: Activate when clicked
        RT_chart2 = alt.Chart(HS_pos).mark_circle(size=100).encode(
            x='rt',
            y='mz',
            color=color_click,
            opacity=opacity,
            tooltip=["mz","rt","unique_homologues","hs_number", "intens_mean", "formulas"]
        ).properties(
            width=500,
            height=500  
        ).interactive(
        ).add_params(
            selection2
        )
            
        # Combine both RT layers
        RT_chart = RT_chart1 + RT_chart2

        # Horizontally concatenate mz vs. KMD and mz vs RT and adjust font size
        HS_concat = alt.hconcat(HS_comb, RT_chart).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        )

    # HS_concat.save(os.path.join(results_folder, f'{output_name}_HS.html'))
    if not HS_pos.empty:
        # your altair chart building code...
        HS_concat = alt.hconcat(HS_comb, RT_chart).configure_axis(
            labelFontSize=20,
            titleFontSize=20
        )
        return HS_concat.to_html()  # <<<<< THIS IS KEY
    else:
        return "<p>No data available for KMD plot.</p>"


# Plotting annotated MS2 spectra (diagnostic fragments and fragment mass differences)
def MS2_spectra_plotter(
        Df_FeatureData,
        idx,
        diffs,
        results_folder,
        output_name,
        font_size = 16
        ):

    fig = go.Figure()
    for n, peak in enumerate(Df_FeatureData['mz_peaks'][idx]):
        fig.add_trace(
            go.Scatter(x = (Df_FeatureData['mz_peaks'][idx][n], Df_FeatureData['mz_peaks'][idx][n]), y = (Df_FeatureData['intens_peaks'][idx][n], 0), 
            mode = 'lines+text', 
            line = dict(color = "black"),
            text = [np.round(Df_FeatureData['mz_peaks'][idx][n], 4)],
            textposition="top center")
            )
    
    for n, peak in enumerate(Df_FeatureData['mz_peaks_diagnostic'][idx]):
        fig.add_trace(
            go.Scatter(x = (Df_FeatureData['mz_peaks_diagnostic'][idx][n], Df_FeatureData['mz_peaks_diagnostic'][idx][n]), y = (Df_FeatureData['intens_peaks_diagnostic'][idx][n], 0),
            mode = 'lines+text',
            line=dict(color="blue", width = 3),
            text = [Df_FeatureData['formula_diagnostic'][idx][n]],
            textposition="bottom right",
            textfont=dict(color="blue"))
            )
    
    for diff in diffs:
        for n, peak in enumerate(Df_FeatureData[f'mz_{diff}'][idx]):
            fig.add_trace(
                go.Scatter(x = (Df_FeatureData[f'mz_{diff}'][idx][n], Df_FeatureData[f'mz_{diff}'][idx][n]), y = (Df_FeatureData[f'intens_{diff}'][idx][n], 0), 
                mode = 'lines+text', 
                line=dict(color="indianred", width = 3),
                text = [np.round(Df_FeatureData[f'mz_{diff}'][idx][n], 4)],
                textposition="top center",
                textfont=dict(color="indianred"))
                )
            
    if type(Df_FeatureData['frag_idx'][idx]) == np.ndarray:

        u, c = np.unique(Df_FeatureData['frag_idx'][idx], return_counts = True)
        prefix_c = np.zeros(len(Df_FeatureData['frag_idx'][idx]))
        for x in range(len(u)):
            prefix_c[np.where(u[x] == Df_FeatureData['frag_idx'][idx])[0]] = np.arange(0, c[x])
        
        for n, peak in enumerate(Df_FeatureData['frag_idx'][idx]):
            fig.add_trace(
                go.Scatter(x = (Df_FeatureData['mz_peaks'][idx][peak], Df_FeatureData['mz_peaks'][idx][peak]),
                           y = (Df_FeatureData['intens_peaks'][idx][peak], 0),
                mode = 'text', 
                line=dict(color="indianred",width = 3, dash ='dash'),
                text = [int(prefix_c[n])* '<br>' + '<br>' + Df_FeatureData['new_formulas'][idx][n]],
                textposition="bottom right",
                textfont=dict(color="indianred"))
                )
            
    fig.update_layout(showlegend=False,
                      template='plotly', #plotly_white, simple_white
                      title={'text': f'm/z = {np.round(Df_FeatureData["mz_msms"][idx],4)} | RT = {np.round(Df_FeatureData["rt_msms"][idx]/60, 2)} | Intensity = {int(Df_FeatureData["intensity"][idx])}'},
                      xaxis_title="m/z",
                      yaxis_title="Counts",
                      xaxis_range=[np.min(Df_FeatureData['mz_peaks'][idx]) - 10, Df_FeatureData["mz"][idx] + 10],
                      font = dict(size = font_size)
                      )
    fig.write_html(os.path.join(results_folder, f'{output_name}_Spec_mz_{str(np.round(Df_FeatureData["mz"][idx], 4))}_intens_{str(np.round(Df_FeatureData["intens_mean"][idx], 0))}.html'))



# m/z vs. RT plot for MSMS alignment validation
def mz_RT_MSMS(
        Df_FeatureData, 
        Df_MS2RawData, 
        idx_in_features, 
        idx_in_MS2RawData,
        Results_folder
        ):
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=Df_FeatureData['rt']/60, y=Df_FeatureData['mz'],
                        mode='markers',
                        marker=dict(color="Navy", size = 15)))
    fig.add_trace(go.Scatter(x=Df_MS2RawData['rt']/60, y=Df_MS2RawData['mz'],
                        mode='markers',
                        marker=dict(color="Green", size = 10)))

    fig.add_trace(go.Scatter(x=Df_MS2RawData['rt'][idx_in_MS2RawData]/60, y=Df_MS2RawData['mz'][idx_in_MS2RawData],
                        mode='markers',
                        marker=dict(color="Orange", size = 8)))
    fig.update_layout(showlegend=False,
                        template='plotly',
                        title={'text': f'{len(idx_in_MS2RawData)} of {len(Df_MS2RawData)} MS2 spectra | {len(np.unique(idx_in_features))} of {len(Df_FeatureData)} Features'},
                        xaxis_title="RT",
                        yaxis_title="m/z",
                        font = dict(size = 20))
    fig.write_html(os.path.join(Results_folder, 'plots/RT_mz_MSMS.html'))


def mass_spec_3d_html(exp, 
                      rt_start=0, 
                      rt_end=10, 
                      mz_min=100, 
                      mz_max=1000, 
                      int_thresh=1000, 
                      ms_level=1) -> None:
    """
    Function to generate 3D representation of an mzML file in a given time window.
    
    Parameters:
    path (str): Path to mzML file
    rt_start (float): Start of the time window in minutes
    rt_end (float): End of the time window in minutes
    mz_min (float): Minimum m/z value
    mz_max (float): Maximum m/z value
    int_thresh (float): Intensity threshold
    ms_level (int): MS level of the spectra to be plotted

    Returns:
    None: The function generates an interactive 3D plot of the mass spectrum using plotly and saves it as an html file.
    """
    rt_array = np.linspace(exp[0].getRT(), exp[exp.getNrSpectra()-1].getRT(), exp.getNrSpectra())
    idx_start = (np.abs(rt_array - rt_start*60)).argmin()
    idx_end = (np.abs(rt_array - rt_end*60)).argmin()

    def getMS1RawData(exp, start_spec, end_spec, ms_level):

        ms1_spec_mz = []
        ms1_spec_intens = []
        for n in range(start_spec, end_spec):
            if exp[n].getMSLevel() == ms_level:
                ms1_spec_mz.append(exp[n].get_peaks()[0])
                ms1_spec_intens.append(exp[n].get_peaks()[1])

        return ms1_spec_mz, ms1_spec_intens

    spec_mz, spec_intens = getMS1RawData(exp, idx_start, idx_end, ms_level)

    rt_dim = []
    t_vec = np.linspace(rt_start, rt_end, len(spec_mz))
    for n in range(len(spec_mz)):
        rt_dim.append(t_vec[n]*np.ones(len(spec_mz[n])))

    def remove_noise(spec_mz, spec_intens, t_dim, int_thresh):

        spec_mz_fil = [None]*len(spec_mz)
        spec_intens_fil = [None]*len(spec_mz)
        t_dim_fil = [None]*len(spec_mz)
        for n in range(len(spec_mz)):
            noise_idx = spec_intens[n] > int_thresh
            spec_mz_fil[n] = spec_mz[n][noise_idx]
            spec_intens_fil[n] = spec_intens[n][noise_idx]
            t_dim_fil[n] = t_dim[n][noise_idx]
        return spec_mz_fil, spec_intens_fil, t_dim_fil

    spec_mz_fil, spec_intens_fil, t_dim_fil = remove_noise(spec_mz, spec_intens, rt_dim, int_thresh)

    def min_max(spec_mz_fil, spec_intens_fil, t_dim_fil, mz_min,  mz_max):

        mz = [None]*len(spec_mz)
        intens = [None]*len(spec_mz)
        t_dim = [None]*len(spec_mz)
        for n in range(len(spec_mz)):
            idx = np.logical_and(spec_mz_fil[n] > mz_min, spec_mz_fil[n] < mz_max)
            mz[n] = spec_mz_fil[n][idx]
            intens[n] = spec_intens_fil[n][idx]
            t_dim[n] = t_dim_fil[n][idx]
        return mz, intens, t_dim

    spec_mz_fil, spec_intens_fil, t_dim_fil = min_max(spec_mz_fil, spec_intens_fil, t_dim_fil, mz_min, mz_max)

    # PLOTLY 3D MASS SPEC DATA
    x = np.concatenate(spec_mz_fil).ravel()
    y = np.concatenate(t_dim_fil).ravel()
    z = np.concatenate(spec_intens_fil).ravel()

    # Idea from: https://community.plotly.com/t/how-to-plot-3d-interactive-stem-plot-in-plotly/75918/2
    stemsx = []
    stemsy = []
    stemsz = []
    for xs, ys in zip(x.flatten(), y.flatten()):
        stemsx.extend([xs, xs, None])
        stemsy.extend([ys, ys, None])
    for zs in z:
        stemsz.extend([0, zs, None])

    fig = go.Figure(go.Scatter3d(x = stemsx, y = stemsy, z = stemsz, 
                                 mode= "lines", line=dict(color="indigo", width=7)))
    
    fig.update_layout(template='simple_white', 
                    scene = dict(xaxis_title='m/z',
                                 yaxis_title='RT (min)',
                                 zaxis_title='Counts'),
                                 font = dict(size = 16))
    fig.update_scenes(aspectmode='manual', aspectratio=dict(x=2, y=1.5, z=1))
    fig.write_html('3DMassSpec.html')
    fig.show()

    # 2D respresentation
    figure = px.scatter(x=x, y=y, color=np.log10(z), color_continuous_scale= 'viridis')


def feature_map_plotter(fm):

    """ 
    Plotting a 3D plot of features from a oms.FeatureMap()
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for feature in fm:
        color = next(ax._get_lines.prop_cycler)['color']
        # chromatogram data is stored in the subordinates of the feature
        for i, sub in enumerate(feature.getSubordinates()):
            retention_times = [x[0] for x in sub.getConvexHulls()[0].getHullPoints()]
            intensities = [int(y[1]) for y in sub.getConvexHulls()[0].getHullPoints()]
            mz = sub.getMetaValue('MZ')
            ax.plot(retention_times, intensities, zs = mz, zdir = 'x', color = color)
            #if i == 0:
            #    ax.text(mz,retention_times[0], max(intensities)*1.02, feature.getMetaValue('label'), color = color)

    ax.set_ylabel('time (s)')
    ax.set_xlabel('m/z')
    ax.set_zlabel('intensity (cps)')
    plt.show()


def show_data_reduction(df, 
                        reduction_measure, 
                        flip=True, 
                        n=200):

    """
    Function to visualizes how many features in df are kept as a threshold is varied.
    """

    if flip == True:
        variation = np.flip(np.linspace(np.min(reduction_measure), np.max(reduction_measure), n))
    else:
        variation = np.linspace(np.min(reduction_measure), np.max(reduction_measure), n)

    fraction = np.zeros(len(variation))
    for i, r in enumerate(variation):
        fraction[i] = len(df[reduction_measure < r]) / len(df) *100


    fig = plt.figure(figsize = (4,4))
    plt.scatter(variation, fraction, color = 'darkcyan', alpha = 0.8)
    plt.xlabel('threshold')
    plt.ylabel('% features considered')
    fig.tight_layout()
    plt.grid()
    plt.show()