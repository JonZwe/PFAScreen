import os
import time
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from find_peaks import find_peaks
from integrate_peak import integrate_peak
from match_isotope_patterns_raw import match_isotope_patters_raw
from chromatogram_utils import get_eic
from get_adduct_data import get_adduct_data
from ms_preprocessing_oms import mzml_to_exp, get_ms2_spectrum
from generate_html_with_plots import generate_html_with_plots


def eic_based_suspect_screening(path_output,
                                path_suspect,
                                paths_samples,
                                adduct = '[M+H]+',                      # Format: [M+H]+ [M+Na]+ [M+K]+ [M-H]- [M+Cl]- [M]-
                                topN = 3,                               # if multiple peaks are detected report only topN of the highest peaks, np.inf to report all    
                                extraction_window = 0.004,              # Extration window for EIC (Da)
                                smoothing_peak_detection = False,       # Set to True if smoothing should be done prior to peak detection
                                prominence_peak_detection = 10000,      # Threshold that peak must exceed for being detection (Most important parameter!)
                                baseline_percent = 90,                  # Percentage of intensity used to estimate the baseline (no effect on peak detection)
                                percent_of_max_peak_integration = 3,    # Percentage of peak height used for integration (np.trapz)
                                smoothing_peak_integration = True,      # Set to True if peak should be smoothed prior to integration
                                mz_tol_ms2 = 0.005,                     # mz tolerance for searching after MS2 spectrum 
                                rt_tol_ms2 = 10,                        # rt tolerance for searching after MS2 spectrum 
                                dpi_plots = 80                          # Dpi for saving figures (consider computational performance!
                                ):

    # TODO: improve peak detection! Estimate noise level, and remove it, or remove constant noise level!
    # TODO: Consider directly deleting peaks with intensity = 0, or show only topX hits
    # TODO: Make summary in table, color code in table, Peak: gaussian similarity

    samples = [os.path.basename(p) for p in paths_samples]

    # Directory containing plots
    os.makedirs(path_output, exist_ok=True)
    file_name = f'EIC-based_SS_{samples[0].split(".")[0]}_{adduct}_peak_prom_{prominence_peak_detection}_top{topN}'
    plot_dir = os.path.join(path_output, file_name)
    os.makedirs(plot_dir, exist_ok=True)
    plt.switch_backend('Agg')

    # read suspect list
    df_suspects = pd.read_csv(path_suspect, encoding='ISO-8859-1')

    # NOTE: REMOVE!
    #df_suspects = df_suspects[df_suspects['mix'] == 'FluoroMix'].reset_index(drop=True)

    # calculate adduct mass difference
    adduct_mass_diff, element, polarity, element_bool = get_adduct_data(adduct)

    df_suspects['mz'] = df_suspects['exact_mass'] + adduct_mass_diff

    start = time.time()

    # read in mzML files
    exps = [mzml_to_exp(s) for s in paths_samples]

    cmap = plt.get_cmap('Dark2')  # or 'tab20', 'Set1', etc.
    cols = [cmap(i % cmap.N) for i in range(len(samples))]

    name_col = []
    formula_col = []
    isotope_score_col = []
    rts_col = []
    mz_col = []
    ppm_deviation_col = []
    peak_area_col = []
    peak_height_col = []
    fwhm_col = []
    ms2_col = []

    plots_counter, temporary_counter = 0, 0

    for n in tqdm(range(len(df_suspects)), desc="Processing suspects"):
        #print(f'Finished: {(n+1)/len(df_suspects)*100:.2f}%...')

        rts, ints, ints_smooth, peak_indices, peak_width, bl = find_peaks(exps[0], 
                                                                            mz = df_suspects['mz'][n], 
                                                                            extraction_window = extraction_window, 
                                                                            smoothing = smoothing_peak_detection, 
                                                                            prominence = prominence_peak_detection, 
                                                                            height=None,
                                                                            baseline_percent=baseline_percent,
                                                                            ms_level=1)
        # update plots counter                         
        plots_counter = temporary_counter

        n_peaks_orig = len(peak_indices)

        # retrieve indices of topN peaks according to peak heights
        idxs_topN = np.argsort(ints[peak_indices])[-topN:][::-1]

        peak_indices = peak_indices[idxs_topN]
        peak_width = tuple(arr[idxs_topN] for arr in peak_width)

        #print(f'Reported top {len(peak_indices)} from {n_peaks_orig}')

        for m in range(len(peak_indices)):

            # append data from suspect screening
            rts_col.append(rts[peak_indices[m]])
            mz_col.append(df_suspects['mz'][n])
            name_col.append(df_suspects['name'][n])
            formula_col.append(df_suspects['formula'][n])

            # perform isotope pattern matching
            isotope_dict = match_isotope_patters_raw(exps[0],
                                                    formula = df_suspects['formula'][n],
                                                    element = element,
                                                    element_bool = element_bool,
                                                    rt = rts[peak_indices[m]])

            isotope_score_col.append(isotope_dict['iso_score'])
            ppm_deviation_col.append(isotope_dict['ppm_error'])

            peak_range = np.arange(peak_width[2][m], peak_width[3][m]).astype(int)

            # perform peak integration and FWHM calculation
            peak_area, peak_height, left_idx, right_idx, fwhm = integrate_peak(x = rts[peak_range], 
                                                                            y = ints[peak_range], 
                                                                            smoothing = smoothing_peak_integration,
                                                                            percent_of_max = percent_of_max_peak_integration)
            
            peak_area_col.append(peak_area)
            peak_height_col.append(peak_height)
            fwhm_col.append(fwhm)

            # get MS2 spectrum if present
            mz_ms2, rt_ms2, mz_array_ms2, ints_array_ms2 = get_ms2_spectrum(exp = exps[0], 
                                                                            mz = df_suspects['mz'][n], 
                                                                            rt = rts[peak_indices[m]], 
                                                                            mz_tol = mz_tol_ms2, 
                                                                            rt_tol = rt_tol_ms2)
            if np.isnan(mz_ms2) == False:
                ms2_col.append(True)
            else:
                ms2_col.append(False)

            # ============= Plotting ======================
            plots_counter += 1

            # ================= EIC large =================
            fig, ax = plt.subplots(figsize=(6,4))

            rts_all, ints_all = [], []
            for exp in exps[1:]:
                rts_s, ints_s, _ = get_eic(exp, mass = df_suspects['mz'][n], rt = rts[peak_indices[m]], rt_width = max(rts[peak_range]) - min(rts[peak_range]), 
                                extraction_window = extraction_window, ms_level = 1)
                rts_all.append(rts_s)
                ints_all.append(ints_s)
            #rts_b, ints_b, _ = get_eic(exps[1], mass = df_suspects['mz'][n], rt = rts[peak_indices[m]], rt_width = max(rts[peak_range]) - min(rts[peak_range]), 
            #                       extraction_window = extraction_window, ms_level = 1)
            
            plt.axhline(y = bl, color = 'fuchsia', linestyle = '-', alpha = 0.7)

            for i in range(len(rts_all)):
                plt.plot(rts_all[i], ints_all[i], color=cols[i+1], alpha=0.8, linewidth=2)

            plt.plot(rts[peak_range], ints[peak_range], color=cols[0], alpha=0.8, linewidth=2)

            rts_orig, ints_orig, _ = get_eic(exps[0], mass = df_suspects['mz'][n], rt = rts[peak_indices[m]], rt_width = max(rts[peak_range]) - min(rts[peak_range]), 
                                extraction_window = extraction_window, ms_level = 1)

            if ints_smooth is not None:
                plt.plot(rts[peak_range], ints_smooth[peak_range], color='mediumseagreen', alpha=0.5, linewidth=1)

            plt.plot(rts_orig, ints_orig, color='mediumblue', alpha=0.5, linewidth=2)

            plt.plot(rts[peak_range][left_idx], ints[peak_range][left_idx], marker = 'o', markersize = 5, color='purple')
            plt.plot(rts[peak_range][right_idx],ints[peak_range][right_idx], marker = 'o', markersize = 5, color='purple')

            ax.fill(rts[peak_range][left_idx:right_idx+1], ints[peak_range][left_idx:right_idx+1], color='lightblue',alpha=0.5)

            plt.title(f'{df_suspects["name"][n]} | m/z = {df_suspects["mz"][n]:.4f}')
            plt.xlabel('Retention time (s)')
            plt.ylabel('Intensity')
            plt.tight_layout()

            plt.savefig(os.path.join(plot_dir, f"eic_large_{plots_counter-1}.png"), dpi=dpi_plots)
            plt.close()

            # ================== EIC complete ==================
            plt.figure(figsize=(6,4))
            plt.plot(rts, ints, color='mediumblue', alpha=0.5, linewidth=2)
            plt.plot(rts[peak_indices[m]], ints[peak_indices[m]], marker = 'o', markersize = 4, color='red')
            plt.xlabel('Retention time (s)')
            plt.ylabel('Intensity')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, f"eic_complete_{plots_counter-1}.png"), dpi=dpi_plots)
            plt.close()

            # ================== EIC thumbnail ==================
            fig=plt.figure(figsize=(8,3))
            ax=fig.add_subplot(1,1,1)
            plt.axis('off')
            for i in range(len(rts_all)):
                plt.plot(rts_all[i], ints_all[i], color=cols[i+1], alpha=0.7, linewidth=3)
            plt.plot(rts_orig, ints_orig, color=cols[0], alpha=0.9, linewidth=3)

            extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
            plt.savefig(os.path.join(plot_dir, f"eic_thumb_{plots_counter-1}.png"), dpi=50, bbox_inches=extent)
            plt.close()

            # ================== Isotopes ==================
            plt.figure(figsize=(6,4))
            col_exp = 'royalblue'
            plt.stem(isotope_dict['mz_arr_theo'], isotope_dict['ints_arr_theo'], col_exp, markerfmt=" ", basefmt=" ")
            _, stemlines, _ = plt.stem(isotope_dict['mz_arr_theo'], isotope_dict['ints_arr_theo'], col_exp, markerfmt=" ", basefmt=" ", label='_nolegend_')
            plt.setp(stemlines, color = col_exp, linewidth= 5)

            for i, txt in enumerate(np.round(isotope_dict['mz_arr_theo'], 4)):
                plt.annotate(txt, (isotope_dict['mz_arr_theo'][i],isotope_dict['ints_arr_theo'][i]), color = col_exp, rotation = 20, fontsize=7)

            plt.stem(isotope_dict['mz_arr_exp'], -isotope_dict['ints_arr_exp'], 'purple', markerfmt=" ", basefmt=" ")
            _, stemlines, _ = plt.stem(isotope_dict['mz_arr_exp'], -isotope_dict['ints_arr_exp'], col_exp, markerfmt=" ", basefmt=" ", label='_nolegend_')
            plt.setp(stemlines, color = 'purple', linewidth= 5)

            plt.stem(isotope_dict['mz_arr_exp_score'], -isotope_dict['ints_arr_exp_score'], 'green', markerfmt=" ", basefmt=" ")
            _, stemlines, _ = plt.stem(isotope_dict['mz_arr_exp_score'], -isotope_dict['ints_arr_exp_score'], col_exp, markerfmt=" ", basefmt=" ", label='_nolegend_')
            plt.setp(stemlines, color = 'green', linewidth= 2)

            plt.legend(['theor', 'exper'])
            plt.title(f'Isotope score = {np.round(isotope_dict["iso_score"], 2)*100}%')
            plt.axhline(y = 0, color = 'k', linestyle = '-') 
            plt.xlabel('m/z')
            plt.ylabel('Intensity')
            plt.savefig(os.path.join(plot_dir, f"isotopes_{plots_counter-1}.png"), dpi=dpi_plots)
            plt.close()

            # ================== Structure ==================
            try:
                Draw.MolToFile(Chem.MolFromSmiles(df_suspects['SMILES'][n]),os.path.join(plot_dir, f"structure_{plots_counter-1}.png"))   
            except ValueError:
                plt.imgsave(os.path.join(plot_dir, f"structure_{plots_counter-1}.png"), np.zeros((100,100,3), dtype=np.uint8))

            # ================== MS2 spectrum ================
            if (np.isnan(mz_ms2) == False) and (len(mz_array_ms2) > 0):
                plt.figure(figsize=(6,4))
                plt.stem(mz_array_ms2, ints_array_ms2, 'Black', markerfmt=" ", basefmt=" ")
                markerline, stemlines, baseline = plt.stem(mz_array_ms2, ints_array_ms2, 'Black',markerfmt=" ", basefmt=" ")
                plt.setp(stemlines, color = 'Black', linewidth= 0.5)

                idx_5_percent = ints_array_ms2/np.max(ints_array_ms2) > 0.05
                for i, txt in enumerate(np.round(mz_array_ms2[idx_5_percent], 4)):
                    plt.annotate(txt, (mz_array_ms2[idx_5_percent][i],ints_array_ms2[idx_5_percent][i]), color = 'Black', rotation = 20, fontsize=7)

                plt.ticklabel_format(axis = 'y', style = 'sci', scilimits=(0,0), useMathText=True)
                plt.title(f'MS2: mz = {np.round(mz_ms2, 4)}, RT = {np.round(rt_ms2, 0)} s, Cpd = {df_suspects["name"][n]}')
                plt.xlabel('m/z')
                plt.ylabel('Counts (-)')
                plt.ylim(ymin = 0)
                plt.savefig(os.path.join(plot_dir, f"ms2_spectrum_{plots_counter-1}.png"), dpi=dpi_plots)
                plt.close()
            
        temporary_counter = plots_counter

    plt.switch_backend('TkAgg')

    df_results = pd.DataFrame({'Name': name_col, 
                            'Formula': formula_col,
                            'Isotope score (%)': np.array(isotope_score_col)*100,
                            'm/z': mz_col,
                            'Peak area':peak_area_col,
                            'Peak height':peak_height_col,
                            'ppm deviation': ppm_deviation_col,
                            'RT (min)': np.array(rts_col)/60,
                            'RT (s)': rts_col,
                            'FWHM (s)': fwhm_col,
                            'MS2': ms2_col,
                            f'{adduct}':[adduct_mass_diff]*len(name_col)
                                })

    df_results['m/z'] = df_results['m/z'].round(4)
    df_results['Isotope score (%)'] = df_results['Isotope score (%)'].round(1)
    df_results['RT (min)'] = df_results['RT (min)'].round(2)
    df_results['RT (s)'] = df_results['RT (s)'].round(2)
    df_results['FWHM (s)'] = df_results['FWHM (s)'].round(1)
    df_results['Peak area'] = df_results['Peak area'].astype(int)
    df_results['Peak height'] = df_results['Peak height'].astype(int)
    df_results['ppm deviation'] = df_results['ppm deviation'].round(1)
    df_results[f'{adduct}'] = df_results[f'{adduct}'].round(4)

    df_results = df_results.sort_values(by=['m/z', 'Peak area'], ascending=[True, False])

    df_results.insert(loc=0, column='ID', value=np.arange(1, len(df_results)+1))

    print(f'{df_results["Name"].nunique()} unique suspects found from {len(df_suspects)}')
    print(f'{(np.unique(df_results["Name"].values, return_counts = True)[1] > 1).sum()} suspects detected more than once')

    end = time.time()

    print(f'Took {np.round((end - start)/60, 2)} minutes')

    # Generate HTML
    generate_html_with_plots(
        df_results,
        plot_dir,
        table_name=f'EIC-based suspect screening: {file_name}',
        output_html=os.path.join(path_output, f'{file_name}.html')
    )

    df_results.to_csv(os.path.join(path_output, f'{file_name}.csv'), index=False)


if __name__ == "__main__":

    path_output = r"C:\Users\wxf439\OneDrive\PYTHON\current_dev_codes\ionization_prediction_qnts\manual_rt_all_mixes_neg"

    path_suspect = r"C:\Users\wxf439\OneDrive\PYTHON\current_dev_codes\ionization_prediction_qnts\manual_rt_all_mixes_neg\qnts_mixes_mod_suspectscreen.csv"

    paths_samples = [
        r"C:\Users\wxf439\Documents\MS_raw_data\FluorineID\20251107_first_qNTS_neg\038_FluoroMix_50_ugL_neg_dda_5_inj.mzML",
        r"C:\Users\wxf439\Documents\MS_raw_data\FluorineID\20251107_first_qNTS_neg\037_FluoroMix_10_ugL_neg_dda_5_inj.mzML",
        r"C:\Users\wxf439\Documents\MS_raw_data\FluorineID\20251107_first_qNTS_neg\036_FluoroMix_5_ugL_neg_dda_5_inj.mzML",
        r"C:\Users\wxf439\Documents\MS_raw_data\FluorineID\20251107_first_qNTS_neg\028_Blank_MeOH_neg_dda_5_inj.mzML"
    ]

    eic_based_suspect_screening(path_output,
                                path_suspect,
                                paths_samples,
                                adduct = '[M-H]-',
                                topN = 3,
                                extraction_window = 0.004,
                                smoothing_peak_detection = False,
                                prominence_peak_detection = 10000,
                                baseline_percent = 90,
                                percent_of_max_peak_integration = 3,
                                smoothing_peak_integration = True,
                                mz_tol_ms2 = 0.005,
                                rt_tol_ms2 = 10,
                                dpi_plots = 80
                                )