import pyopenms as oms
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict

"""
# Collection of reading/conversion/preprocessing function based on pyOpenMS (for PFAScreen)
# Currently inludes:
# mzml_to_exp, ms1_feature_finding, feature_alignment_oms, feature_map_to_df, load_exps, average_isotopes, all_ms2_spectra_to_df
# Jonathan Zweigle, 06/2025
"""

def mzml_to_exp(path) -> oms.MSExperiment:

    """
    Read mzML file and convert it to OpenMS experiment object
    """
    exp = oms.MSExperiment()
    oms.MzMLFile().load(path, exp)
    exp.sortSpectra(True)

    return exp


def multiple_mzml_to_exps(path_list):
     exps = [mzml_to_exp(path) for path in tqdm(path_list, desc="Loading mzML files")]
     return exps


def get_spectrum(exp, rt, ms_level=1):
    # Takes always the closest spectrum to the given RT in forward direction
    rt_array = np.array([spec.getRT() for spec in exp])
    idx = int(np.argmin(np.abs(rt_array - rt)))
    if exp[idx].getMSLevel() != ms_level:
        while exp[idx].getMSLevel() != ms_level:
            idx += 1
        spec = exp[idx]
    else:
        spec = exp[idx]
    return spec.get_peaks()[0], spec.get_peaks()[1]


def ms1_feature_finding(exp, 
                        mass_error_ppm = 10,                     # 'Allowed mass deviation (in ppm).'
                        intensity_threshold = 500,               # 'Intensity threshold below which peaks are removed as noise.'
                        chrom_peak_snr = 3,                      # 'Minimum intensity above noise_threshold_int (signal-to-noise) a peak should have to be considered an apex.'
                        reestimate_mt_sd = 'true',               # 'Enables dynamic re-estimation of m/z variance during mass trace collection stage.'
                        quant_method = 'area',                   # "Method of quantification for mass traces. For LC 'area' is recommended, 'median' for direct injection data. 'max_height' simply uses the most intense peak in the trace."
                        trace_termination_criterion = 'outlier', # "Termination criterion for the extension of mass traces. In 'outlier' mode, trace extension cancels if a predefined number of consecutive outliers are found (see trace_termination_outliers parameter). In 'sample_rate' mode, trace extension in both directions stops if ratio of found peaks versus visited spectra falls below the 'min_sample_rate' threshold."
                        trace_termination_outliers = 5,          # 'Mass trace extension in one direction cancels if this number of consecutive spectra with no detectable peaks is reached.'
                        min_sample_rate = 0.5,                   # 'Minimum fraction of scans along the mass trace that must contain a peak.'
                        min_trace_length = 3,                    # 'Minimum expected length of a mass trace (in seconds).'
                        max_trace_length = -1,                   # 'Maximum expected length of a mass trace (in seconds). Set to a negative value to disable maximal length check during mass trace detection.'
                        isotope_model = 'none',
                        score_by_elements = 'false',
                        elements = 'CFHNOPSClBr',
                        remove_single_traces = 'false',
                        report_convex_hulls = 'false'
                        ) -> oms.FeatureMap:
    
    """
    Preprocessing: perform ROI, feature detection, and componentization from raw data
    using pyOpenMS. Main algorithms: MassTraceDetection, ElutionPeakDetection, FeatureFindingMetabo.
    """

    # 1)  Mass trace detection
    mass_traces = []
    mtd = oms.MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue("mass_error_ppm", float(mass_error_ppm))           # set according to your instrument mass error
    mtd_params.setValue("noise_threshold_int", float(intensity_threshold)) # adjust to noise level in your data
    mtd_params.setValue("chrom_peak_snr", float(chrom_peak_snr))
    mtd_params.setValue("reestimate_mt_sd", reestimate_mt_sd)
    mtd_params.setValue("quant_method", quant_method)
    mtd_params.setValue("trace_termination_criterion", trace_termination_criterion)
    mtd_params.setValue("trace_termination_outliers", trace_termination_outliers)
    mtd_params.setValue("min_sample_rate", min_sample_rate)
    mtd_params.setValue("min_trace_length", float(min_trace_length))
    mtd_params.setValue("max_trace_length", float(max_trace_length))
    mtd.setParameters(mtd_params)
    mtd.run(exp, mass_traces, 0)

    # 2) Elution peak detection (deconvolution)
    mass_traces_split = []
    mass_traces_final = []
    epd = oms.ElutionPeakDetection()
    epd_params = epd.getDefaults()

    epd_params.setValue("width_filtering", "fixed") # removes mass traces outside min_fwhm of 1 and max_fwhm of 60 s
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)

    if epd.getParameters().getValue("width_filtering") == "auto":
        epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
    else:
        mass_traces_final = mass_traces_split

    # 3) Feature detection (isotope reduction)
    fm = oms.FeatureMap()
    feat_chrom = [] # feature chromatograms
    ffm = oms.FeatureFindingMetabo()
    ffm_params = ffm.getDefaults()
    ffm_params.setValue("isotope_filtering_model", isotope_model)      # metabolites (2% RMS), metabolites (5% RMS)
    ffm_params.setValue("remove_single_traces", remove_single_traces)  # set false to keep features with only one mass trace
    ffm_params.setValue("mz_scoring_by_elements", score_by_elements)
    ffm_params.setValue("report_convex_hulls", report_convex_hulls)
    ffm_params.setValue("elements", elements)
    ffm.setParameters(ffm_params)
    ffm.run(mass_traces_final, fm, feat_chrom)

    fm.setUniqueIds()

    return fm


def feature_alignment_oms(feature_maps,
                          sample_names,
                          mz_distance_ppm = 10,
                          max_num_peaks_considered = -1,
                          pairfinder_mz_unit = "ppm",
                          superimposer_rt_pair_distance_fraction = 0.1,
                          superimposer_mz_pair_max_distance = 0.5,
                          superimposer_num_used_points = 2000,
                          superimposer_scaling_bucket_size = 0.005,
                          grouper_mz_tolerance = None,
                          grouper_mz_unit = None,
                          grouper_rt_tolerance = None,
                          apply_transformed_rt = True) -> tuple:

    """
    Feature alignment by pyOpenMS.

    Main algorithm: MapAlignmentAlgorithmPoseClustering
    Secondary algorithm: FeatureGroupingAlgorithmKD

    Defaults are chosen to preserve current behavior for existing calls.
    Additional optional parameters expose the most important alignment
    and grouping settings without changing legacy usage.
    """

    # (works well if you have a pooled QC for example)
    # use the largest feature map as reference
    ref_index = feature_maps.index(sorted(feature_maps, key=lambda x: x.size())[-1])

    aligner = oms.MapAlignmentAlgorithmPoseClustering()

    # parameter optimization
    aligner_par = aligner.getDefaults()
    aligner_par.setValue("max_num_peaks_considered", int(max_num_peaks_considered))  # infinite if -1
    aligner_par.setValue("pairfinder:distance_MZ:max_difference", float(mz_distance_ppm))  # Never pair features with larger m/z distance
    aligner_par.setValue("pairfinder:distance_MZ:unit", pairfinder_mz_unit)

    # Pose clustering superimposer defaults from OpenMS docs/comments.
    aligner_par.setValue("superimposer:rt_pair_distance_fraction", float(superimposer_rt_pair_distance_fraction))
    aligner_par.setValue("superimposer:mz_pair_max_distance", float(superimposer_mz_pair_max_distance))
    aligner_par.setValue("superimposer:num_used_points", int(superimposer_num_used_points))
    aligner_par.setValue("superimposer:scaling_bucket_size", float(superimposer_scaling_bucket_size))

    # for p in aligner.getParameters().keys():
    #    print(p, aligner.getParameters().getDescription(p))

    # aligner.getParameters().getValue(b'superimposer:mz_pair_max_distance')
    
    # Defaults:
    # b'superimposer:rt_pair_distance_fraction' = 0.1
    # b'superimposer:mz_pair_max_distance' = 0.5
    # b'superimposer:num_used_points' = 2000
    # b'superimposer:scaling_bucket_size' = 0.005 (The scaling of the retention time interval is being hashed into buckets of this size during pose clustering.  A good choice for this would be a bit smaller than the error you would expect from repeated runs.)
    # 

    aligner.setParameters(aligner_par)
    aligner.setReference(feature_maps[ref_index])

    for feature_map in feature_maps[:ref_index] + feature_maps[ref_index + 1 :]:
        trafo = oms.TransformationDescription()  # save the transformed data points
        aligner.align(feature_map, trafo)
        transformer = oms.MapAlignmentTransformer()
        transformer.transformRetentionTimes(feature_map, trafo, bool(apply_transformed_rt))

    feature_grouper = oms.FeatureGroupingAlgorithmKD()
    grouper_par = feature_grouper.getDefaults()
    if grouper_mz_tolerance is not None:
        grouper_par.setValue("distance_MZ:max_difference", float(grouper_mz_tolerance))
    if grouper_mz_unit is not None:
        grouper_par.setValue("distance_MZ:unit", grouper_mz_unit)
    if grouper_rt_tolerance is not None:
        grouper_par.setValue("distance_RT:max_difference", float(grouper_rt_tolerance))
    feature_grouper.setParameters(grouper_par)

    consensus_map = oms.ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()

    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, oms.ColumnHeader())
        file_description.filename = sample_names[i]
        file_description.size = feature_map.size()
        file_descriptions[i] = file_description

    feature_grouper.group(feature_maps, consensus_map)
    consensus_map.setColumnHeaders(file_descriptions)
    consensus_map.setUniqueIds()

    df = consensus_map.get_df()

    # NOTE: This is in principle an optional step, but it is required to
    # later retrieve the isotope information from all feature maps for averaging

    # annotate original feature IDs into consensusmap (Axel Walter)
    fnames = [value.filename for value in consensus_map.getColumnHeaders().values()]
    ids = [[] for _ in fnames]
    for cf in consensus_map:
        fids = {f.getMapIndex(): f.getUniqueId() for f in cf.getFeatureList()}
        for i, fname in enumerate(fnames):
            if i in fids.keys():
                ids[i].append(str(fids[i]))
            else:
                ids[i].append(pd.NA)
    for i, f in enumerate(fnames):
        df[f"{fnames[i]}_IDs"] = ids[i] 

    return df, consensus_map


def feature_map_to_df(feature_map) -> pd.DataFrame:

    """
    Retain Pandas DataFrame from oms.FeatureMap()
    This is needed as the isotopes are not retained 
    during usual conversion to df using fm.get_df()
    """
    # initialize arrays and empty lists
    n_features = feature_map.size()
    mz, rt, mz_area, peak_width = (np.zeros(n_features) for _ in range(4))
    feature_id, mzs_isotopes, ints_isotopes = ([None]*n_features for _ in range(3))

    for n, feature in enumerate(feature_map):
        mz[n] = feature.getMZ()
        rt[n] = feature.getRT()
        mz_area[n] = feature.getIntensity()
        peak_width[n] = feature.getWidth()
        feature_id[n] = feature.getUniqueId()

        # retrieve isotopes into lists
        mzs_iso, ints_iso = [], []
        for m in range(len(feature.getMetaValue('masstrace_centroid_mz'))):
            mzs_iso.append(feature.getMetaValue('masstrace_centroid_mz')[m])
            ints_iso.append(feature.getMetaValue('masstrace_intensity')[m])
        
        mzs_isotopes[n] = mzs_iso
        ints_isotopes[n] = ints_iso

    df = pd.DataFrame(data = {'mz': mz, 
                              'rt': rt, 
                              'mz_area': mz_area, 
                              'peak_width': peak_width,
                              'unique_id': feature_id, 
                              'mzs_isotopes': mzs_isotopes, 
                              'ints_isotopes': ints_isotopes})

    df['unique_id'] = df['unique_id'].astype(str)

    return df


def load_exps(df_sample_info, 
              sample_names, 
              files_column = 'file_path', 
              sample_column = 'file_name') -> list:
    """
    Load several experiments from a DataFrame containing sample information.
    """

    sample_paths = df_sample_info[files_column][df_sample_info[sample_column].isin(sample_names)]

    exps = [mzml_to_exp(sample) for sample in sample_paths]

    return exps


def average_isotopes(mzs_arrays, 
                     ints_arrays, 
                     mz_tolerance=0.02, 
                     iso_spacing=1.0) -> tuple:
    """
    Average m/z and intensity per isotope position only if spacing is ~1 Da apart.
    
    Parameters:
        mzs_arrays: list of np.arrays of m/z values
        ints_arrays: list of np.arrays of intensity values
        mz_tolerance: how close m/z values must be to consider them equal
        iso_spacing: expected delta between isotopes (~1 Da)
    
    Returns:
        (avg_mzs, avg_ints): lists of averaged values
    """

    mz_by_pos = defaultdict(list)
    int_by_pos = defaultdict(list)

    for mzs, ints in zip(mzs_arrays, ints_arrays):
        min_len = min(len(mzs), len(ints))

        for i in range(min_len):
            # Only include if this isotope spacing looks valid (~1 Da)
            if i == 0 or abs(mzs[i] - mzs[i-1] - iso_spacing) < mz_tolerance:
                mz_by_pos[i].append(mzs[i])
                int_by_pos[i].append(ints[i])

    avg_mzs = [np.mean(mz_by_pos[i]) for i in sorted(mz_by_pos)]
    avg_ints = [np.mean(int_by_pos[i]) for i in sorted(int_by_pos)]

    return avg_mzs, avg_ints


def average_isotopes_df_align_fm_dfs(df_alignment, 
                                     fm_dfs, 
                                     sample_names, 
                                     mz_tol = 0.02, 
                                     iso_spacing = 1) -> tuple:

    """
    Average isotopes from feature maps based on a DataFrame alignment.
    NOTE: This is a major bottleneck! It is very slow currently!
    Consider only averaging top N isotopes!
    """

    mzs_isotopes = []
    ints_isotopes = []

    for i, idx in tqdm(enumerate(df_alignment.index), desc="Averaging isotopes"):
        prec_ints = []
        mzs_arrays = []
        ints_arrays = []

        for n, sample in enumerate(sample_names):
            df = fm_dfs[n][fm_dfs[n]['unique_id'] == df_alignment[f'{sample_names[n]}_IDs'].iloc[i]]

            if not df.empty:
                row = df.iloc[0]
                prec_ints.append(row['mz_area'])
                mzs_arrays.append(np.array(row['mzs_isotopes']))

                ints = np.array(row['ints_isotopes'])
                ints_arrays.append(ints / np.max(ints))

        if mzs_arrays and ints_arrays:
            mzs_avg, ints_avg = average_isotopes(mzs_arrays, ints_arrays, mz_tolerance=mz_tol, iso_spacing=iso_spacing)
            ints_avg = list(np.array(ints_avg) / np.max(np.array(ints_avg)))
        else:
            mzs_avg = [np.nan]
            ints_avg = [np.nan]

        mzs_isotopes.append(mzs_avg)
        ints_isotopes.append(ints_avg)

    return mzs_isotopes, ints_isotopes


def ms2_spectra_to_df(exp, 
                      rel_noise = 0.03) -> pd.DataFrame:

    """
    Extract all MS2 spectra an OpenMS experiment
    and convert them to a Pandas DataFrame.
    NOTE: rel_noise is not a good choise for high intensity spectra,
    as it will remove too many peaks. 
    It should have an option e.g. using topN approach.
    """

    # NOTE: Here a fix is needed for DIA data!
    # In case of DIA, no spectra are saved currently

    prec_mz, prec_rt, prec_ints, mzs_arr_ms2, ints_arr_ms2 = [], [], [], [], []
    for spec in exp:
        if spec.getMSLevel() == 2:
            try:
                prec_mz.append(spec.getPrecursors()[0].getMZ())
                prec_rt.append(spec.getRT())
                prec_ints.append(spec.getPrecursors()[0].getIntensity())

                mzs = spec.get_peaks()[0]
                ints = spec.get_peaks()[1]

                if len(ints) > 0:

                    idx = ints/np.max(ints) > rel_noise

                    mzs_arr_ms2.append(mzs[idx])
                    ints_arr_ms2.append(ints[idx])

                else: 
                    mzs_arr_ms2.append(mzs)
                    ints_arr_ms2.append(ints)
            except IndexError:
                continue

    df_ms2 = pd.DataFrame(data = {'prec_mz':prec_mz, 
                                  'prec_rt':prec_rt, 
                                  'prec_ints':prec_ints, 
                                  'mzs_arr_ms2':mzs_arr_ms2, 
                                  'ints_arr_ms2':ints_arr_ms2})

    return df_ms2


def all_ms2_spectra_to_df(exps, 
                          rel_noise = 0.03) -> pd.DataFrame:

    """
    Extract all MS2 spectra from a list of OpenMS experiments
    and convert them to a Pandas DataFrame.
    NOTE: rel_noise is not a good choise for high intensity spectra,
    as it will remove too many peaks. It should have an option e.g. using topN approach.
    """

    dfs = []
    for exp in exps:

        prec_mz, prec_rt, prec_ints, mzs_arr_ms2, ints_arr_ms2 = [], [], [], [], []
        for spec in exp:
            if spec.getMSLevel() == 2:
                prec_mz.append(spec.getPrecursors()[0].getMZ())
                prec_rt.append(spec.getRT())
                prec_ints.append(spec.getPrecursors()[0].getIntensity())

                mzs = spec.get_peaks()[0]
                ints = spec.get_peaks()[1]

                if len(ints) > 0:

                    idx = ints/np.max(ints) > rel_noise

                    mzs_arr_ms2.append(mzs[idx])
                    ints_arr_ms2.append(ints[idx])

                else: 
                    mzs_arr_ms2.append(mzs)
                    ints_arr_ms2.append(ints)

        df = pd.DataFrame(data = {'prec_mz':prec_mz, 'prec_rt':prec_rt, 'prec_ints':prec_ints, 
                                  'mzs_arr_ms2':mzs_arr_ms2, 'ints_arr_ms2':ints_arr_ms2})
        dfs.append(df)
    
    df_ms2 = pd.concat(dfs, ignore_index=True)

    return df_ms2


def get_ms2_spectrum(exp, 
                     mz, 
                     rt, 
                     mz_tol = 0.005, 
                     rt_tol = 10):

    """
    Extract specific MS2 spectrum from an oms.MSExperiment object
    requires mz, rt, mz_tol, rt_tol
    It takes the MS2 spectrum closest to the specified rt (s)
    """
    
    # get ms2 mz, rt, intens
    precursor_mz_arr = []
    precursor_rt_arr = []
    scan_idx = []
    for n, spec in enumerate(exp):
        if spec.getMSLevel() == 2: 
            precursor_mz_arr.append(spec.getPrecursors()[0].getMZ())
            precursor_rt_arr.append(spec.getRT())
            scan_idx.append(n)

    precursor_mz_arr = np.array(precursor_mz_arr)
    precursor_rt_arr = np.array(precursor_rt_arr)
    scan_idx = np.array(scan_idx)

    idx_mz = np.where(np.abs(precursor_mz_arr - mz) < mz_tol)[0]

    if len(idx_mz) > 0:

        idx_close_rt = idx_mz[np.argmin(np.abs(precursor_rt_arr[idx_mz] - rt))]

        if (precursor_rt_arr[idx_close_rt] > rt - rt_tol) and (precursor_rt_arr[idx_close_rt] < rt + rt_tol): 
        
            idx_scan_target = int(scan_idx[idx_close_rt])

            rt_ms2 = exp[idx_scan_target].getRT()
            mz_ms2 = exp[idx_scan_target].getPrecursors()[0].getMZ()
            mz_array_ms2 = exp[idx_scan_target].get_peaks()[0]
            ints_array_ms2 = exp[idx_scan_target].get_peaks()[1]

        else:
            rt_ms2 = np.nan
            mz_ms2 = np.nan
            mz_array_ms2 = np.nan
            ints_array_ms2 = np.nan

    else:
        rt_ms2 = np.nan
        mz_ms2 = np.nan
        mz_array_ms2 = np.nan
        ints_array_ms2 = np.nan

    return mz_ms2, rt_ms2, mz_array_ms2, ints_array_ms2


def get_all_eics(exp, 
                 fm):

    """
    Extract all EICs from a given experiment and feature map.
    Works only if the feature map has convex hulls (setting during feature detection).
    Currenly df_dict is the major bottleneck for performance.
    Could be improved if necessary.
    """

    df_exp = exp.get_df(long='True')
    df_dict = {(rt, mz): inty for rt, mz, inty in zip(df_exp['RT'], df_exp['mz'], df_exp['inty'])}

    rts, ints = [], []
    max_mz_dev = []
    
    for f in fm:
        if not f.getConvexHulls() or len(f.getConvexHulls()) == 0:
            rts.append(np.array([f.getRT()]))
            ints.append([np.nan])
            max_mz_dev.append(0.0)
            continue
            
        arr = f.getConvexHulls()[0].getHullPoints()
        keys = [tuple(row) for row in arr]
        inty = [df_dict.get(key, np.nan) for key in keys]

        rts.append(arr[:, 0])
        ints.append(inty)
        max_mz_dev.append(np.max(arr[:, 1]) - np.min(arr[:, 1]))

    return rts, ints, max_mz_dev


def get_all_eics_original(exp, fm):

    """
    Extract all EICs from a given experiment and feature map.
    Currently df_dict is the major bottleneck for performance.
    Could be improved if necessary.
    ORIGINAL UNTOUCHED VERSION - BACKUP
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


def fill_feature_map(mz, rt, intensity):

    feature_map = oms.FeatureMap()
    for n in range(len(mz)):
        feature = oms.Feature()
        feature.setMZ(mz[n])
        feature.setCharge(1)
        feature.setRT(rt[n])
        feature.setIntensity(intensity[n])

        feature_map.push_back(feature)

    return feature_map


# rts, ints = get_all_eics(exp, fm)
# plt.figure()
# for rt, inty in zip(rts, ints):
#     plt.plot(rt, inty, alpha=0.5)
# plt.show()

# from eic_utils import get_eic_start_stop
# def get_all_eics_2(exp, fm, extraction_window=0.005):
#
#    """
#    Different implementation
#    NOTE: Much lower performance than get_all_eics
#    """
#
#    df_fm = fm.get_df()
#
#    rts, ints = [], []
#    for n in range(len(df_fm)):
#        rt, i, _ = get_eic_start_stop(exp,
#                                   mass = df_fm['mz'].iloc[n],
#                                   rt_start = df_fm['RTstart'].iloc[n],
#                                   rt_stop = df_fm['RTend'].iloc[n],
#                                   extraction_window = extraction_window,
#                                   ms_level = 1)
#        rts.append(rt)
#        ints.append(i)
#
#    return rts, ints