
import time
from tqdm import tqdm
import pandas as pd

from ms_preprocessing_oms import (mzml_to_exp, 
                                  ms1_feature_finding, 
                                  feature_alignment_oms, 
                                  feature_map_to_df,
                                  average_isotopes_df_align_fm_dfs, 
                                  ms2_spectra_to_df, 
                                  get_all_eics)

from ms2_utils import append_ms2_spec_to_df_align

def run_oms_pipeline(params, 
                     sample_names):

    """
    Run the whole OpenMS pipeline for feature finding, alignment, isotopes, and MS2 data alignment.
    This is based on the parameters.yaml file which specifies all parameters.
    # NOTE: Fix issues with isotopes averaging etc. 
    # Add alignment parameters.
    """

    # Complete pyOpenMS pipeline for feature finding and alignment and isotopes and MS2 data alignment

    start = time.time()

    # Convert config object to dict if needed (for backward compatibility)
    if hasattr(params, 'to_dict'):
        params = params.to_dict()

    # parameters
    mass_error_ppm = params['mass_error_ppm']
    intensity_threshold = params['intensity_threshold']
    min_trace_length = params['min_trace_length']
    isotope_model = params['isotope_model']
    score_by_elements = params['score_by_elements']
    elements = params['elements']
    remove_single_traces = params['remove_single_traces']
    report_convex_hulls = params['report_convex_hulls']
    polarity = params['polarity']

    mz_tol_ms2 = params['mz_tol_ms2']
    rt_tol_ms2 = params['rt_tol_ms2']

    # NOTE: THAT IS NOT GOOD! removes valuable information in high intensity spectra!
    # Use topN approach or similar
    rel_noise_ms2 = params['rel_noise_ms2']
    ms2_avg_tol = params['ms2_avg_tol']
    # NOTE: Add alignment parameters! Move them out of feature_alignment_oms

    df_samples = pd.read_csv(params['path_sample_file'])
    #sample_names = df_samples['sample'] # [os.path.basename(sample).replace(".mzML", "") for sample in df_samples['files']]

    # load mzML files and perform feature detection
    # NOTE: Due to memory issues one file is loaded at a time
    print('Loading files and performing feature detection...')
    fms = []
    dfs_ms2 = []
    sample_eics_dict = {}
    for path, sample_name in tqdm(zip(df_samples['files'], sample_names), desc="Loading files and feature detection"):
        
        exp = mzml_to_exp(path)

        fm = ms1_feature_finding(
            exp,
            mass_error_ppm=mass_error_ppm,
            intensity_threshold=intensity_threshold,
            min_trace_length=min_trace_length,
            isotope_model=isotope_model,
            score_by_elements=score_by_elements,
            elements=elements,
            remove_single_traces=remove_single_traces,
            report_convex_hulls=report_convex_hulls
        )
        fms.append(fm)

        df_ms2 = ms2_spectra_to_df(exp, rel_noise=rel_noise_ms2)
        dfs_ms2.append(df_ms2)

        # Extract EICs before alignment if convex hulls are reported
        if report_convex_hulls == "true":
            print('Extracting EICs...')
            
            rts, ints, max_mz_dev = get_all_eics(exp, fm)

            # Get unique IDs in the same order as the EIC data
            feature_ids = [str(feature.getUniqueId()) for feature in fm]
            
            # Create dictionary using the correct order
            sample_eics = {}
            for i in range(len(feature_ids)):
                unique_id = feature_ids[i]
                sample_eics[unique_id] = {
                    'rt': rts[i],
                    'intensity': ints[i], 
                    'max_mz_dev': max_mz_dev[i]
                }
            
            sample_eics_dict[sample_name] = sample_eics

    df_ms2_all = pd.concat(dfs_ms2, ignore_index=True)


    # perform alignment
    print('Starting alignment...')
    df_alignment, consensus_map = feature_alignment_oms(fms, 
                                                        sample_names)

    df_alignment = df_alignment.rename(columns={'RT': 'rt'})

    # Map EIC data to alignment table if available
    if sample_eics_dict:
        print('Mapping EIC data to alignment table...')
        
        for sample_name in sample_names:
            eic_rt_col = f"{sample_name}_EIC_rt"
            eic_int_col = f"{sample_name}_EIC_intensity"
            eic_mzdev_col = f"{sample_name}_EIC_max_mz_dev"
            
            # Initialize columns with None
            df_alignment[eic_rt_col] = None
            df_alignment[eic_int_col] = None
            df_alignment[eic_mzdev_col] = None
            
            # Fill EIC data where IDs match
            id_col = f"{sample_name}_IDs"
            if id_col in df_alignment.columns:
                print(f'Processing {sample_name}...')
                matches = 0
                total_features = 0
                
                for idx, row in df_alignment.iterrows():
                    if pd.notna(row[id_col]):
                        total_features += 1
                        unique_id = str(row[id_col])
                        
                        if sample_name in sample_eics_dict and unique_id in sample_eics_dict[sample_name]:
                            eic_data = sample_eics_dict[sample_name][unique_id]
                            
                            df_alignment.at[idx, eic_rt_col] = eic_data['rt']
                            df_alignment.at[idx, eic_int_col] = eic_data['intensity']
                            df_alignment.at[idx, eic_mzdev_col] = eic_data['max_mz_dev']
                            matches += 1
        print(f'Added EIC columns for {len(sample_names)} samples to alignment table')

    # retrieve df from FeatureMaps (to have isotope information)
    fm_dfs = [feature_map_to_df(fm) for fm in fms]

    # average isotope patterns
    print('Start appending isotopes...')
    # NOTE MAJOR BOTTLENECK! VERY SLOW CURRENTLY! CONSIDER USING ONLY TOP 5 ISOTOPES!
    # e.g., 6 min for 10000 features
    mzs_isotopes, ints_isotopes = average_isotopes_df_align_fm_dfs(df_alignment, 
                                                                   fm_dfs, 
                                                                   sample_names)

    df_alignment['mzs_isotopes'] = mzs_isotopes
    df_alignment['ints_isotopes'] = ints_isotopes

    # set all row only containing M+0 isotope to nan, get the index of rows where the list length == 1
    #idx = df_alignment[df_alignment['mzs_isotopes'].apply(lambda x: isinstance(x, list) and len(x) == 1)].index
    #df_alignment.loc[idx, ['mzs_isotopes', 'ints_isotopes']] = np.nan

    # load all MS2 data from all mzML files

    print('Append MS2 spectra')
    # average MS2 and append them to df_alignment within given tolerance
    ms2_specs_mz, ms2_specs_ints =  append_ms2_spec_to_df_align(df_alignment, 
                                                                df_ms2_all, 
                                                                mz_tol = mz_tol_ms2, 
                                                                rt_tol = rt_tol_ms2, 
                                                                tol_average = ms2_avg_tol)
    df_alignment['mzs_ms2'] = ms2_specs_mz
    df_alignment['ints_ms2'] = ms2_specs_ints

    print(f'Done! Preprocessing took {(time.time() - start)/60:.2f} mins')

    # THIS NEED TO BE FIXES WITHOUT ALL EXPS LOADED!
    #df_ms2_all = all_ms2_spectra_to_df(exps, 
    #                                   rel_noise = rel_noise_ms2)

    # remove unnecessary columns
    # Drop only '_IDs' columns whose base name exists in sample_names
    cols_to_drop = [col for col in df_alignment.columns if col.endswith('_IDs') and col[:-4] in set(sample_names)]
    df_alignment = df_alignment.drop(columns=cols_to_drop + ['sequence'])
    df_alignment.reset_index(drop=True, inplace=True)


    # NOTE: COMPONENTIZATION NEEDS TO BE IMPLEMENTED!

    df_alignment['adduct'] = '[M-H]-' if polarity == 'neg' else '[M+H]+'

    # This is only a temporary solution

    # This needs to be added into the feature detection part 
    # Ensure that correct adducts are collected into df_alignment
    # Otherwise only [M-H]- or [M+H]+ adducts are considered!!

    # adducts_pos = [b"H:+:0.4",b"Na:+:0.2",b"NH4:+:0.2",b"H-1O-1:+:0.1",b"H-3O-2:+:0.1"]
    # adducts = [b"H-1:-:0.5", b"Cl-1:-:0.4", b"C2H4O2:-:0.1"]

    # print(fm.size())

    # mfd = oms.MetaboliteFeatureDeconvolution()
    # mdf_par = mfd.getDefaults()
    # mdf_par.setValue("negative_mode", "true")
    # mdf_par.setValue("potential_adducts", adducts)
    # mdf_par.setValue("charge_min", -2, "Minimal possible charge")
    # mdf_par.setValue("charge_max", 0, "Maximal possible charge")
    # mdf_par.setValue("charge_span_max", 2)
    # mdf_par.setValue("max_neutrals", 1)
    # mdf_par.setValue("retention_max_diff", 3.0)
    # mdf_par.setValue("retention_max_diff_local", 3.0)
    # mfd.setParameters(mdf_par)
    # fm_adduct = oms.FeatureMap()
    # mfd.compute(fm, fm_adduct, oms.ConsensusMap(), oms.ConsensusMap())

    # print(fm_adduct.size())

    # feature_maps_adducts = []
    # for feature_map in feature_maps:
    #     mfd = oms.MetaboliteFeatureDeconvolution()
    #     mdf_par = mfd.getDefaults()
    #     mdf_par.setValue(
    #         "potential_adducts",
    #         [
    #             b"H:+:0.4",
    #             b"Na:+:0.2",
    #             b"NH4:+:0.2",
    #             b"H-1O-1:+:0.1",
    #             b"H-3O-2:+:0.1",
    #         ],
    #     )
    #     mfd.setParameters(mdf_par)
    #     feature_map_adduct = oms.FeatureMap()
    #     mfd.compute(feature_map, feature_map_adduct, oms.ConsensusMap(), oms.ConsensusMap())
    #     feature_maps_adducts.append(feature_map_adduct)
    # feature_maps = feature_maps_adducts

    return df_alignment