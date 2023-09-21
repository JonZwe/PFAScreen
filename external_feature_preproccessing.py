#%%
import os
import glob
import numpy as np
import pandas as pd
from pyopenms import MSExperiment, MzMLFile
from blank_correc import blank_correc               # Blank correction by m/z, RT, and FoldChange
from ismembertol import ismembertol                 # Array comparison within an absolute tolerance
import plotting                                     # Plotting routines

def external_feature_preproccessing(
        sample_path,
        sample_feature_list_path, 
        blank_feature_list_path,
        mass_tolerance_MSMS_assignment,
        RT_tolerance_MSMS_assignment,
        mass_tolerance_blank_correction,
        RT_tolerance_blank_correction,
        fold_change_blank_correction
        ):

    path = os.path.dirname(sample_feature_list_path)
    sample_file = os.path.basename(sample_feature_list_path)
    blank_file = os.path.basename(blank_feature_list_path)

    # ====================================================================================================
    # Folder generation 
    Results_folder = f'PFAScreen_{sample_file.rsplit(".", 1)[0]}'
    Plots_folder = os.path.join(Results_folder, 'plots')
        
    # check if folder exists, if yes overwrite it
    if not os.path.exists(Results_folder):
        os.makedirs(Results_folder)
    if not os.path.exists(Plots_folder):
        os.makedirs(Plots_folder)
        
    # ====================================================================================================
    # Generate list with concatenated file names
    if not blank_feature_list_path:
        samples = [sample_feature_list_path]
        sample_names = [sample_file]
    else:
        samples = [sample_feature_list_path, blank_feature_list_path]
        sample_names = [sample_file, blank_file]

    # Perform feature finding
    Dfs = []
    for s, sample in enumerate(samples):
        Df = pd.read_excel(sample)
        Dfs.append(Df)

    # %%
    # Optional blank correction if desired
    # ====================================================================================================
    if blank_file != '':
        idx_not_in_blank, _, _ = blank_correc(
                                    Dfs[0]['m/z'], Dfs[0]['RT'], Dfs[0]['m/z intens'], 
                                    Dfs[1]['m/z'], Dfs[1]['RT'], Dfs[1]['m/z intens'], 
                                    mass_tolerance_blank_correction, 
                                    RT_tolerance_blank_correction, 
                                    fold_change_blank_correction
                                    )
        print(f'{len(Dfs[0])} features before blank correction')
        Df_FeatureData = Dfs[0][idx_not_in_blank].reset_index(drop = True)
        print(f'{len(Df_FeatureData)} features after blank correction')
    else:
        Df_FeatureData = Dfs[0]

    Df_FeatureData['RT'] = Df_FeatureData['RT']*60 # NOTE: maybe change RT directly to minutes earlier?
    # %%
    # Get MS2 precursor masses, RT, TIC, and raw spectra from OpenMS MSExperiment and put in DataFrame
    # ====================================================================================================

    exp = MSExperiment()
    MzMLFile().load(sample_path, exp)

    def getMS2RawData(exp):
        precursor_mz = []
        precursor_RT = []
        precursor_intens = []
        MS2_spec_mz = []
        MS2_spec_intens = []
        for spec in exp:
            if spec.getMSLevel() == 2: # NOTE: CONSIDER MULTIPLE CEs HERE!!! OR MERGE SPECTRA!
                precursor_mz.append(spec.getPrecursors()[0].getMZ())
                precursor_RT.append(spec.getRT())
                precursor_intens.append(spec.getPrecursors()[0].getIntensity())
                MS2_spec_mz.append(spec.get_peaks()[0])
                MS2_spec_intens.append(spec.get_peaks()[1])

        precursor_mz = np.array(precursor_mz)
        precursor_RT = np.array(precursor_RT)
        precursor_intens = np.array(precursor_intens)

        return precursor_mz, precursor_RT, precursor_intens, MS2_spec_mz, MS2_spec_intens

    precursor_mz, precursor_RT, precursor_intens, MS2_spec_mz, MS2_spec_intens = getMS2RawData(exp)

    Df_MS2RawData = pd.DataFrame(data = {'m/z': precursor_mz, 
                                         'RT': precursor_RT, 
                                         'intens': precursor_intens, 
                                         'MS2SpecMZ': MS2_spec_mz, 
                                         'MS2SpecIntens': MS2_spec_intens})
    
    # NOTE: REMOVE! ONLY FOR TEST PURPOSES
    import matplotlib.pyplot as plt
    all_MSMS_intens = np.concatenate(Df_MS2RawData['MS2SpecIntens'].values).ravel() # concatenate all MSMS intensities
    plt.hist(np.ma.log10(all_MSMS_intens), bins = 1000)
    plt.title('MS2 Intensity log histogram')
    plt.xlabel('log(Intensity)')
    plt.ylabel('Counts')
    plt.savefig(os.path.join(Results_folder, 'plots/MSMS_histo.png'))


    # check if MSMS_files foldes exists, if True append all MS/MS spectra to Df_MS2RawData
    if os.path.exists(os.path.join(path, 'MSMS_files')):
        MSMS_files = glob.glob(os.path.join(path, 'MSMS_files', '*.mzML'))
        MS2_Dfs = []
        for MSMS_file in MSMS_files:
            exp_MSMS = MSExperiment()
            MzMLFile().load(MSMS_file, exp_MSMS)
            precursor_mz, precursor_RT, precursor_intens, MS2_spec_mz, MS2_spec_intens = getMS2RawData(exp_MSMS)
            MS2_Dfs.append(pd.DataFrame(data = {'m/z':precursor_mz, 
                                                'RT':precursor_RT, 
                                                'intens':precursor_intens, 
                                                'MS2SpecMZ':MS2_spec_mz, 
                                                'MS2SpecIntens':MS2_spec_intens}))
        Df_MS2RawData = pd.concat([Df_MS2RawData] + MS2_Dfs)
        Df_MS2RawData = Df_MS2RawData.reset_index(drop = True)

    # Merge MS2 spectra to MS1 features
    # ====================================================================================================

    _, _, M_bool_mz = ismembertol(Df_FeatureData['m/z'], Df_MS2RawData['m/z'], mass_tolerance_MSMS_assignment)
    _, _, M_bool_RT = ismembertol(Df_FeatureData['RT'], Df_MS2RawData['RT'], RT_tolerance_MSMS_assignment)

    M_combined = np.logical_and(M_bool_mz, M_bool_RT)

    idx_in_MS2RawData =np.where(M_combined)[0]
    idx_in_features = np.where(M_combined)[1]


    plotting.mz_RT_MSMS(
                Df_FeatureData, 
                Df_MS2RawData, 
                idx_in_features, 
                idx_in_MS2RawData, 
                Results_folder
                )

    # add index of MS1 features to MS2 spectra
    Df_MS2RawData['idx_MS1'] = np.nan
    Df_MS2RawData.loc[idx_in_MS2RawData, 'idx_MS1'] = idx_in_features

    # keep only MS2 spectra that have a MS1 feature
    Df_MS2RawData = Df_MS2RawData.dropna(subset=['idx_MS1'])

    # NOTE: ÜBERPRÜFEN! IST DAS KORREKT, UND: WENN MS/MS PRECURSOR MASSEN NICHT IDENTISCH SIND, LÄUFT CODE DANN OHNE DIESEN SCHRITT?
    Df_MS2RawData = Df_MS2RawData.sort_values('intens').drop_duplicates('idx_MS1').sort_index()
    # Df_MS2RawData = Df_MS2RawData.sort_values('intens', ascending=False).drop_duplicates('m/z').sort_index()

    print(f'{len(Df_FeatureData)} MS1 Features')
    print('Evaluation successfully finished!')

    return Df_FeatureData, Df_MS2RawData, idx_in_features

    
    #Df_FeatureData, Df_MS2RawData, idx_in_features = external_feature_preproccessing(
    #        sample_feature_list_path = r'C:\Users\Jonathan\OneDrive\PYTHON_CODES\PyOpenMS\FindPFAS_Workflow_GUI\external_featurefinding\11_Paper_4525_ddMS2_40eV_1it_Inj2_external_feature_list.xlsx', 
    #        blank_feature_list_path = r'C:\Users\Jonathan\OneDrive\PYTHON_CODES\PyOpenMS\FindPFAS_Workflow_GUI\external_featurefinding\09_Blank_Neg_ddMS2_Inj10_external_feature_list.xlsx', 
    #        sample_mzML_path = r'C:\Users\Jonathan\OneDrive\PYTHON_CODES\PyOpenMS\FindPFAS_Workflow_GUI\external_featurefinding\11_Paper_4525_ddMS2_40eV_1it_Inj2.mzML',
    #        mass_tolerance_MSMS_assignment = 0.005,
    #        RT_tolerance_MSMS_assignment = 0.2*60,
    #        mass_tolerance_blank_correction = 0.005,
    #        RT_tolerance_blank_correction = 0.1*60,
    #        fold_change_blank_correction = 10
    #        )