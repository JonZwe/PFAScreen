#%%
# Script to perform PFAS feature prioritization
import os
import numpy as np
import pandas as pd
from MS2_differences_fragments import MS2_differences_fragments # FindPFAS: search in MS/MS raw spectra for PFAS
from KMD_analysis import KMD_analysis                           # Kendrick mass defect analysis
from suspect_screening import suspect_screening                 # Basic suspect screening
import plotting                                                 # Plotting routines

def PFAS_feature_prioritization(
        Df_FeatureData,
        Df_MS2RawData,
        idx_in_features,
        Results_folder,
        diffs = ['C2F4', 'C2F4', 'HF'],
        number_of_fragments = 1,
        mass_tolerance_FindPFAS = 0.002,
        intensity_threshold = 5,
        diffs_KMD = ['CF2'],
        mass_tolerance_KMD = 0.002,
        n_homologues = 3,
        adducts = 1,
        tol_suspect = 0.002,
        save_MSMS_spectra = False,
        mC_cutoff = 0,
        ):
    
    # considers only MS2 data
    Df_FindPFAS = MS2_differences_fragments(
                    Df_MS2RawData['m/z'].to_numpy(),
                    Df_MS2RawData['RT'].to_numpy(),
                    Df_MS2RawData['intens'].to_numpy(),
                    Df_MS2RawData['MS2SpecMZ'].to_list(),
                    Df_MS2RawData['MS2SpecIntens'].to_list(),
                    Df_MS2RawData['idx_MS1'].to_numpy(),
                    diffs = diffs, 
                    number_of_fragments = number_of_fragments, # factor of two for diagnostic fragments
                    mass_tolerance = mass_tolerance_FindPFAS, 
                    intensity_threshold = intensity_threshold,
                    adducts = adducts
                    )

    Df_FeatureData = pd.concat([Df_FeatureData, Df_FindPFAS], axis = 1) # NOTE: Works only of only unique MS1 data is present!!! (see line before FindPFAS)


    # write zeros at columns where MSMS spectra are present but not diffs or dias were found
    idx_MSMS_all = np.unique(idx_in_features) # indices with MSMS
    idx_no_hit = Df_FeatureData.index[Df_FeatureData['n_diffs'].isnull()].to_numpy()
    idx_no_MSMS_hit = idx_no_hit[np.in1d(idx_no_hit, idx_MSMS_all)]
    Df_FeatureData.loc[idx_no_MSMS_hit, 'n_diffs'] = 0
    Df_FeatureData.loc[idx_no_MSMS_hit, 'n_dias'] = 0

    #%%
    # MD/C-m/C, and MD
    # ==============================================================================================
    def calc_MDC_mC(mz, intens_C12, intens_C13):
        # estimate number of carbons per molecule
        C = intens_C13/intens_C12/0.011145
        # calculate mass defect
        MD = mz - np.round_(mz, decimals = 0)
        # calculate m/C and MD/C
        MDC = MD/C
        mC = mz/C

        # specify constants
        mC_CF2 = 49.996806
        MDC_CF2 = -0.003194
        m_CF = -8.40596e-05
        m_CHF = -0.0005237

        # transform data to origin and rotate by the slope of the CHF line
        mC_sr =  (mC - mC_CF2) * np.cos(-m_CHF) - (MDC - MDC_CF2) * np.sin(-m_CHF)
        MDC_sr = (mC - mC_CF2) * np.sin(-m_CHF) + (MDC - MDC_CF2) * np.cos(-m_CHF)

        # calculate radial distance from the CF2 location (with ellipse factor lambda)
        r_CF2 = np.sqrt( (mC_sr/3000)**2 + (MDC_sr)**2 )

        return C, MD, MDC, mC, mC_sr, MDC_sr, r_CF2
        
    Df_FeatureData['C'], Df_FeatureData['MD'], Df_FeatureData['MD/C'], Df_FeatureData['m/C'], _, _, Df_FeatureData['r_CF2'] = calc_MDC_mC(
                                                                                                                                Df_FeatureData['m/z'], 
                                                                                                                                Df_FeatureData['m/z intens'], 
                                                                                                                                Df_FeatureData['m/z+1 intens']
                                                                                                                                )
    Df_FeatureData = Df_FeatureData[Df_FeatureData['m/C'] >= mC_cutoff]

    #%%
    # KMD analysis test
    # ==============================================================================================

    Df_KMD = KMD_analysis(
                mz_vec = Df_FeatureData['m/z'], 
                RT_vec = Df_FeatureData['RT'],
                diffs = diffs_KMD, 
                hs_tol = mass_tolerance_KMD,
                n_min = n_homologues
                )

    Df_FeatureData = pd.concat([Df_FeatureData, Df_KMD], axis=1) 

    print(f'{len(np.unique(Df_FeatureData[Df_FeatureData["min Homologues"] == True]["HS Number"]))} HS detected')

    # Suspect screening
    # ==================================================================================

    Df_susp_list = pd.read_excel('suspect_list.xlsx')

    Df_suspect_screening = suspect_screening(
                                tol_suspect,
                                adducts,
                                Df_FeatureData['m/z'],
                                Df_susp_list,
                                Df_FeatureData
                                )

    Df_suspect_screening.set_index(Df_FeatureData.index, inplace = True) # NOTE: check if this always works!
    Df_FeatureData = pd.concat([Df_FeatureData, Df_suspect_screening], axis = 1)

    # =======================================================================================
    # Prepare to save as Excel file with conditional formatting

    report_excel_list = ['m/z', 'm/z intens', 'RT',  # specify columns that should be saved
                         'n_diffs', 'n_dias', 
                         'C', 'MD', 'MD/C', 'm/C', 
                         'KMD', 'HS Number', 'Unique Homologues', 
                         'hit_in_list', 'FORMULA', 'SMILES']

    Df_FeatureData_Excel = Df_FeatureData[report_excel_list]
    Df_FeatureData_Excel.insert(3, "RT (min)", np.array(Df_FeatureData_Excel['RT']/60), True)
    Df_FeatureData_Excel = Df_FeatureData_Excel.sort_values(by=['m/C'], ascending = False)

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(os.path.join(Results_folder, f'Results_{Results_folder}.xlsx'), engine="xlsxwriter")

    # Convert the dataframe to an XlsxWriter Excel object.
    Df_FeatureData_Excel.to_excel(writer, sheet_name="Sheet1")

    # Get the xlsxwriter workbook and worksheet objects.
    # workbook = writer.book
    worksheet = writer.sheets["Sheet1"]

    # Get the dimensions of the dataframe.
    (max_row, max_col) = Df_FeatureData_Excel.shape

    column_settings = [{"header": column} for column in Df_FeatureData_Excel.columns]

    column_settings.insert(0, {"header" : 'Index'}) # include also index at first position

    # Add the Excel table structure. Pandas will add the data.
    worksheet.add_table(0, 0, max_row, max_col, {"columns": column_settings})

    # Make the columns wider for clarity.
    worksheet.set_column(0, max_col - 1, 12)

    # Apply a conditional format to the required cell range.# NOTE: DOES NOT WORK YET
    worksheet.conditional_format(1, 3, max_row, 3, {"type": "3_color_scale"})

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    # ========================================================================================
    # Plotting routines

    plotting.mz_RT(
                Df_FeatureData,
                Results_folder
                )

    plotting.MDC_mC_plot(
                Df_FeatureData,
                Results_folder
                )
    
    plotting.mC_histogram(
        Df_FeatureData,
        Results_folder
        )

    plotting.KMD_plot(
                Df_FeatureData, 
                Results_folder,
                mC_limit = 0
                )
    
    if save_MSMS_spectra == True:
        for idx in Df_FindPFAS.index:
            plotting.MS2_spectra_plotter(
                        Df_FeatureData,
                        idx,
                        diffs,
                        Results_folder,
                        font_size = 20
                        )
    
    print('Evaluation successfully finished!')

    return Df_FeatureData, Df_FindPFAS