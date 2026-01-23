import os
import glob
import pandas as pd
import numpy as np
import xml.etree.ElementTree as ET
from utils import match_mz_rt_intens
from pyteomics import mgf
import matplotlib.pyplot as plt
from ms2_utils import merge_spectra

"""
Collection of reading/writing functions (for PFAScreen)
Currently including:
df_to_excel_csv_pickle
write_mgf
write_extended_mgf
msdial_peaktable_to_df
msdial_aligmenttable_to_df
masshunter_cef_to_df
Jonathan Zweigle, 06/2025
"""

def df_to_excel_csv_pickle(df, 
                           output_path):

    # Create a Pandas Excel writer using XlsxWriter as the engine.
    writer = pd.ExcelWriter(f'{output_path}.xlsx', engine = "xlsxwriter")

    # Convert the dataframe to an XlsxWriter Excel object.
    df.to_excel(writer, sheet_name = "Sheet1")

    # Get the xlsxwriter workbook and worksheet objects.
    # workbook = writer.book
    worksheet = writer.sheets["Sheet1"]

    # Get the dimensions of the dataframe.
    (max_row, max_col) = df.shape

    column_settings = [{"header": column} for column in df.columns]

    column_settings.insert(0, {"header" : 'Index'}) # include also index at first position

    # Add the Excel table structure. Pandas will add the data.
    worksheet.add_table(0, 0, max_row, max_col, {"columns": column_settings})

    # Make the columns wider for clarity.
    worksheet.set_column(0, max_col - 1, 12)

    # Close the Pandas Excel writer and output the Excel file.
    writer.close()

    # Save results additionally as .csv file
    df.to_csv(f'{output_path}.csv')
    # And as pickle
    df.to_pickle(f'{output_path}.pkl')


def write_mgf(df, 
              output_path, 
              polarity = 'pos'):
    """
    Generate standard MGF file from DataFrame containing the alignment table data
    NOTE: The given df should already be filtered to only contain features with associated MS2 data
    """

    with open(output_path, 'w') as f:
        for n in range(len(df)):
            f.write("BEGIN IONS\n")
            f.write(f"TITLE=MS/MS of feature {df.index[n]} at RT of {df['rt'].iloc[n]}\n")
            f.write(f"PEPMASS={df['mz'].iloc[n]}\n")
            if polarity == 'pos':
                f.write(f"CHARGE=1\n")
            else:
                f.write(f"CHARGE=-1\n")
            f.write(f"Num peaks={len(df['mzs_ms2'].iloc[n])}\n")
            for mz, intensity in zip(df['mzs_ms2'].iloc[n], df['ints_ms2'].iloc[n]):
                f.write(f"{mz} {intensity}\n")
            f.write("END IONS\n\n")


def write_extended_mgf(df, 
                       output_path, 
                       polarity = 'pos'):
    """
    Generate extended MGF file from DataFrame containing the alignment table data
    This is in the style of the SIRIUS export file from MZmine 4
    NOTE: The given df should already be filtered to only contain features with associated MS2 data
    """

    with open(output_path, 'w') as f:
        for n in range(len(df)):
            f.write("BEGIN IONS\n")
            f.write(f"FEATURE_ID={df.index[n]}\n")
            f.write(f"MSLEVEL=1\n")
            f.write(f"RTINSECONDS={df['rt'].iloc[n]}\n")
            f.write(f"PEPMASS={df['mz'].iloc[n]}\n")
            if polarity == 'pos':
                f.write(f"CHARGE=1\n")
            else:
                f.write(f"CHARGE=-1\n")
            f.write(f"SPECTYPE='CORRELATED MS')\n")
            f.write(f"FILENAME=dummy.mzML\n")
            f.write(f"SCANS=dummy\n")
            f.write(f"Num peaks={len(df['mzs_isotopes'].iloc[n])}\n")
            for mz, intensity in zip(df['mzs_isotopes'].iloc[n], df['ints_isotopes'].iloc[n]):
                f.write(f"{mz} {intensity}\n")
            f.write("END IONS\n\n")
            f.write("BEGIN IONS\n")
            f.write(f"FEATURE_ID={df.index[n]}\n")
            f.write(f"MSLEVEL=2\n")
            f.write(f"RTINSECONDS={df['rt'].iloc[n]}\n")
            f.write(f"PEPMASS={df['mz'].iloc[n]}\n")
            if polarity == 'pos':
                f.write(f"CHARGE=1\n")
            else:
                f.write(f"CHARGE=-1\n")
            f.write(f"SCANS=dummy\n")
            f.write(f"Num peaks={len(df['mzs_ms2'].iloc[n])}\n")
            for mz, intensity in zip(df['mzs_ms2'].iloc[n], df['ints_ms2'].iloc[n]):
                f.write(f"{mz} {intensity}\n")
            f.write("END IONS\n\n")


def string_colon_spectrum_as_list(data_str):
    """
    Convert a string of colon-separated pairs into a list of floats.
    This is how spectra are saved in MSDIAL and other software.

    Parameter
    ----------
    data_str : str
        A string of colon-separated pairs, e.g. "1.0:2.0 3.0:4.0" 

    Returns
    -------
    mz_list : list
        The list of mass-to-charge ratios
    ints_list : list
        The list of intensities
    """
    # split the string of colon pairs into a list of strings
    pairs = [item.split(":") for item in data_str.split()]

    # convert the to lists and append
    mz_list = [float(p[0]) for p in pairs]
    ints_list = [float(p[1]) for p in pairs]

    return mz_list, ints_list


# This script contains functions to convert MSDIAL peaktable and alignment table files into Pandas DataFrames.
# Currently it is suitable for version 5.5.250221 of MSDIAL.

def msdial_peaktable_to_df(path):
    """
    Function to read and convert peaktable file (.txt) from MSDIAL 
    to a Pandas DataFrame. This function also removes unused columns
    and converts MS1 isotopes into lists.
    MSDIAL version: 5.5.250221

    Parameters
    ----------
    path : str
        The path to the peaktable file from MSDIAL.

    Returns
    -------
    df : pd.DataFrame
        A Pandas DataFrame containing the cleaned data.
    """
    # read the peaktable file
    df_msdial = pd.read_csv(path, sep='\t')

    # adjust here if the column names become different
    header = ['Precursor m/z', 
              'RT (min)', 
              'RT left(min)',
              'RT right (min)',
              'Height',
              'Area',
              'Gaussian similarity',
              'Adduct',
              'Isotope',
              'MS1 isotopes']

    # remove unused columns
    df = df_msdial[header]

    # perform operation only for features with MS1 isotopes
    df = df[df['Isotope'] > 0].reset_index(drop=True)

    # loop through the MS1 isotopes and put them into a list format
    mzs_isotopes, ints_isotopes = zip(*[
        string_colon_spectrum_as_list(x) if pd.notna(x) else (np.nan, np.nan)
        for x in df['MS1 isotopic spectrum']
    ])

    # put the list into the DataFrame columns
    df['mzs_isotopes'] = list(mzs_isotopes)
    df['ints_isotopes'] = list(ints_isotopes)

    # convert RT to seconds
    df['rt'] = df['RT (min)']*60
    df['rt_left'] = df['RT left(min)']*60
    df['rt_right'] = df['RT right (min)']*60
    df['rt_width'] = df['rt_right'] - df['rt_left']

    # remove the columns that are not needed anymore
    df = df.drop(columns=['RT left(min)', 'RT right (min)', 'RT (min)', 'MS1 isotopes', 'Isotope'])

    # rename the columns to the PFAScreen standard
    df.rename(columns={'Precursor m/z': 'mz', 
                       'Area': 'mz_area',
                       'Adduct': 'adduct',
                       'Height': 'mz_height',
                       'Gaussian similarity': 'gaussian_similarity'}, 
                       inplace=True)
    return df


def msdial_alignmenttable_to_df(path, 
                                skiprows=4):
    """
    Function to read alignment table from MSDIAL and convert it into a pandas DataFrame
    MSDIAL version: 5.5.250221
    Parameters
    ----------
    path : str
        Path to the alignment table file
    skiprows : int, optional
        Number of rows to skip at the beginning of the file, defaults to 4
    Returns
    -------
    df : pandas.DataFrame
        DataFrame containing the alignment table data
    """
    # read aligment table file
    df = pd.read_csv(path, sep='\t', skiprows=skiprows)

    # generate lists of mz and intensity values from the MS1 isotopic spectrum, and MS2 spectrum
    mzs_isotopes, ints_isotopes = zip(*[
        string_colon_spectrum_as_list(x) if pd.notna(x) else (np.nan, np.nan)
        for x in df['MS1 isotopic spectrum']
    ])

    mzs_ms2, ints_ms2 = zip(*[
        string_colon_spectrum_as_list(x) if pd.notna(x) else (np.nan, np.nan)
        for x in df['MS/MS spectrum']
    ])

    # put the list into the DataFrame columns
    df['mzs_isotopes'] = list(mzs_isotopes)
    df['ints_isotopes'] = list(ints_isotopes)
    df['mzs_ms2'] = list(mzs_ms2)
    df['ints_ms2'] = list(ints_ms2)

    # rename the columns to the PFAScreen standard
    df.rename(columns={'Average Mz': 'mz', 
                       'Average Rt(min)': 'rt',
                       'Adduct type': 'adduct',
                       'MS/MS assigned': 'msms_assigned',
                       'S/N ratio average': 'snr',
                       'Spectrum reference file name': 'spectrum_file'},
                       inplace=True)

    # convert RT to seconds
    df['rt'] = df['rt']*60

    # convert msms_assigned to boolean
    df['msms_assigned'] = df['msms_assigned'].astype(bool)

    # remove the columns that are not needed anymore
    # added currently: 'Metabolite name', 
    df = df.drop(columns=['Alignment ID', 'Post curation result', 
                          'Fill %', 'Reference RT', 'Reference m/z', 'Formula', 'Ontology',
                          'INCHIKEY', 'SMILES', 'Annotation tag (VS1.0)', 'RT matched',
                          'm/z matched', 'MS/MS matched', 'Comment', 'Manually modified for quantification',
                          'Manually modified for annotation', 'Isotope tracking parent ID', 
                          'Isotope tracking weight number', 'RT similarity', 'm/z similarity', 
                          'Simple dot product', 'Weighted dot product', 'Reverse dot product', 
                          'Matched peaks count', 'Matched peaks percentage', 
                          'Total score', 'MS1 isotopic spectrum','MS/MS spectrum'],  
                          errors='ignore') # allows also reading of older versions
    return df

def msdial_49_peaktable_to_df(path_peaktable):

    """
    Old function for reading MSDIAL version 4.9 peaktables, has no clear logic and should not be used
    NOTE: Fix pandas warning messages
    """

    df = pd.read_csv(path_peaktable, sep='\t')

    df = df[['PeakID','Precursor m/z', 'RT (min)', 'Area','Isotope', 'Comment']]

    isotope_nr = [int(iso[-1]) for iso in df['Isotope']]  
    df.loc[:,'isotope_nr'] = isotope_nr

    df_isotopes = df[df['isotope_nr'] > 0]

    isotope_id = [int(comment.replace('isotope of ', '')) for comment in df_isotopes['Comment']]
    df_isotopes.loc[:,'isotope_id'] = isotope_id

    df.loc[:, 'mz+1'] = np.nan
    df.loc[:, 'mz+1_area'] = np.nan
    df.loc[:, 'mz+2'] = np.nan
    df.loc[:, 'mz+2_area'] = np.nan
    df.loc[:, 'mz+3'] = np.nan
    df.loc[:, 'mz+3_area'] = np.nan
    df.loc[:, 'mz+4'] = np.nan
    df.loc[:, 'mz+4_area'] = np.nan
    df.loc[:, 'mz+5'] = np.nan
    df.loc[:, 'mz+5_area'] = np.nan
    for n in range(len(df_isotopes)):
        if df_isotopes['isotope_nr'].iloc[n] == 1:
            idx = df['PeakID'] == df_isotopes['isotope_id'].iloc[n]
            df.loc[idx, 'mz+1'] = df_isotopes['Precursor m/z'].iloc[n]
            df.loc[idx, 'mz+1_area'] = df_isotopes['Area'].iloc[n]
        elif df_isotopes['isotope_nr'].iloc[n] == 2:
            idx = df['PeakID'] == df_isotopes['isotope_id'].iloc[n]
            df.loc[idx, 'mz+2'] = df_isotopes['Precursor m/z'].iloc[n]
            df.loc[idx, 'mz+2_area'] = df_isotopes['Area'].iloc[n]
        elif df_isotopes['isotope_nr'].iloc[n] == 3:
            idx = df['PeakID'] == df_isotopes['isotope_id'].iloc[n]
            df.loc[idx, 'mz+3'] = df_isotopes['Precursor m/z'].iloc[n]
            df.loc[idx, 'mz+3_area'] = df_isotopes['Area'].iloc[n]
        elif df_isotopes['isotope_nr'].iloc[n] == 4:
            idx = df['PeakID'] == df_isotopes['isotope_id'].iloc[n]
            df.loc[idx, 'mz+4'] = df_isotopes['Precursor m/z'].iloc[n]
            df.loc[idx, 'mz+4_area'] = df_isotopes['Area'].iloc[n]
        elif df_isotopes['isotope_nr'].iloc[n] == 5:
            idx = df['PeakID'] == df_isotopes['isotope_id'].iloc[n]
            df.loc[idx, 'mz+5'] = df_isotopes['Precursor m/z'].iloc[n]
            df.loc[idx, 'mz+5_area'] = df_isotopes['Area'].iloc[n]

    df.rename(columns={'Precursor m/z': 'mz', 'RT (min)': 'rt', 'Area': 'mz_area'}, inplace=True)

    return df


def msdial_49_aligmenttable_to_df(path_alignmenttable, 
                                  path_peaklists, 
                                  mz_tol = 0.005, 
                                  rt_tol = 0.2, 
                                  rel_int_tol = 0.005):
    """
    Old function for combine isotope data MSDIAL version 4.9 peaktables to the alignmenttable, has no clear logic and should not be used
    Very, slow and prone to small errors in isotope intensity
    NOTE: Fix pandas warning messages
    """

    df_align = pd.read_csv(path_alignmenttable, skiprows=4, sep='\t')

    path_peaklists = os.path.join(path_peaklists, '*.txt')
    for f in glob.glob(path_peaklists):

        df = msdial_49_peaktable_to_df(f)

        #sample = os.path.basename(f).split('.')[0]
        sample =  os.path.basename(f).replace(".txt", "")

        # write isotopes in alignment table
        idx_in_a, idx_in_b, _ = match_mz_rt_intens(df['mz'].values, 
                                                   df['rt'].values,
                                                   df['mz_area'].values,
                                                   df_align['Average Mz'].values, 
                                                   df_align['Average Rt(min)'].values,
                                                   df_align[sample].values, 
                                                   mz_tol = mz_tol, 
                                                   rt_tol = rt_tol, 
                                                   rel_int_tol = rel_int_tol)

        # NOTE: change in the future to circument pandas warning messages
        df_align[f'mz_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+1_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+2_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+3_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+4_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+5_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz_area_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+1_area_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+2_area_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+3_area_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+4_area_{sample}'] = pd.Series(dtype=str)
        df_align[f'mz+5_area_{sample}'] = pd.Series(dtype=str)

        df_align[f'mz_{sample}'].iloc[idx_in_a] = df['mz'].iloc[idx_in_b]
        df_align[f'mz+1_{sample}'].iloc[idx_in_a] = df['mz+1'].iloc[idx_in_b]
        df_align[f'mz+2_{sample}'].iloc[idx_in_a] = df['mz+2'].iloc[idx_in_b]
        df_align[f'mz+3_{sample}'].iloc[idx_in_a] = df['mz+3'].iloc[idx_in_b]
        df_align[f'mz+4_{sample}'].iloc[idx_in_a] = df['mz+4'].iloc[idx_in_b]
        df_align[f'mz+5_{sample}'].iloc[idx_in_a] = df['mz+5'].iloc[idx_in_b]
        df_align[f'mz_area_{sample}'].iloc[idx_in_a] = df['mz_area'].iloc[idx_in_b]
        df_align[f'mz+1_area_{sample}'].iloc[idx_in_a] = df['mz+1_area'].iloc[idx_in_b]
        df_align[f'mz+2_area_{sample}'].iloc[idx_in_a] = df['mz+2_area'].iloc[idx_in_b]
        df_align[f'mz+3_area_{sample}'].iloc[idx_in_a] = df['mz+3_area'].iloc[idx_in_b]
        df_align[f'mz+4_area_{sample}'].iloc[idx_in_a] = df['mz+4_area'].iloc[idx_in_b]
        df_align[f'mz+5_area_{sample}'].iloc[idx_in_a] = df['mz+5_area'].iloc[idx_in_b]

        df_align[f'm+1/m_{sample}'] = df_align[f'mz+1_area_{sample}']/df_align[f'mz_area_{sample}']
        df_align[f'm+2/m_{sample}'] = df_align[f'mz+2_area_{sample}']/df_align[f'mz_area_{sample}']
        df_align[f'm+3/m_{sample}'] = df_align[f'mz+3_area_{sample}']/df_align[f'mz_area_{sample}']
        df_align[f'm+4/m_{sample}'] = df_align[f'mz+4_area_{sample}']/df_align[f'mz_area_{sample}']
        df_align[f'm+5/m_{sample}'] = df_align[f'mz+5_area_{sample}']/df_align[f'mz_area_{sample}']

    samples = [os.path.basename(f).replace(".txt", "") for f in glob.glob(path_peaklists)]

    # NOTE: Change to median or max?!
    df_align['m'] = df_align[[f'mz_{sample}' for sample in samples]].mean(axis=1, skipna=True)
    df_align['m+1'] = df_align[[f'mz+1_{sample}' for sample in samples]].mean(axis=1, skipna=True)
    df_align['m+2'] = df_align[[f'mz+2_{sample}' for sample in samples]].mean(axis=1, skipna=True)
    df_align['m+3'] = df_align[[f'mz+3_{sample}' for sample in samples]].mean(axis=1, skipna=True)
    df_align['m+4'] = df_align[[f'mz+4_{sample}' for sample in samples]].mean(axis=1, skipna=True)
    df_align['m+5'] = df_align[[f'mz+5_{sample}' for sample in samples]].mean(axis=1, skipna=True)

    df_align['m/m'] = 1
    df_align['m+1/m'] = df_align[[f'm+1/m_{sample}' for sample in samples]].median(axis=1, skipna=True)
    df_align['m+2/m'] = df_align[[f'm+2/m_{sample}' for sample in samples]].median(axis=1, skipna=True)
    df_align['m+3/m'] = df_align[[f'm+3/m_{sample}' for sample in samples]].median(axis=1, skipna=True)
    df_align['m+4/m'] = df_align[[f'm+4/m_{sample}' for sample in samples]].median(axis=1, skipna=True)
    df_align['m+5/m'] = df_align[[f'm+5/m_{sample}' for sample in samples]].median(axis=1, skipna=True)

    #df_align['m+1/m_std'] = df_align[[f'm+1/m_{sample}' for sample in samples]].std(axis=1, skipna=True)

    mzs_isotopes = []
    ints_isotopes = []
    for n in range(len(df_align)):
        if pd.notna(df_align['m+1/m'].iloc[n]):
            mzs_isotopes.append(df_align[['m', 'm+1', 'm+2', 'm+3', 'm+4', 'm+5']].iloc[n].dropna().to_list())
            ints_isotopes.append(df_align[['m/m', 'm+1/m', 'm+2/m', 'm+3/m', 'm+4/m', 'm+5/m']].iloc[n].dropna().to_list())
        else:
            mzs_isotopes.append([df_align['m'].iloc[n]])
            ints_isotopes.append([1])

    df_align['mzs_isotopes'] = mzs_isotopes
    df_align['ints_isotopes'] = ints_isotopes

    mzs_ms2, ints_ms2 = zip(*[
        string_colon_spectrum_as_list(x) if pd.notna(x) else (np.nan, np.nan)
        for x in df_align['MS/MS spectrum']
    ])

    df_align['mzs_ms2'] = list(mzs_ms2)
    df_align['ints_ms2'] = list(ints_ms2)

    cols_to_drop = []
    for i in range(1, 6):
        cols_to_drop += [f'm+{i}/m_{sample}' for sample in samples]

    for i in range(1, 6):
        cols_to_drop += [f'mz+{i}_area_{sample}' for sample in samples]

    for i in range(1, 6):
        cols_to_drop += [f'mz+{i}_{sample}' for sample in samples]

    # Drop them from df_align
    df_align = df_align.drop(columns=cols_to_drop)
    df_align = df_align.drop(columns='SMILES')

    # rename the columns to the PFAScreen standard
    df_align.rename(columns={'Average Mz': 'mz', 
                             'Average Rt(min)': 'rt',
                             'Adduct type': 'adduct',
                             'MS/MS assigned': 'msms_assigned',
                             'S/N ratio average': 'snr',
                             'Spectrum reference file name': 'spectrum_file'},
                              inplace=True)
    
    return df_align


def get_isotope_candidates(mzs, 
                           intensities, 
                           precursor_mz,
                           ppm_tol=5, 
                           z_max=3,
                           n_carbons=2, 
                           n_halogen_steps=3):
    """
    Keep peaks that look like isotope candidates:
      - C13 M+1 (1.003355/z)
      - Halogen-type M+2n (1.997/z steps)
    Returns (mzs_filtered, intensities_filtered)
    """
    def ppm_to_da(mz, ppm):
        return mz * ppm / 1e6

    mzs = np.asarray(mzs)
    intensities = np.asarray(intensities)

    # ensure sorted order by m/z
    order = np.argsort(mzs)
    mzs, intensities = mzs[order], intensities[order]

    keep_idx = set()

    def find_near(target_mz):
        diff = np.abs(mzs - target_mz)
        idx = np.argmin(diff)
        if diff[idx] <= ppm_to_da(target_mz, ppm_tol):
            return idx
        return None

    for z in range(1, z_max+1):
        # --- carbon M+1, M+2, etc. (up to n_carbons)
        for k in range(1, n_carbons+1):
            target = precursor_mz + k * (1.003355 / z)
            idx = find_near(target)
            if idx is not None:
                keep_idx.add(idx)

        # --- halogen-type M+2, +4, +6 etc.
        for k in range(1, n_halogen_steps+1):
            target = precursor_mz + k * (1.997 / z)
            idx = find_near(target)
            if idx is not None:
                keep_idx.add(idx)

    # always include the monoisotopic if present
    idx0 = find_near(precursor_mz)
    if idx0 is not None:
        keep_idx.add(idx0)

    # output filtered arrays
    keep_idx = sorted(list(keep_idx))
    return mzs[keep_idx], intensities[keep_idx]


def has_M1_isotope(masses, delta=1.003355, tol=0.005):
    """
    masses: list of m/z values (e.g., [195.066177, 199.061157])
    delta: theoretical M+1 spacing
    tol: allowed error tolerance
    """
    masses = sorted(masses)
    mono = masses[0]
    
    # expected M+1 m/z
    expected = mono + delta
    
    # check if ANY peak is within tolerance of expected
    return any(abs(m - expected) <= tol for m in masses[1:])


def mzmine_to_df(config,
                 path, 
                 name):
    
    # workes for mzmine 4.8.0

    path_full_feature_table_csv = os.path.join(path, f'{name}_full_feature_table.csv')

    df_alignment = pd.read_csv(path_full_feature_table_csv)
    
    def count_spectra(filepath):
        with open(filepath) as f:
            return sum(1 for line in f if line.startswith("BEGIN IONS"))

    path_sirus_mgf = os.path.join(path, f'{name}_sirius.mgf')

    n_spectra = count_spectra(path_sirus_mgf)

    columns = ['feature_id', 'mz_prec', 'rt_prec', 'intens_prec', 'mzs_ms1', 'ints_ms1', 'mzs_ms2', 'ints_ms2']
    df_spectra = pd.DataFrame(np.nan, index=range(n_spectra), columns=columns)

    df_spectra = df_spectra.astype({
        'mzs_ms1': 'object', 'ints_ms1': 'object', 
        'mzs_ms2': 'object', 'ints_ms2': 'object'})

    with mgf.MGF(path_sirus_mgf) as spectra:
        for n, spec in enumerate(spectra):
                df_spectra.at[n, 'feature_id'] = int(spec['params']['feature_id'])
                df_spectra.at[n, 'mz_prec'] = spec['params']['pepmass'][0]
                df_spectra.at[n, 'rt_prec'] = float(spec['params']['rtinseconds'])
                df_spectra.at[n, 'intens_prec'] = float(spec['params']['feature_ms1_height'])

                if spec['params'].get('mslevel') == '1':
                    df_spectra.at[n, 'mzs_ms1'] = spec['m/z array']
                    df_spectra.at[n, 'ints_ms1'] = spec['intensity array']
                elif spec['params'].get('mslevel') == '2':
                    df_spectra.at[n, 'mzs_ms2'] = spec['m/z array']
                    df_spectra.at[n, 'ints_ms2'] = spec['intensity array']


    df_alignment[['feature_id', 'mzs_isotopes', 'ints_isotopes', 'mzs_ms2', 'ints_ms2', 'prec_mz']] = np.nan
    df_alignment = df_alignment.astype({
        'mzs_isotopes': 'object', 'ints_isotopes': 'object', 
        'mzs_ms2': 'object', 'ints_ms2': 'object'})

    df_spectra['feature_id'] = df_spectra['feature_id'].astype(int)

    for n in range(len(df_alignment)):
        
        df_curr = df_spectra[df_spectra['feature_id'] == df_alignment['id'].iloc[n]]

        if not df_curr.empty:

            df_alignment.at[n, 'prec_mz'] = df_curr['mz_prec'].iloc[0]

            # note, it happens that the isotope spacing is not 1! Either fix it here, or make sure PFAScreen handles it
            mzs_isotopes, ints_isotopes = get_isotope_candidates(df_curr['mzs_ms1'].dropna().iloc[0], 
                                                                                        df_curr['ints_ms1'].dropna().iloc[0], 
                                                                                        df_curr['mz_prec'].iloc[0],
                                                                                        ppm_tol=5, 
                                                                                        z_max=3,
                                                                                        n_carbons=2, 
                                                                                        n_halogen_steps=3)
            if len(mzs_isotopes) > 1:
                df_alignment.at[n, 'mzs_isotopes'] = list(mzs_isotopes)
                df_alignment.at[n, 'ints_isotopes'] = list(ints_isotopes)
            
            if len(df_curr['mzs_ms2'].dropna()) > 0:
                mzs_ms2, ints_ms2 = merge_spectra(df_curr['mzs_ms2'].dropna().to_list(),
                                                                df_curr['ints_ms2'].dropna().to_list(),
                                                                tolerance=0.005)
                df_alignment.at[n, 'mzs_ms2'] = mzs_ms2
                df_alignment.at[n, 'ints_ms2'] = ints_ms2

    print(f'{df_alignment["mzs_ms2"].notna().sum()} features from {len(df_alignment)} have MS2 spectra associated.')
    print(f'{df_alignment["mzs_isotopes"].notna().sum()} features from {len(df_alignment)} have isotopes')

    df_alignment["has_M1"] = df_alignment["mzs_isotopes"].apply(
    lambda x: has_M1_isotope(x) if isinstance(x, list) else False
    )
    df_alignment.loc[df_alignment["has_M1"] == False,
                    ["mzs_isotopes", "ints_isotopes"]] = np.nan
    
    print(f'After [M+1] isotope check: {df_alignment["mzs_isotopes"].notna().sum()} isotope patterns remain.')
    
    # Rename columns to PFAScreen standard
    # Not a particularily good solution, as file names can differ from sample names (but ok for now)
    df_sample_list = pd.read_csv(config['path_sample_file'])
    renaming_cols = [f'datafile:{os.path.split(s)[1]}:area' for s in df_sample_list['files'].tolist()]
    sample_names = df_sample_list['sample'].tolist()
    df_alignment = df_alignment.rename(columns=dict(zip(renaming_cols, sample_names)))

    polarity = config['polarity']
    # NOTE: COMPONENTIZATION NEEDS TO BE IMPLEMENTED!
    df_alignment['adduct'] = '[M-H]-' if polarity == 'neg' else '[M+H]+'

    return df_alignment

    
def masshunter_cef_to_df(path): 
    """
    Function to read in MassHunter .CEF files from MolecularFeatureExtraction to a pandas DataFrame
    NOTE: Not extensively tested.
    """
    tree = ET.parse(path)   # parse .CEF file
    root = tree.getroot()   # generate elementtree
    
    # loot to get data from elementtree in numpy matrix
    data = np.zeros((len(root[0]), 6))
    for n in range(len(root[0])):
        data[:,0][n] = root[0][n][0].attrib.get('m')
        data[:,1][n] = root[0][n][0].attrib.get('rt')
        data[:,2][n] = root[0][n][0].attrib.get('v')
        if root[0][n][3][0].attrib.get('p') == '-':
            for s in range(len(root[0][n][3][3])):
                if root[0][n][3][3][s].attrib.get('s') == 'M-H':
                    data[:,3][n] = root[0][n][3][3][s].attrib.get('x')  # M-H
                    data[:,4][n] = root[0][n][3][3][s].attrib.get('v')  # volume M
            for m in range(len(list(root[0][n][3][3]))):
                if root[0][n][3][3][m].attrib.get('s') == 'M-H+1':
                    data[:,5][n] = root[0][n][3][3][m].attrib.get('v') # volume M+1               
        elif root[0][n][3][0].attrib.get('p') == '+':
            for s in range(len(root[0][n][3][3])):
                if root[0][n][3][3][s].attrib.get('s') == 'M+H':
                    data[:,3][n] = root[0][n][3][3][s].attrib.get('x')  # M+H
                    data[:,4][n] = root[0][n][3][3][s].attrib.get('v')  # volume M
            for m in range(len(list(root[0][n][3][3]))):
                if root[0][n][3][3][m].attrib.get('s') == 'M+H+1':
                    data[:,5][n] = root[0][n][3][3][m].attrib.get('v') # volume M+1
    # create a DataFrame
    df = pd.DataFrame(data, columns=['mass', 'RT', 'Vol',  'M+-H', 'Vol_M+-H', 'Vol_M+-H+1'])
    df = df[['M+-H', 'RT', 'Vol_M+-H', 'Vol_M+-H+1']]
    df.rename({'M+-H': 'm/z', 'Vol_M+-H': 'm/z intens','Vol_M+-H+1':'m/z+1 intens'}, axis=1, inplace=True)

    return df