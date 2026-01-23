
import json
import numpy as np
from itertools import compress
import pandas as pd

def massbank_to_df(path):

    """
    Load JSON MassBank File into pandas dataframe
    Note: In this function, only HRMS spectra with at least
    two peaks are kept
    """

    with open(path,encoding="utf8") as f:
        data = json.load(f)

    hrms = ['CE-ESI-TOF','ESI-ITFT','ESI-ITTOF','ESI-QTOF','ESI-TOF',\
            'LC-ESI-ITFT','LC-ESI-ITTOF','LC-ESI-QFT', 'LC-ESI-QTOF',\
                'LC-ESI-TOF', "LC-ESI-ITFT"]

    bool_keep = [None]*len(data)
    for n in range(len(data)):
        if data[n]['metaData'][9]['value'] in hrms and \
            data[n]['metaData'][10]['value'] == 'MS2':# and \
            bool_keep[n] = True
        else:
            bool_keep[n] = False

    data = list(compress(data, bool_keep)) # keep only HRMS MS2 spectra

    bool_keep_2 = [None]*len(data)
    for n in range(len(data)):
        if data[n]['spectrum'].count(':') > 1:
            bool_keep_2[n] = True
        else:
            bool_keep_2[n] = False

    data = list(compress(data, bool_keep_2)) # keep only spectra with at least 2 peaks
    # =============================================================================
    # read out precursor m/z's
    prec_mz = np.zeros((len(data)))
    for n in range(len(data)):
        for m in range(len(data[n]['metaData'])):
            if data[n]['metaData'][m]['name'] == 'precursor m/z':
                prec_mz[n] = float(data[n]['metaData'][m]['value'])

    # read out compound names
    cpd_names = [None]*len(data)
    for n in range(len(data)):
        cpd_names[n] = data[n]['compound'][0]['names'][0]['name']

    # read out CEs
    CEs = [None]*len(data)
    for n in range(len(data)):
        for m in range(len(data[n]['metaData'])):
            if data[n]['metaData'][m]['name'] == 'collision energy':
                CEs[n] = data[n]['metaData'][m]['value']

    # read out sum formulas
    sum_formula = [None]*len(data)
    for n in range(len(data)):
        for m in range(len(data[n]['compound'][0]['metaData'])):
            if data[n]['compound'][0]['metaData'][m]['name'] == 'molecular formula':
                sum_formula[n] = data[n]['compound'][0]['metaData'][m]['value']

    # read out SMILES
    smiles = [None]*len(data)
    for n in range(len(data)):
        for m in range(len(data[n]['compound'][0]['metaData'])):
            if data[n]['compound'][0]['metaData'][m]['name'] == 'SMILES':
                smiles[n] = data[n]['compound'][0]['metaData'][m]['value']

    # read out instrument
    instrument = [None]*len(data)
    for n in range(len(data)):
        instrument[n] = data[n]['metaData'][9]['value']


    # read out sum formulas
    fragment_formulas = [None]*len(data)
    for n in range(len(data)):
        try:
            fragment_formulas[n] = [[None for x in range(len(data[n]['annotations']))] for x in range(2)]
        except KeyError:
            fragment_formulas[n] = [[None],[None]]

    # read out sum annotations
    for n in range(len(data)):
        try:
            for m in range(len(data[n]['annotations'])):
                fragment_formulas[n][0][m] = data[n]['annotations'][m]['name']
                fragment_formulas[n][1][m] = float(data[n]['annotations'][m]['value'])
        except KeyError:
            fragment_formulas[n][0] = None
            fragment_formulas[n][1] = None

    # =============================================================================
    # split spectra and convert them to arrays
    spec = [None]*len(data)
    for n in range(len(spec)):
        spec[n] = data[n]['spectrum'].split()

    specs = [None]*len(spec)
    for n in range(len(spec)):
        specs[n] = [sub.split(':') for sub in spec[n]]
        
    specs_array = [None]*len(specs)
    for n in range(len(specs)):
        specs_array[n] = np.array(specs[n]).astype(float)

    # split mz and intensity in two separate arrays
    spec_mz = [None]*len(specs_array)
    spec_intens = [None]*len(specs_array)
    for n in range(len(specs)):    
        spec_mz[n] = specs_array[n][:,0]
        spec_intens[n] = specs_array[n][:,1]


    spec_MD = [None]*len(specs_array)
    for n in range(len(specs)):
        spec_MD[n] = spec_mz[n]-np.round(spec_mz[n])


    Df = pd.DataFrame(data = {'prec_mz': prec_mz,
                            'cpd_names': cpd_names,
                            'CEs': CEs,
                            'sum_formula': sum_formula,
                            'smiles': smiles,
                            'instrument': instrument,
                            'spec_mz': spec_mz,
                            'spec_MD': spec_MD,
                            'spec_intens': spec_intens,
                            'fragment_formulas': fragment_formulas
                            })
    return Df