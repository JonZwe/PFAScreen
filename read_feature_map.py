import numpy as np
import pandas as pd

def read_feature_map(feature_map):

    '''
    Function to read data from a pyOpenMS FeatureMap object into a pandas DataFrame
    Data includes: m/z, RT, and m/z and peak area of detected M+1 and M+2 isotopes and
    chromatograpic peak width
    '''

    data = np.zeros((feature_map.size(), 8))
    for n, feature in enumerate(feature_map):
        data[:,0][n] = feature.getMZ() # m/z
        data[:,1][n] = feature.getMetaValue('masstrace_centroid_mz')[1] # m/z + 1
        # check if more than two isotopes are present
        if len(feature.getMetaValue('masstrace_centroid_mz')) > 2: 
            data[:,2][n] = feature.getMetaValue('masstrace_centroid_mz')[2] # m/z + 2
            data[:,5][n] = feature.getMetaValue('masstrace_intensity')[2]   # peak area of M+2
        else:
            data[:,2][n] = np.nan
            data[:,5][n] = np.nan
        data[:,3][n] = feature.getIntensity() # peak area of M
        data[:,4][n] = feature.getMetaValue('masstrace_intensity')[1] # peak area of M+1
        data[:,6][n] = feature.getRT()    # RT
        data[:,7][n] = feature.getWidth() # Chromatographic peak width (s)

    Df = pd.DataFrame(data, columns=['mz', 'mz+1', 'mz+2', 'mz_area', 'mz+1_area', 'mz+2_area', 'rt', 'peak_width'])

    return Df