# Function to read data from pyOpenMS FeatureMap to pandas DataFrame

import numpy as np
import pandas as pd

def read_feature_map(feature_map):

    data = np.zeros((feature_map.size(), 7))
    for n, feature in enumerate(feature_map):
        data[:,0][n] = feature.getMZ() # m/z
        data[:,1][n] = feature.getMetaValue('masstrace_centroid_mz')[1] # m/z + 1
        # check if more than two isotopes are present
        if len(feature.getMetaValue('masstrace_centroid_mz')) > 2: 
            data[:,2][n] = feature.getMetaValue('masstrace_centroid_mz')[2] # m/z + 2
            data[:,5][n] = feature.getMetaValue('masstrace_intensity')[2]   # intensity + 2
        else:
            data[:,2][n] = np.nan
            data[:,5][n] = np.nan
        data[:,3][n] = feature.getIntensity() # intensity M
        data[:,4][n] = feature.getMetaValue('masstrace_intensity')[1] # intensity + M+1
        data[:,6][n] = feature.getRT() # RT

    Df = pd.DataFrame(data, columns=['m/z', 'm/z+1', 'm/z+2', 'm/z intens', 'm/z+1 intens', 'm/z+2 intens', 'RT'])

    return Df