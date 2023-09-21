# Function to read data from featureXML files

import xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import json

def read_featureXML(path):

    tree = ET.parse(path)   # parse featureXML
    root = tree.getroot()
    
    data = np.zeros((len(root[0]), 5))
    for n in range(len(root[0])):
        for m in range(len(root[0][n])):
            if root[0][n][m].attrib.get('name') == 'masstrace_centroid_mz':
                data[:,0][n] = json.loads(root[0][n][m].attrib.get('value'))[0] # m/z
                data[:,1][n] = json.loads(root[0][n][m].attrib.get('value'))[1] # m/z + 1
            if root[0][n][m].attrib.get('name') == 'masstrace_intensity':
                data[:,2][n] = json.loads(root[0][n][m].attrib.get('value'))[0] # intens m/z
                data[:,3][n] = json.loads(root[0][n][m].attrib.get('value'))[1] # intens m/z + 1
            if root[0][n][m].attrib.get('name') == 'masstrace_centroid_rt':
                data[:,4][n] = json.loads(root[0][n][m].attrib.get('value'))[0] # RT of m/z
    
    Df = pd.DataFrame(data, columns=['m/z', 'm/z+1', 'm/z intens',  'm/z+1 intens', 'RT'])
    
    return Df