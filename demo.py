import os
import matplotlib.pyplot as plt
from run_oms_pipeline import run_oms_pipeline
from screening import pfascreen
from params import get_config
from utils import filter_df, fold_change_filter
from visualization import feature_overview_plot, scatter_plot, kmd_plot
from generate_html import generate_full_html
from readers_writers import msdial_alignmenttable_to_df
from networks import molecular_network, mass_difference_network

"""Demo function for PFAScreen workflow"""
# Load configuration
config, sample_names, output_folder = get_config()

# Run complete OpenMS pipeline (feature finding, alignment, isotopes, and MS2 data alignment)
df = run_oms_pipeline(config, sample_names)
#df = msdial_alignmenttable_to_df(r"D:\MS_raw_data\DeltaPFAS\UFZ_RPLC_Orbitrap_neg\mzML\Area_1_2025_08_31_15_51_03.txt")

# Perform PFAScreen analysis
df['adduct'] = '[M-H]-' if config['polarity'] == 'neg' else '[M+H]+'

df_pfascreen = pfascreen(df,
                         sample_names,
                         config,
                         adducts=1
                         )
# df_pfascreen = pd.read_pickle('df_pfascreen_demo.pkl')

# perform blank filtering (fold change 5)
df_blank_fil = fold_change_filter(df_pfascreen, sample_names, 'Blank', 3)

# Filter for specific m/C range
df_fil = filter_df(df_blank_fil, 'm/C', 25, 300)
df_fil = filter_df(df_fil, 'MD', -0.25, 0.1)


scatter_plot(df_fil, 'm/C', 'MD/C', col='KMD')
kmd_plot(df_fil, 'CF2', hs_tol=0.005, n_min=4)
feature_overview_plot(df_fil, 7, sample_names)

molecular_network(df_fil, 
                  col_highlight='compound_names',
                  colormap_column='rt',
                  score_cutoff=0.8, 
                  max_links=15)

mass_difference_network(df_fil['mz'].values, 
                        diffs=['C2F4', 'O', 'O2'], 
                        mz_tol=0.005, 
                        network_n=1, 
                        k=0.1)

# Generate HTML report
plt.switch_backend('Agg')
generate_full_html(df_fil, sample_names, os.path.join(output_folder, "Results_new.html"))
plt.switch_backend('TkAgg')