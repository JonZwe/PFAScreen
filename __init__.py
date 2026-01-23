"""
PFAScreen - PFAS screening and analysis package
"""

__version__ = "1.0.0"
__author__ = "Your Name"

# Import main functions for easy access
from .screening import pfascreen  # renamed from pfascreen.py to screening.py
from .params import get_config
from .run_oms_pipeline import run_oms_pipeline
from .pfascreen_utils.utils import filter_df, fold_change_filter
from .pfascreen_utils.plotting import (
    kmd_plot, scatter_plot, feature_overview_plot
)
from .pfascreen_utils.generate_html import generate_full_html
from .networks import molecular_network, mass_difference_network
from .readers_writers import msdial_alignmenttable_to_df

__all__ = [
    'pfascreen',
    'get_config', 
    'run_oms_pipeline',
    'filter_df',
    'fold_change_filter',
    'kmd_plot',
    'scatter_plot',
    'feature_overview_plot',
    'generate_full_html',
    'molecular_network',
    'mass_difference_network',
    'msdial_alignmenttable_to_df'
]
