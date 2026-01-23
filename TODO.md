
# TODO for PFAScreen (2025/11/10)

## Create clear order
- Create a package!
- Sort in few modules that make sense (e.g., spectra methods, chromatogram methods, filtering methods, df_methods)
- Or: Methods applicable to raw data (e.g, oms.Experiment)
- Make screening.py clearer
    - Make its more reusable, exclude several functions that are inside
- Develop wrapper method which applies methods to alignment_df (like df.apply())
    - Ensure correct nan handling!
- Change all string concat to f-strings
- Give default parameters in all functions, and their data types + documentation
- Change naming convention to: pfascreen_df, samples_df
- add: if __name__ == __main__ statements (where needed)

## Issues to fix
- Make weird data types into dicts (e.g., dia_frags etc.), maybe also spectra?
- update get_adduct_data() for handling more complex adducts
- add alignment parameters
- add alignment possibility from openMS via individual feature lists to FeatureMap() (in feature_detection_comparis folder!)
- some isotopes patterns are still very wrong (is it the suspect list, happend for salts). Likely its wrong adduct assignment
- use networkx for a proper combine_dia_diffs.py
- Repair the wrong diagnostic fragments function (currently uses rounding!)
- dia/diff -> if other element -> error warning!

- PFAScreen 1: 
    - change openpyxl version to 3.1.0 in requirements.txt
    - why mz in Df_FeatureData instead of m/z?
    - blank_correc has an error if fold_change is not given!
    - remove multiple mzML option
    - BUG in line 57 in MS2 extractor?

## Additions: Would be nice
- Add spectral library matching option + metfrag
- Add CCS option (one simple function)
- Use Dash app (instead of tkinter) for GUI purposes & vizualization (static matplotlib, dynamic plotly)
- AG Grid: High-Performance React Grid, Angular Grid, JavaScript Grid (?)
- Analogue search (molecular networking with database, but no accurate mass match first)
- make KMD network, with KMD or m/C as color!
- in network plots: Allow direct accessibility of features
- mzmine & MSDIAL CLI to directly run from python?
- Add MS2 filtering options (md_cone, intensity, etc.)
- show all info in feature overview plot (MSMS, suspect hit, etc.)
- Most frequent diffs with m/C > 25!
- Suspect hits for whole HS group?
- Summary txt that summarizes the parameters