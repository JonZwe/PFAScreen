# Function to perform mass trace detection, elution peak detection, feature detection, and adduct finding
# with the FeatureFinderMetabo algorithm from OpenMS (via PyOpenMS)
from pyopenms import (MSExperiment, MzMLFile, MassTraceDetection, ElutionPeakDetection, FeatureMap, 
                      FeatureFindingMetabo, FeatureXMLFile, MetaboliteFeatureDeconvolution, ConsensusMap)

# NOTE: Several parameters could be investigated in more detail

def MS1_feature_finding(sample, 
                        mass_error_ppm, 
                        noise_threshold_int,
                        isotope_model
                        ):
    
    exp = MSExperiment()
    MzMLFile().load(sample, exp)
    exp.sortSpectra(True)

    # 1)  Mass trace detection
    mass_traces = []
    mtd = MassTraceDetection()
    mtd_params = mtd.getDefaults()
    mtd_params.setValue("mass_error_ppm", float(mass_error_ppm))           # set according to your instrument mass error
    mtd_params.setValue("noise_threshold_int", float(noise_threshold_int)) # adjust to noise level in your data
    mtd.setParameters(mtd_params)
    mtd.run(exp, mass_traces, 0)

    # 2) Elution peak detection (deconvolution)
    mass_traces_split = []
    mass_traces_final = []
    epd = ElutionPeakDetection()
    epd_params = epd.getDefaults()

    # CHECK THE INFLUENCE OF THIS PARAMETER!!!
    epd_params.setValue("width_filtering", "fixed") # removes mass traces outside min_fwhm of 1 and max_fwhm of 60 s
    epd.setParameters(epd_params)
    epd.detectPeaks(mass_traces, mass_traces_split)

    if epd.getParameters().getValue("width_filtering") == "auto":
        epd.filterByPeakWidth(mass_traces_split, mass_traces_final)
    else:
        mass_traces_final = mass_traces_split

    # 3) Feature detection (isotope reduction)
    fm = FeatureMap()
    feat_chrom = [] # feature chromatograms
    ffm = FeatureFindingMetabo()
    ffm_params = ffm.getDefaults()
    ffm_params.setValue("isotope_filtering_model", isotope_model) # metabolites (2% RMS), metabolites (5% RMS)
    ffm_params.setValue("remove_single_traces", "true")           # set false to keep features with only one mass trace
    ffm_params.setValue("mz_scoring_by_elements", "false")
    ffm_params.setValue("report_convex_hulls", "false")
    ffm.setParameters(ffm_params)
    ffm.run(mass_traces_final, fm, feat_chrom)

    fm.setUniqueIds()
    # fm.setPrimaryMSRunPath(["ms_data.mzML".encode()])
    # fh = FeatureXMLFile()
    # fh.store("FM_" + sample_names[s]  + ".featureXML", fm)

    '''
    # 4) Adduct finding (optional)
    mfd = MetaboliteFeatureDeconvolution()
    mdf_par = mfd.getDefaults()

    # understand how adducts should be specified!
    ionization_mode = 'neg'
    if ionization_mode == 'pos':
        mdf_par.setValue("potential_adducts", [b"H:+:0.4",b"Na:+:0.2",b"NH4:+:0.2", b"H-1O-1:+:0.1", b"H-3O-2:+:0.1"])
    else:
        mdf_par.setValue("potential_adducts", [b"H-1:-:0.9", b"H-2O-1:0:0.05", b"CH2O2:0:0.5", b"CH3COOH:0:0.4", b"Cl:-0:0.1"])

    mdf_par.setValue("charge_min", 1, "Minimal possible charge")
    mdf_par.setValue("charge_max", 1, "Maximal possible charge")
    mdf_par.setValue("charge_span_max", 1)
    mdf_par.setValue("max_neutrals", 1)
    mfd.setParameters(mdf_par)
    feature_map_MFD = FeatureMap()
    cons_map0 = ConsensusMap()
    cons_map1 = ConsensusMap()
    mfd.compute(fm, feature_map_MFD, cons_map0, cons_map1)

    mfd = MetaboliteFeatureDeconvolution()
    mdf_par = mfd.getDefaults()
    mdf_par.setValue("potential_adducts", [b"H-1:-:0.9", b"H-2O-1:0:0.05", b"CH2O2:0:0.5", b"CH3COOH:0:0.4", b"Cl:-0:0.1"])
    mdf_par.setValue("charge_min", 1, "Minimal possible charge")
    mdf_par.setValue("charge_max", 1, "Maximal possible charge")
    mdf_par.setValue("charge_span_max", 1)
    mdf_par.setValue("max_neutrals", 1)
    mdf_par.setValue("retention_max_diff", 3.0)
    mdf_par.setValue("retention_max_diff_local", 3.0)
    mfd.setParameters(mdf_par)
    feature_map_MFD = FeatureMap()
    cons_map0 = ConsensusMap()
    cons_map1 = ConsensusMap()
    mfd.compute(fm, feature_map_MFD, cons_map0, cons_map1)

    '''
    return exp, fm