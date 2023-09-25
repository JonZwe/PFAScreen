# Function to perform mass trace detection, elution peak detection, feature detection, and adduct finding
# with the FeatureFinderMetabo algorithm from OpenMS (via PyOpenMS)
from pyopenms import (MSExperiment, MzMLFile, MassTraceDetection, ElutionPeakDetection, FeatureMap, 
                      FeatureFindingMetabo, FeatureXMLFile, MetaboliteFeatureDeconvolution, ConsensusMap)

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

    return exp, fm
