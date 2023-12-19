# %%
# PFAScreen: open-source tool to perform PFAS-specific
# non-targeted screening on LC- and GC-HRMS raw data
import os
import tkinter as tk
from tkinter import ttk, Button, IntVar, BooleanVar, Radiobutton, filedialog, font
from PIL import ImageTk, Image
from PFAS_feature_prioritization import PFAS_feature_prioritization
from feature_preproccessing import feature_preproccessing
from external_feature_preproccessing import external_feature_preproccessing
from EIC_extractor import EIC_extractor
from MS1_extractor import MS1_extractor
from MS2_extractor import MS2_extractor
from EIC_correlator import EIC_correlator

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

def browse_button_sample_mzML():
    global sample_path_var 
    sample_path = filedialog.askopenfilename()
    sample_path_var = r'{}'.format(sample_path)
    label_sample_path.configure(text = os.path.basename(sample_path_var)[0:20])

def browse_button_blank_mzML():
    global blank_path_var
    blank_path = filedialog.askopenfilename()
    blank_path_var = r'{}'.format(blank_path)
    label_blank_path.configure(text = os.path.basename(blank_path_var)[0:20])

def browse_button_sample_features_csv():
    global sample_path_csv_var
    sample_path_csv = filedialog.askopenfilename()
    sample_path_csv_var = r'{}'.format(sample_path_csv)
    label_sample_external_path.configure(text = os.path.basename(sample_path_csv)[0:20])

def browse_button_blank_features_csv():
    global blank_path_csv_var
    blank_path_csv = filedialog.askopenfilename()
    blank_path_csv_var = r'{}'.format(blank_path_csv)
    label_blank_external_path.configure(text = os.path.basename(blank_path_csv)[0:20])

if not 'blank_path_var' in globals():
    blank_path_var = ''

def RunFeatureFinding():
    global Df_FeatureData, Df_MS2RawData, idx_in_features, exps
    Df_FeatureData, Df_MS2RawData, idx_in_features, exps = feature_preproccessing(
                                                                sample_path_var, 
                                                                blank_path_var,
                                                                mass_error_ppm = float(entryMassError.get()),
                                                                intensity_threshold = float(entryIntThreshold.get()),
                                                                isotope_model_bool = var_IsotopeModel.get(),
                                                                mass_tolerance_MSMS_assignment = float(entryMSMSMassTol.get()),
                                                                RT_tolerance_MSMS_assignment = float(entryMSMSRTTol.get()),
                                                                mass_tolerance_blank_correction = float(entryBlankMassTol.get()),
                                                                RT_tolerance_blank_correction = float(entryBlankRTTol.get()),
                                                                fold_change_blank_correction = float(entryBlankFoldChange.get())
                                                                )
if not 'blank_path_csv_var' in globals():
    blank_path_csv_var = ''    

def RunExternalFeatures():
    global Df_FeatureData, Df_MS2RawData, idx_in_features
    Df_FeatureData, Df_MS2RawData, idx_in_features = external_feature_preproccessing(
                                                        sample_path_var,
                                                        sample_path_csv_var, 
                                                        blank_path_csv_var,
                                                        mass_tolerance_MSMS_assignment = float(entryMSMSMassTol.get()),
                                                        RT_tolerance_MSMS_assignment = float(entryMSMSRTTol.get()),
                                                        mass_tolerance_blank_correction = float(entryBlankMassTol.get()),
                                                        RT_tolerance_blank_correction = float(entryBlankRTTol.get()),
                                                        fold_change_blank_correction = float(entryBlankFoldChange.get())
                                                        )
def RunPFASPrioritization():
    global Df_FeatureData_Final, Df_FindPFAS
    Df_FeatureData_Final, Df_FindPFAS = PFAS_feature_prioritization(
                                            Df_FeatureData,
                                            Df_MS2RawData,
                                            idx_in_features,
                                            Results_folder = f'PFAScreen_{os.path.basename(sample_path_var).rsplit(".", 1)[0]}',
                                            diffs = entryMassDiffs.get().split(' '),
                                            number_of_fragments = float(entryNFrags.get()),
                                            mass_tolerance_FindPFAS = float(entryMassTolFrags.get()),
                                            intensity_threshold = float(entryIntThresholdMSMS.get()),
                                            diffs_KMD = entryKMDDiffs.get().split(' '),
                                            mass_tolerance_KMD = float(entryKMDMassTol.get()),
                                            n_homologues = float(entryNumberOfHS.get()),
                                            adducts = int(var_AdductsSuspectScreen.get()),
                                            tol_suspect = float(entrySuspectMassTol.get()),
                                            save_MSMS_spectra = var_save_MSMS_spectra.get(),
                                            mC_range = entrymCCutoff.get().split(' '),
                                            MDC_range = entryMDCCutoff.get().split(' '),
                                            MD_range = entryMDRange.get().split(' ')
                                            )

def RunEIC():
    EIC_extractor(
        exps,
        mz_list = entryEICMass.get().split(','),
        extraction_width = float(entryEICwidth.get()),
        rep_unit = str(entryHSrepunit.get()),
        n_rep_unit = float(entryHSnrepunit.get())
        )
    
def RunMS1_extractor():
    MS1_extractor(
        exps[0],
        RT = float(entryMS1RT.get()),
        formula = str(entryMS1Formula.get())
        )
    
if not 'Df_FindPFAS' in globals():
    Df_FindPFAS = '' # check if FindPFAS is there

def RunMS2_extractor():
    MS2_extractor(
        Df_MS2RawData,
        Df_FindPFAS,
        prec_mz = float(entryprecMS2.get())
        )
    
def RunEIC_correlator():
    EIC_correlator(
        exps[0],
        Df_FeatureData,
        mz_interest = float(entryEIC_mz_interest.get()),
        extraction_width = float(entryEICwidth_corr.get()),
        RT_interest = float(entryEICRT_corr.get()),
        RT_width = float(entryEICRTwidth_corr.get()),
        R2_threshold = float(entryR2_threshold.get())
        )
    
# %% ==============================================================
# DEFAULT PARAMETERS
# FeatureFinding
mass_error_ppm = 10
intensity_threshold_int = 1000

mass_tolerance_MSMS_assignment = 0.005 # HIGHLY IMPORTANT
RT_tolerance_MSMS_assignment = 0.2*60  # Consider peak width

mass_tolerance_blank_correction = 0.002
RT_tolerance_blank_correction = 0.1*60
fold_change_blank_correction = 5

# MD/C-m/C
mC_range = [0, float('inf')]
MDC_range = [-0.5, 0.5]
MD_range = [-0.5, 0.5]

# KMD analysis
diffs_KMD = ['CF2']
mass_tolerance_KMD = 0.002
n_homologues = 3

# MS2 differences & diagnostic fragments
diffs = ['CF2', 'C2F4', 'HF']
number_of_fragments = 1 # factor of 2 for diagnostic fragments
mass_tolerance_FindPFAS = 0.002
intensity_threshold = 1000
plot_MS2_spectra = True
save_HTML = True

# Suspect Screening
ionization_mode = 'neg' # specify ionization mode
tol_suspect = 0.002

# EIC extractor
mz_list = []
extraction_width = 0.005
rep_unit = ''
n_rep_unit = 3

# MS1 spectrum extractor
RT = ''
formula = ''

# MS2 spectrum extractor
prec_mz = ''

# EIC correlator
mz_interest = ''
RT_interest = ''
RT_width = 20
R2_threshold = 0.95



# %% ==============================================================
# Create GUI via tkinter
root = tk.Tk(className='PF\u0394Screen')
root.title('PF\u0394Screen')
root.geometry('1000x600')
root.resizable(width=True, height=True)
style = ttk.Style(root)

image = Image.open(os.path.join(script_dir, 'logo.jpg')).resize((143, 90))
photo = ImageTk.PhotoImage(image)
labelLogo = ttk.Label(root, image=photo)
labelLogo.place(x = 820, y = 20)
labelLogo.configure(font=('Courier', 18, 'bold'))

s = ttk.Style()
s.configure('Red.TLabelframe.Label', font=('Courier', 13, 'bold'))
label_frame = ttk.LabelFrame(root, text='FeatureFinding', style = "Red.TLabelframe")
label_frame.place(x=30, y=10, width=300, height=580)

label_frame = ttk.LabelFrame(root, text='PFASPrioritization', style = "Red.TLabelframe")
label_frame.place(x=350, y=10, width=300, height=580)

label_frame = ttk.LabelFrame(root, text='RawDataVisualization', style = "Red.TLabelframe")
label_frame.place(x=710, y=160, width=280, height=380)

label_sample_path = ttk.Label(root)
label_sample_path.place(x = 200, y = 40) 

label_blank_path = ttk.Label(root)
label_blank_path.place(x = 200, y = 70)

label_sample_external_path = ttk.Label(root)
label_sample_external_path.place(x = 200, y = 100) 

label_blank_external_path = ttk.Label(root)
label_blank_external_path.place(x = 200, y = 120) 


# ==============================================================
myfont = font.Font(size=10, weight='bold')
ButtonSamplePath = Button(text="Browse Sample.mzML", bg='#E7ECF7', command = browse_button_sample_mzML)
ButtonSamplePath.place(x = 40, y = 40, height = 30, width = 160)
ButtonSamplePath['font'] = myfont

ButtonBlankPath = Button(text="Browse Blank.mzML", bg='#E7ECF7', command = browse_button_blank_mzML)
ButtonBlankPath.place(x = 40, y = 70, height = 30, width = 160)
ButtonBlankPath['font'] = myfont

ButtonSamplePathcsv = Button(text="Browse SampleFeatures.csv", command = browse_button_sample_features_csv)
ButtonSamplePathcsv.place(x = 40, y = 100, height = 20, width = 160)

ButtonBlankPathcsv = Button(text="Browse BlankFeatures.csv", command = browse_button_blank_features_csv)
ButtonBlankPathcsv.place(x = 40, y = 120, height = 20, width = 160)


labelTitleMSMS = ttk.Label(root,  text='Peak finding')
labelTitleMSMS.place(x = 100, y = 150)
labelTitleMSMS.configure(font=('Courier', 13))

labelMassError = ttk.Label(root, text='Mass error (ppm)')
labelMassError.place(x = 50, y = 180) 
entryMassError = ttk.Entry(root)
entryMassError.place(x = 200, y = 180, height = 20, width = 50)
entryMassError.insert('end', mass_error_ppm)

labelIntThreshold = ttk.Label(root, text='Intensity threshold')
labelIntThreshold.place(x = 50, y = 200) 
entryIntThreshold = ttk.Entry(root)
entryIntThreshold.place(x = 200, y = 200, height = 20, width = 50)
entryIntThreshold.insert('end', intensity_threshold_int)


var_IsotopeModel = BooleanVar()
var_IsotopeModel.set(True)
CheckButtonIsotopeModel = tk.Checkbutton(root, text='Isotope model',variable = var_IsotopeModel, onvalue=1, offvalue=0)
CheckButtonIsotopeModel.place(x = 50, y = 220)


labelTitleMSMS = ttk.Label(root,  text='MS/MS alignment')
labelTitleMSMS.place(x = 100, y = 250)
labelTitleMSMS.configure(font=('Courier', 13))

labelMSMSMassTol = ttk.Label(root, text='Mass tolerance (Da)')
labelMSMSMassTol.place(x = 50, y = 280)
entryMSMSMassTol = ttk.Entry(root)
entryMSMSMassTol.place(x = 200, y = 280, height = 20, width = 50)
entryMSMSMassTol.insert('end', mass_tolerance_MSMS_assignment)

labelMSMSRTTol = ttk.Label(root, text='RT tolerance (s)')
labelMSMSRTTol.place(x = 50, y = 310)
entryMSMSRTTol = ttk.Entry(root)
entryMSMSRTTol.place(x = 200, y = 310, height = 20, width = 50)
entryMSMSRTTol.insert('end', RT_tolerance_MSMS_assignment)


labelTitleBlankCorr = ttk.Label(root,  text='Blank correction')
labelTitleBlankCorr.place(x = 100, y = 340) 
labelTitleBlankCorr.configure(font=('Courier', 13))

labelBlankMassTol = ttk.Label(root, text='Mass tolerance (Da)')
labelBlankMassTol.place(x = 50, y = 370)
entryBlankMassTol = ttk.Entry(root)
entryBlankMassTol.place(x = 200, y = 370, height = 20, width = 50)
entryBlankMassTol.insert('end', mass_tolerance_blank_correction)

labelBlankRTTol = ttk.Label(root, text='RT tolerance (s)')
labelBlankRTTol.place(x = 50, y = 400)
entryBlankRTTol = ttk.Entry(root)
entryBlankRTTol.place(x = 200, y = 400, height = 20, width = 50)
entryBlankRTTol.insert('end', RT_tolerance_blank_correction)

labelBlankFoldChange = ttk.Label(root, text='Fold change')
labelBlankFoldChange.place(x = 50, y = 430)
entryBlankFoldChange = ttk.Entry(root)
entryBlankFoldChange.place(x = 200, y = 430, height = 20, width = 50)
entryBlankFoldChange.insert('end', fold_change_blank_correction)

buttonRunFeatureFinding = Button(root, text = 'Run FeatureFinding', bg='#E7ECF7', command = RunFeatureFinding)
buttonRunFeatureFinding.place(x = 100, y = 490, height = 50, width = 150)

buttonRunExternalFeatures = Button(root, text = 'Run ExternalFeatureTable', command = RunExternalFeatures)
buttonRunExternalFeatures.place(x = 100, y = 550, height = 25, width = 150)


# PFAS feature prioritization
# ==============================================================
labelTitleMDCmC = ttk.Label(root,  text='MD/C-m/C & MD filtering')
labelTitleMDCmC.place(x = 400, y = 40) 
labelTitleMDCmC.configure(font=('Courier', 13))

labelmCCutoff = ttk.Label(root, text='m/C range')
labelmCCutoff.place(x = 370, y = 65) 
entrymCCutoff = ttk.Entry(root)
entrymCCutoff.place(x = 520, y = 65, height = 20, width = 50)
entrymCCutoff.insert('end', mC_range)

labelMDCCutoff = ttk.Label(root, text='MD/C range')
labelMDCCutoff.place(x = 370, y = 85) 
entryMDCCutoff = ttk.Entry(root)
entryMDCCutoff.place(x = 520, y = 85, height = 20, width = 50)
entryMDCCutoff.insert('end', MDC_range)

labelMDRange = ttk.Label(root, text='MD range')
labelMDRange.place(x = 370, y = 105) 
entryMDRange = ttk.Entry(root)
entryMDRange.place(x = 520, y = 105, height = 20, width = 50)
entryMDRange.insert('end', MD_range)


labelTitleKMD = ttk.Label(root,  text='Kendrick mass defect')
labelTitleKMD.place(x = 400, y = 125) 
labelTitleKMD.configure(font=('Courier', 13))

labelKMDDiffs = ttk.Label(root, text='KMD difference')
labelKMDDiffs.place(x = 370, y = 150) 
entryKMDDiffs = ttk.Entry(root)
entryKMDDiffs.place(x = 520, y = 150, height = 20, width = 50)
entryKMDDiffs.insert('end', diffs_KMD)

labelKMDMassTol = ttk.Label(root, text='KMD mass tolerance (Da)')
labelKMDMassTol.place(x = 370, y = 170) 
entryKMDMassTol = ttk.Entry(root)
entryKMDMassTol.place(x = 520, y = 170, height = 20, width = 50)
entryKMDMassTol.insert('end', mass_tolerance_KMD)

labelNumberOfHS = ttk.Label(root, text='Number of homologues')
labelNumberOfHS.place(x = 370, y = 190) 
entryNumberOfHS = ttk.Entry(root)
entryNumberOfHS.place(x = 520, y = 190, height = 20, width = 50)
entryNumberOfHS.insert('end', n_homologues)



labelTitleMS2 = ttk.Label(root,  text='MS2 differences and DFs')
labelTitleMS2.place(x = 400, y = 220) 
labelTitleMS2.configure(font=('Courier', 13))

labelMassDiffs = ttk.Label(root, text='Mass differences')
labelMassDiffs.place(x = 370, y = 250) 
entryMassDiffs = ttk.Entry(root)
entryMassDiffs.place(x = 520, y = 250, height = 20, width = 100)
entryMassDiffs.insert('end', diffs)

labelNFrags = ttk.Label(root, text='Number of fragments')
labelNFrags.place(x = 370, y = 280) 
entryNFrags = ttk.Entry(root)
entryNFrags.place(x = 520, y = 280, height = 20, width = 50)
entryNFrags.insert('end', number_of_fragments)

labelMassTolFrags = ttk.Label(root, text='Mass tolerance Frags (Da)')
labelMassTolFrags.place(x = 370, y = 310) 
entryMassTolFrags = ttk.Entry(root)
entryMassTolFrags.place(x = 520, y = 310, height = 20, width = 50)
entryMassTolFrags.insert('end', mass_tolerance_FindPFAS)

labelIntThresholdMSMS = ttk.Label(root, text='Intensity threshold MS/MS')
labelIntThresholdMSMS.place(x = 370, y = 340) 
entryIntThresholdMSMS = ttk.Entry(root)
entryIntThresholdMSMS.place(x = 520, y = 340, height = 20, width = 50)
entryIntThresholdMSMS.insert('end', intensity_threshold)

var_save_MSMS_spectra = BooleanVar()
SaveMSMS_true = Radiobutton(root, text="yes", variable = var_save_MSMS_spectra, value = True)
SaveMSMS_true.place(x = 370, y = 360)
SaveMSMS_false = Radiobutton(root, text="no", variable = var_save_MSMS_spectra, value = False)
SaveMSMS_false.place(x = 420, y = 360)
var_save_MSMS_spectra.set(False)

labelHTML = ttk.Label(root, text='Save HTML MSMS spectra?')
labelHTML.place(x = 470, y =363) 



labelTitleSuspect = ttk.Label(root, text='Suspect screening')
labelTitleSuspect.place(x = 400, y = 390) 
labelTitleSuspect.configure(font=('Courier', 13))

var_AdductsSuspectScreen = IntVar()
M_min_H = Radiobutton(root, text="M-H", variable = var_AdductsSuspectScreen, value = 1)
M_min_H.place(x = 370, y = 420)
M_plus_H = Radiobutton(root, text="M+H", variable = var_AdductsSuspectScreen, value = 2)
M_plus_H.place(x = 450, y = 420)
M_plus = Radiobutton(root, text="M+", variable = var_AdductsSuspectScreen, value = 3)
M_plus.place(x = 530, y = 420)
var_AdductsSuspectScreen.set(1)

labelSuspectMassTol = ttk.Label(root, text='Suspect Mass Tol (Da)')
labelSuspectMassTol.place(x = 370, y = 450) 
entrySuspectMassTol = ttk.Entry(root)
entrySuspectMassTol.place(x = 520, y = 450, height = 20, width = 50)
entrySuspectMassTol.insert('end', tol_suspect)


buttonRunPFASPrioritization = Button(root, text = 'Run PFASPrioritization', bg = '#F7EFEA', command = RunPFASPrioritization)
buttonRunPFASPrioritization.place(x = 400, y = 520, height = 50, width = 150)

# Raw Data Visualization
# ==============================================================

buttonRunEIC = Button(root, text = 'EICExtractor', bg='#F1EAF7', command = RunEIC)
buttonRunEIC.place(x = 880, y = 235, height = 40, width = 80)

labelEICMass = ttk.Label(root, text='EIC m/z')
labelEICMass.place(x = 730, y = 210) 
entryEICMass = ttk.Entry(root)
entryEICMass.place(x = 810, y = 210, height = 20, width = 150)
entryEICMass.insert('end', mz_list)

labelEICwidth = ttk.Label(root, text='EIC width (Da)')
labelEICwidth.place(x = 730, y = 230) 
entryEICwidth = ttk.Entry(root)
entryEICwidth.place(x = 810, y = 230, height = 20, width = 50)
entryEICwidth.insert('end', extraction_width)

labelHSrepunit = ttk.Label(root, text='Rep. unit')
labelHSrepunit.place(x = 730, y = 250) 
entryHSrepunit = ttk.Entry(root)
entryHSrepunit.place(x = 810, y = 250, height = 20, width = 50)
entryHSrepunit.insert('end', rep_unit)

labelHSnrepunit = ttk.Label(root, text='n')
labelHSnrepunit.place(x = 730, y = 270) 
entryHSnrepunit = ttk.Entry(root)
entryHSnrepunit.place(x = 810, y = 270, height = 20, width = 50)
entryHSnrepunit.insert('end', n_rep_unit)

# ==============================================================
buttonRunMS1_extractor = Button(root, text = 'MS1Extractor', bg='#F1EAF7', command = RunMS1_extractor)
buttonRunMS1_extractor.place(x = 880, y = 310, height = 40, width = 80)

labelMS1RT = ttk.Label(root, text='RT (min)')
labelMS1RT.place(x = 730, y = 310)
entryMS1RT = ttk.Entry(root)
entryMS1RT.place(x = 780, y = 310, height = 20, width = 90)
entryMS1RT.insert('end', RT)

labelMS1Formula = ttk.Label(root, text='Formula')
labelMS1Formula.place(x = 730, y = 330) 
entryMS1Formula = ttk.Entry(root)
entryMS1Formula.place(x = 780, y = 330, height = 20, width = 90)
entryMS1Formula.insert('end', formula)

# ==============================================================
buttonRunMS2_extractor = Button(root, text = 'MS2Extractor', bg='#F1EAF7', command = RunMS2_extractor)
buttonRunMS2_extractor.place(x = 880, y = 380, height = 30, width = 80)

labelprecMS2 = ttk.Label(root, text='Prec m/z')
labelprecMS2.place(x = 730, y = 380) 
entryprecMS2 = ttk.Entry(root)
entryprecMS2.place(x = 780, y = 380, height = 20, width = 90)
entryprecMS2.insert('end', prec_mz)

# ==============================================================
buttonRunEIC_correlator = Button(root, text = 'EICCorrelator', bg='#F1EAF7', command = RunEIC_correlator)
buttonRunEIC_correlator.place(x = 880, y = 455, height = 40, width = 80)

labelEIC_mz_interest = ttk.Label(root, text='m/z')
labelEIC_mz_interest.place(x = 730, y = 430) 
entryEIC_mz_interest = ttk.Entry(root)
entryEIC_mz_interest.place(x = 810, y = 430, height = 20, width = 150)
entryEIC_mz_interest.insert('end', mz_interest)

labelEICwidth_corr = ttk.Label(root, text='EIC width (Da)')
labelEICwidth_corr.place(x = 730, y = 450) 
entryEICwidth_corr = ttk.Entry(root)
entryEICwidth_corr.place(x = 810, y = 450, height = 20, width = 50)
entryEICwidth_corr.insert('end', extraction_width)

labelEICRT_corr = ttk.Label(root, text='RT (min)')
labelEICRT_corr.place(x = 730, y = 470) 
entryEICRT_corr = ttk.Entry(root)
entryEICRT_corr.place(x = 810, y = 470, height = 20, width = 50)
entryEICRT_corr.insert('end', RT)

labelEICRTwidth_corr = ttk.Label(root, text='RT width (s)')
labelEICRTwidth_corr.place(x = 730, y = 490) 
entryEICRTwidth_corr = ttk.Entry(root)
entryEICRTwidth_corr.place(x = 810, y = 490, height = 20, width = 50)
entryEICRTwidth_corr.insert('end', RT_width)

labelR2_threshold = ttk.Label(root, text='R2 threshold')
labelR2_threshold.place(x = 730, y = 510) 
entryR2_threshold = ttk.Entry(root)
entryR2_threshold.place(x = 810, y = 510, height = 20, width = 50)
entryR2_threshold.insert('end', R2_threshold)


root.mainloop()
