PFΔScreen
========

PFΔScreen is an open-source Python based non-target screening software tool to prioritize potential PFAS features in raw data from liquid- or gas chromatograpy coupled to high-resolution mass spectrometry (LC- or GC-HRMS) measurements with a simple graphical user interface (GUI). pyOpenMS (Python interface to the C++ OpenMS library) is used for feature detection in MS raw data. Optionally, custom feature lists can be included. PFΔScreen uses several techniques for prioritization such as the MD/C-m/C approach, Kendrick mass defect (KMD) analysis and fragment mass differences and diagnostic fragments in the MS2 data. PFΔScreen is easily installable via batch files. Raw mass spectrometric data can be included vendor-independently in the mzML format (data-dependent acquisition with centroided spectra, mzML files can be generated via the MSConvert software tool).

If you use PFΔScreen please cite (more detailed examples and explations can be found here):\
Zweigle, J.; Bugsel, B.; Fabregat-Palau, J., Zwiener, C., PFΔScreen – An open-source tool for automated PFAS feature prioritization in non-target HRMS data.

![alt text](https://github.com/JonZwe/PFAScreen/logo.jpg?raw=true)

Installation:
--------------------

PFΔScreen can be installed and executed within the standard Python environment. To make installation and use as easily as possible, PFΔScreen can be automatically installed with the Installation.bat file and executed with the Run_PFAScreen.bat file without the need of additional software. In the following, the two steps needed for a simple installation are explained.

### Automated installation

- Download PFΔSScreen by clicking on the green “Code” button and click “Download ZIP”. When downloaded, unzip the folder and move it to a local folder on your computer. 
-	Automatic installation of Python and the required packages with Installation.bat: Navigate into the folder where PFΔSScreen was copied (PFAScreen-main). Double click the Installation.bat file. The Windows command line interface will open, and the Microsoft Store opens automatically if you do not have Python installed on your computer. Click on “Install” and wait until the installation of Python is finished and close the Microsoft Store. Back to the Windows command line press any button to automatically install pip (Package Manager for Python) and in the following all required Python packages. Finally, when the message “Installation successfully finished” pops up, press any button and the installation is completed.

### Package installation 
If Python is already installed, the PFΔScreen dependencies can be installed within the standard Python environment with the following command:

```
pip install -r requirements.txt
```

General explanation
----------------------
To start PFΔScreen, double click the Run_PFAScreen.bat file. Both the GUI and a console window will open. To load a MS raw datafile, click the “Browse Sample.mzML” button and choose the mzML file of a sample and an optional mzML file of a blank control (Browse Blank.mzML). Sample and blank for raw data input in PFΔScreen should have been measured under data-dependent acquisition (ddMS2) with centroided spectra, ideally with one collision energy per precursor. The parameters for feature finding, MS2 alignment and blank correction can be specified and executed by pressing the “Run FeatureFinding” button. In case another feature finding procedure (e.g., from vendor software) is desired, custom feature lists (see external_feature_list.xlsx) together with the respective mzML files can instead be included in PFΔScreen. This is done by the “Browse SampleFeatures.xlsx“ and “Browse BlankFeatures.xlsx“ buttons, which are preprocessed by the “Run ExternalFeatureTable“ button. Note that data evaluation only works when the corresponding mzML files are also given; otherwise MS2 data would be missing. Whenever the FeatureFinding tab is completed, the RawDataVisualization can be used even without PFAS-specific data. 

To perform the PFASPrioritization, appropriate input parameters can be set, and then PFAS-specific data evaluation is performed with the “Run PFASPrioritization“ button. The overall rather short runtime (e.g., less than one minute in case of 4000 spectra per sample for the whole workflow), allows a convenient adjustment of input parameters. Afterwards, MS2 spectra displayed by the RawDataVisualization tool (MS2 extractor), have highlighted fragment mass differences and diagnostic fragments, if some were detected. After executing the PFASPrioritization tab, the PFΔScreen results table (Excel format) and several interactive HTML plots are saved in a folder named after the sample that can be easily inspected, including a MD/C-m/C plot, a m/z vs. RT plot (with and without MS2 raw data), a KMD vs. m/z with linked m/z vs. RT plot (to verify systematic RT-shifts), and a m/C histogram. Data from the results table can be used to visualize EICs (and extrapolate HS with common repeating units such as CF2), MS1 and MS2 spectra. Additionally, a coelution correlation can be performed with the RawDataVisualization tool. Also, the theoretical isotope patterns of suspect hits can be displayed over the experimental isotope patterns (MS1).

Call for Contributions
----------------------

We appreciate help and suggestions from other researchers to improve
PFΔScreen. Please don't hesitate to contact us.
