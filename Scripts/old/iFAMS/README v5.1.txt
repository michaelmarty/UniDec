Welcome to iFAMS (interactive Fourier Analysis for Mass Spectrometry) version 5.1!
(by Sean P. Cleary and James S. Prell, Department of Chemistry and Biochemistry, 1253 University of Oregon, Eugene, OR, 97403-1253)

If you use this program for an article or presentation, please cite the following article:

Cleary, Sean P.; Thompson, Avery M.; Prell, James S. "Fourier Analysis Method for Analyzing Highly Congested Mass Spectra of Ion Populations with Repeated Subunits" Analytical Chemistry, 2016, 88, 6205-6213, doi: 10.1021/acs.analchem.6b01088.

Please contact Sean Cleary at scleary@uoregon.edu with any questions, and check our website (prellgroup.uoregon.edu) for future updates of IFAMS.


This program is designed for use with Python v. 3.5, using the PyQt5 interface, and the PyQtGraph (included in the zip folder) module for all of the plotting.
To ensure that this is installed in your site packages, you can use the PIP installer (if you using something like pycharm to compile your python code) or the conda installer (specifically for for Anaconda).  To install, open a command prompt (for anaconda, open your conda command prompt) and simply type "pip install pyqt5" (on anaconda, this will be "conda install pyqt5" although you may already have it if you are using a later version of Anaconda.) 

After you have downloaded and installed Anaconda, unzip the file IFAMSv4.1.zip, and click IFAMSv4-1.py.


Please note: all input files for iFAMS version 5.1 must be comma-separated values (.csv) or tab-separated (.txt) files!


Included in the zip file is a data file "test.csv", which contains a Nanodisc mass spectrum. Please go through the following steps to familiarize yourself with the program.

1.) After you open IFAMSv5-1.py, you should see a graphical user interface with a small white rectangle in the middle.
2.) Load "test.csv" by clicking "file", then "Load Spectrum". The mass spectrum should load in the top left corner.  
3.) Under "Analysis", click "Fourier Transform". To the right of the mass spectrum, after a couple seconds, the (absolute value of the) Fourier transform should now appear. The x-axis will be in units of z/m. You can zoom in by using you mouse (please note, if you are using a laptop, I find it easier to first right click on the graph itself, go to mouse options, and click "1 button")
4.) You can now export these data as a CSV file by clicking "Export Fourier Spectrum" under "File".
5.) There should be small red dots at the tops of the peaks in the Fourier spectrum. iFAMS is set up to find peaks based on local maxima in the spectrum. The defaults for finding these peak maxima in iFAMS have been chosen such that they are a realistic distance apart for subunits weighing a few hundred Da.
6.) The peaks (representing fundamental frequencies for the charge states present in the Nanodisc mass spectrum) that we are looking for have been identified, but other peaks before and after those peaks have also been found.  We will want to correct that before calculating charge states and sub-unit mass, as the program calculates these values starting with the first peak found in the spectrum.
7.) In order to correct the peak maxima, we will adjust original preset values that the program uses when it finds peaks.  These values can be found in the first column on the left.
8.) In "minimum frequency used", type "0.01". This value tells iFAMS where to begin looking for a local maximum, where the value is associted with frequency.   
9.) In "min FFT peak height %", type "25". This value, the noise threshold, adjusts the minimum % of the tallest peak in the Fourier spectrum window that will be identified as a peak (see below for more information).
10) In "Delta", type "5".  This value adjusts the minimum distance between acceptable local maxima when searching for peaks.  See below for more information.
11.) You may have noticed that "5" is already the default value for "Delta".  Even though this is the default value, you will still need to type a value for all three settings in order to re-find peaks. See note below for more information.
12.) Click "Update these values". You will now see the dots turn blue, and the peaks that we are interested in are now the first peaks in the Fourier spectra. In future spectra, you may need to play with the previous three values to select the peaks you want.  
13.) Under "Analysis", Click "Calculate Charge States and Average Sub-unit Mass". Under Status/information, you should now see "Sub-unit Mass: 677.177 +- 3.516", "Charge States/Rounded From: [12, 13, 14, 15], [12.056603773584914, 13.075471698113217, ...]".
14.) The above information means that iFAMS has identified a repeated subunit with mass ~677 Da (the correct mass is actually 677.9 Da), and that the dominant charge states in the mass spectrum appear to be 12+ through 15+.
15.) Note that the peaks identified as fundamentals in the Fourier spectrum for each charge state now have large green dots. These are the peaks that have been used for the subunit mass and charge state determination.
16) I have found in some spectra that this process of finding peaks can sometimes be difficult to do accurately.  So I have also included in this version of iFAMS a way to manually input the charge states and subunit mass if you know what they are.  To do this, under "analysis" click "manually enter charge states and subunit mass". You will notice that the tabs to input the paramters to find peak maxima have now been changed to "minimum charge state, Maximum charge state, and subunit mass". In this case, you will be enter range of charge states (for example, if you want 12-15+, you would enter 12 into "Minimum charge state" and 15 in "Maximum charge state") along with the subunit mass. Please note, iFAMS will not fit exactly what you enter, but rather find points in the Fourier spectrum that most closely fit what you enter, as the Fourier spectrum is a descrete data set. Often times, this will not be the tops of peaks!) 
17.) We will now attempt to reconstruct the mass spectral envelope functions of each charge state, which include white-noise bands.  First, we need to calculate the noise associated with these functions, which is done using the second column.
18.) In "Minimum Noise Frequency", type "0.2".  In "Maximum Noise Frequency", type "0.3". This represents the range of frequencies over which iFAMS calculates the noise.  In general, you should use a range where you believe there is only white noise.
19.) Click "Calculate Noise Statistics".  You will see values for the average value for both the real and imaginary noise, as well as the standard deviation of each under "Status/information".
20.) Now, under "Analysis", click on "Reconstruct Envelope Functions". Overlaid on top of the mass spectrum, you should now see 4 new spectral envelopes in different colors, indexed by their corresponding charge state, as well as the noise associated with that function, represented as dotted lines.
21.) You can export these data by clicking "Export Envelope Functions".   This will export both the envelope function and the noise as CSV file, with each charge state exported as a separate, labeled file.
22.) Two special notes here.  First, the "Recontruct Envelope Functions" button will not work until after white noise statistics have been calculated.  Second, to increase the speed of the program, the number of data points for the noise is reduced to 1/50 of the points of the envelope functions/Mass spectrum.
23.) You can use iFAMS to analyze the reconstructed envelope functions to determine the average and standard deviation in the number of sub-units for charge state.  To do this, press "Calculate Sub-Unit Statistics" under "Analysis".
24.) A new Window should open up with one column and three variables. In "Input Charge", type "12".  This will be the charge state you wish to analyze.
25.) In "Base Mass", type "49323.8". This value represents the "base" mass, i.e., the mass in the absence of any sub-units.  For example, if you are studying a protein-lipid complex, where the number of sub-units would be the amount of lipids on the protein, this would be the total protein mass.
26.) In "Sub-unit Mass", type "677.177".  This is the value of sub-unit mass you calculated ealier, but you can adjust this as needed.  
27.) Click "Analyze Envelope Function".  This will prompt an "Open File" window.  Locate the appropriate envelope function data, stored where your original data is located.  In this case, pick the file titled "testIFFT12.csv"
28.) iFAMS will now calculate the sub-unit statistics for you.  In the text window at the bottom window, the values should appear.  The text should now read, "the average total sub-units is 124.4...", and "with a stdev of 19.77..."
29.) This version of iFAMS also includes the ability to Fourier Filter your spectra.  To do this, under analysis, type "Fourier Filter"  Please note, in order to use the Fourier Filter, you need to have first calculated the charge states and subunit mass. iFAMS uses these to determine where to apply the filter.
30.) A new window will pop up with 2 input parameters.  The first, "number of overtones used" is how many series of overtone you wish the filter to use.  The second, "zero-frequency data" is how many "windows" you would like to go out when selecting the zero frequency data.  The length of this window is determined by the fundamental frequency, where, for example, if typed 3, the window would be 3/subunitmass.  Type in 3 for each of these and click "calculate"
31.) You should now see three new graphs overlaying your mass spectrum. The first, in red, is the Fourier filtered data.  The second, which is in blue, is the "Low information signal" or basically what was selected using the zero frequency data.  The third, baseline subtracted data, is the mass spectrum minus the low information signal  
 

Now try these same steps with a .csv or .txt mass spectrum file of your own.


The user has some flexibility in controlling the way that IFAMS analyzes the Fourier spectrum. In particular, the user can change the "Minimum Frequency Used", the "Min FFT Peak Height %", and the value of "Delta".
What are these?

Minimum Frequency Used: The is the lowest Fourier spectrum frequency at which iFAMS will try to find a local maximum (a charge state fundamental peak). The default is 10 (corresponding to a m/z comb spacing of 1/10 of the mass spectrum width), which should work
in most cases. This helps to avoid analyzing the (typically very large) peak centered at 0 frequency, which tends to contain little helpful information. If you think (based on the blue dots that appear in the Fourier spectrum) that
frequencies that are too low are being analyzed, try adjusting the Minimum Frequency Used upward until you see the blue dots on those frequencies disappear. See NOTE, below.

Min FFT Peak Height %: This is the minimum % of the tallest peak in the Fourier spectrum window that will be identified as a peak (and get a blue dot). If your Fourier spectrum is very noisy, you may need to adjust this upward to
correctly identify peaks. See NOTE, below.

Delta: This is the number of Fourier spectrum data points to either side of a blue dot that are used to find the Fourier peak centroids, the sub-unit mass, and ultimately the charge states. If you think the blue dots in your
Fourier spectrum are much too close (i.e., IFAMS is finding too many local maxima), try adjusting Delta upward. See NOTE, below.

NOTE: If you adjust any of these three parameters, you must input all the others. For example, if all you want to do is adjust Delta from 5 to 10, you should type "10" below "Minimum Frequency Used" and "10" below "Min FFT Peak Height"
so that those values remain at their defaults. Once you have entered values in all three of the parameter boxes, click "Update These Values". The Fourier spectrum display should change and show you the results using the new parameters.
If you want to recalculate the charges and subunit mass and/or plot the charge-state specific mass spectra, you should click the corresponding buttons again. Please note that this will overwrite the old data in the charge-state
specific mass spectrum output files ("...IFFT...csv").



FAQs:

1.) Can iFAMS load files that are not .csv or .txt? Currently, no.
2.) What happens if a mass spectrum contains many m/z values for which the abundance is 0 (sparse data)? iFAMS interpolates the mass spectral data with a cubic spline. If any abundance values from this interpolation are negative,
they will be replaced with zeroes. The "Warnings" box will let you know for how many m/z values this has been done.
3.) What happens if a mass spectrum contains multiple different abundance values for the same m/z? iFAMS will spread these values out between nearest m/z neighbors of the values in question, and the "Warnings" box will let you know
for how many values this has been done. It is preferable to just start out with m/z values reported with as many decimal places as possible.
4.) What does the Warning "too few charge states" mean? If iFAMS thinks there are fewer than three adjacent charge states from the Fourier data, it will return this message to let you know you may need data with higher signal-to-noise.
It will still perform the Sub-unit Mass and Charge State calculations.
5.) iFAMS got the subunit mass and charge state completely wrong. What's going on? There are a few cases in which iFAMS may return information that you intuitively know cannot be correct (from knowledge of your sample):

5a.) The subunit mass iFAMS finds is very wrong. This most frequently happens when the mass spectral envelopes for each charge state are too narrow. As a rule of thumb, the mass spectrum for each charge state needs
to be at least about 2 times the subunit mass in width, i.e., if you have ions that are 10+ and the subunit mass is about 500 Da, you'll need the mass spectral envelope for the 10+ charge state to be about 1000 m/z wide (it should
have 20 or so peaks in it). If the mass spectral envelope is much narrower than this, the Fourier peaks will be broad and can interfere with each other, resulting in erroneous mass determination.

5b.) The subunit mass iFAMS finds is about half what it should be. This can occur if the "base mass" of each charge-state specific mass comb is very nearly an exact integer multiple of the subunit mass. For example, if you look
at protonated lipid clusters, the "base mass" is just the mass of however many protons are present in the clusters and will be very close to zero as compared to the mass of the lipid. What this does is make every other Fourier
spectrum peak very small (due to interference). Similarly, if you are looking at a protein that weighs almost exactly 10x the mass of the subunit you're interested in, you may see something similar. If you suspect this is what is going on, consider doubling the subunit mass and charge states you obtain using iFAMS. This
outcome is unfortunately a result of the fundamental ambiguity of the peak spacing in such a mass spectrum!

5c.) There are large interferences from chemical noise in your mass spectrum. For example, protonated lipid clusters and salt clusters can appear at low m/z in mass spectra of lipidated proteins. These should be trimmed off of the mass spectrum
before using iFAMS, if possible, by restricting the input data to a range where they are not present. Another option, of course, is to improve sample preparation so that these chemical interferents are eliminated as much as possible. Two-dimensional
Fourier-related methods can also be used, and automated versions of these will be released in future versions of iFAMS. 