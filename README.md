# oscillations_paper


Code to reproduce the figures and analyses presented in: A model for the decrease in amplitude of carbon isotope excursions throughout the Phanerozoic, Aviv Bachan, Kimberly V. Lau, Matthew R. Saltzman, Ellen Thomas, Lee R. Kump, and Jonathan L. Payne. American Journal of Science, Vol. 317, June, 2017, P. 641-676, DOI 10.2475/06.2017.01

Files are as follows:

PhanData.mat & Hettangian_LBC_d13C.mat: contain the Phanerozoic data compilation from Saltzman and Thomas, 2012, and Upper Triassic -- Lower Jurassic data from Bachan et al. 2012, respectively. 

getData.m: Loads the data, cleans it up a bit, and plots it if you wish along with a smoothing spline fitted to the data. Is the main data loading routine which is called by the other files. 

Figure1Plot.m: Makes a plot of the Phanerozoic data with close-up views of various intervals of interest. This is Figure 1 in the paper. 

delCFourier.m: Carries out a variety of spectral analyses on the data which make up Figures 2, 3, and 4 in the paper, and Figure A1 in the supplement. 

BahamasModeFig.ipynb: Plots the data from Swart et al and carries out the regression to back out the fraction of shelf versus pelagic carbonate. This is Figure 6 in the paper. 


Oscillations.m: This file runs the model with internal oscillations presented in Figure 8. 

OscillationsModel.m, OscillationsSense.m, & PlotOscillationsOutput.m produce Figures 9, 10 & 11, as well as the supplemental figures A2 A, B, C, and A3. OscillationsModel contains the model itself. OscillationsSense repeatedly runs the model with forcings of various amplitudes and frequencies. PlotOscillationsOutput plots the output. Note that the runs can take a long time (a few hours). Pre-computed output can be found at: https://drive.google.com/open?id=0BzQ--Chu8wozTVliSnlQZkF1SWM
