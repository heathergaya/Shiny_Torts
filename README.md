# Shiny Torts
Conventional line-transect and enhanced line-transect distance sampling shiny app 

All analysis is processed in a Bayesian framework. 

This app is written to be gopher tortoise survey specific - however, the general framework can be used for a variety of other line-transect surveys.

Longer directions may be found inside the .zip, but here are some quick and dirty directions:

Step 1: Download the ShinyTorts.zip file and unzip

Step 2: Install JAGS (if not already on your computer) from https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/ 

Step 3: Install R-Portable into the folder ShinyTorts/dist. There is no need to open the application once it is installed. R-Portable can be downloaded here: https://sourceforge.net/projects/rportable/ 

Step 4: Open the “BayesLTDS” app (Windows Batch File) to run the shiny app! It will take a little while to install packages and similar when it first opens. You can check on the progress of the app at any point by opening the textfile called “error” in the ShinyTorts/log folder. 

Step 5: Try selecting the various options. Example data can be found in the ShinyTorts/app/shiny folder. 

Step 6: Set the run time to the desired length of time (Note: "0" hours will still take about 10 minutes to run)

Step 7: Hit the “Load Model” button – a little message should pop up saying the model is loaded.

Step 8: Press “Run Model” – again, a message will appear. You can monitor the running of the chains in the “error” file in the ShinyTorts/log folder

Step 9: When it finishes running, a table of parameters will appear at the bottom of the page. You will have to scroll down to see it. 

Step 10: Press “Plot Size Curve” to see a bias-adjusted curve of the burrow diameter.

Step 11: Navigate to ShinyTorts/app/shiny and make sure a “results_Bias_Adjusted_LTDS.csv” appeared, as well as a PDF called “Rplots”.


Example files are provided as examples for vegetation, burrow and transect data. 
 
 The model requires certain headers be present in each type of data for successful processing.
 
 For burrows:
 - Diameter (in cm; for non-tortoise surveys this could represent group size, body weight, etc.)
 - Occ (1 = occupied, 0 = not occupied, NA = unknown; if not measured, make all = 1; for non-tortoise surveys, this could be male/female or any other binary data)
 - Dist	(distance to transect in m)
 - Found (1 = found, 0 = missed; for standard LTDS surveys, all should = 1)
 - Veg (if vegetation measurements were taken they should be included here; all values should range between 0 and 1)
 
 For transects: 
  - "Length" is the only required column (must be in meters)
  
 For veg data (when included):
  - Veg (with values ranging between 0 and 1)
  
  
