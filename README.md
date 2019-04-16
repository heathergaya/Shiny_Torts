# Enhanced_LTDS
Enhanced line-transect distance sampling shiny app 

Three files and two programs are mandatory for the shiny app to run:
  - Shiny_BayesLTDS.R
  - 'Bias_Adjusted_Model_412019.R' 
  - RJMCMC.R 
  - JAGS (https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/)
  - R or Rstudio (https://www.rstudio.com/products/rstudio/download/) 
  
 Additional files are provided as examples for vegetation, burrow and transect data. 
 
 The model requires certain headers be present in each type of data for successful processing.
 
 For burrows:
 - Diameter (in cm)
 - Occ (1 = occupied, 0 = not occupied, NA = unknown; if not measured, make all = 1)
 - Dist	(distance to transect in m)
 - Found (1 = found, 0 = missed; for standard LTDS surveys, all should = 1)
 - Veg (if vegetation measurements were taken; all values should range between 0 and 1)
 
 For transects: 
  - "Length" is the only required column; must be in meters 
  
 For veg data (when included):
  - Veg (with values ranging between 0 and 1)
  
  
