# This is the python scripts used for generating landscape formation

In this repo, it includes:

1. main script
      Watershed.py: generating 2 watershed, and the same time calling functions for calculating and generating files and figures for slopes, drainage area, stream density, curvation, TMR, soil indicator files for Parflow model.
      
2. function: 
      CalandConv.py: called from Watershed.py, calculating the slopes, drainage area, stream density, curvation, TMR
      Soil_parameterInParfow_tclwritting.py: called from Watershed.py, calculating the soil thickness for each soil layer and generating soil indicator file as Parflow input
      pfbwriting.py: called from Soil_parameterInParfow_tclwritting.py, used as a functing to wring pfb file
