# python_ssp
WIP. Run test.py to get an SSP


test.py sets inital values for all the variables. It then runs the following programs:
  - stpars - This code calculates the Teff, logg and Z components of several stars along an isochrone. The only isochrone I have atm is Padova age 10
  - interpall - This creates interpolated spectra for all the stars contained in the stpars output file
  - SSP_model - This combines all the spectra together to create an SSP.


The contained codes and their functions:
stpars.py
  - stpars uses the following required inputs:
    - n_ms = number of main sequence stars
    - n_rg = number of red giant stars
    - feh = iron abundance [Fe/H]
    - afe = [alpha/Fe]
    - age = age of the population (in Gyrs)
  - stpars uses the folowing additional inputs:
    - 
      
