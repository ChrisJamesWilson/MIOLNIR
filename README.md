# python_ssp
WIP. Run test.py to get an SSP


test.py sets inital values for all the variables. It then runs the following programs:
  - stpars - This code calculates the Teff, logg and Z components of several stars along an isochrone. The only isochrone I have atm is Padova age 10
  - interpall - This creates interpolated spectra for all the stars contained in the stpars output file
  - SSP_model - This combines all the spectra together to create an SSP.


The contained codes and their functions:

#stpars.py

stpars.py contains 3 different functions: stpars(), gettracks() and set_stpars_filename().

**stpars(n_ms, n_rg, feh, afe, age, logg_cn = 3, fig = False, iso = 'DARTMOUTH')**
  - stpars uses the following required inputs:
    - n_ms = number of main sequence stars
    - n_rg = number of red giant stars
    - feh = iron abundance [Fe/H]
    - afe = [alpha/Fe]
    - age = age of the population (in Gyrs)
    
  - stpars uses the following additional inputs:
    - logg_cn = stars with logg <= logg_cn
    - fig = logical, whether you want the isochrone plotted or not
    - iso = isochrone model used, currently supports DARTMOUTH or PADOVA

stpars finds the appropriate isochrone for the desired inputs and calculates the Teff and logg of the n_ms stars on the main sequence and n_rg on the red giant branch. It outputs a table of the Teff's, logg's, masses, Hband logL's and their respective phases in a .dat file found in the folder Stellar_pars.

**gettracks(feh, afe, age, iso = 'DARTMOUTH')**

gettracks uses the previously defined parameters to find the correct isochrone. It outputs the file path of the isochrone selected.

**set_stpars_filename(n_ms, n_rg, feh, afe, age, logg_cn = 3)**

This outputs the filename and path for the parameter file given the previously defined inputs.


#interp.py

interp.py contains two functions: interpall() and interpolate()

**interpolate(Teff_new, logg_new, Z_new)**
  - interpolate requires the following required inputs
      - Teff_new = the effective temperature of the new star n Kelvins
      - logg_new = the log(g) of the new star
      - Z = the metallicity of the new star in Z
  
This is the main interpolator. Given the above parameters it will interpolate a spectra from the data in the irtf folder. The spectra will be stored in Stellar_Spectra/ along with a plot of the spectra.
  
**interpall(n_ms, n_rg, feh, afe, age, Z)**

interpall's inputs have already been defined.

interpall uses interpolate() to calculate spectra for all stars that are in the paramater file created in **stpars.py**.
