## -*- coding: utf-8 -*-
"""
Edited on Fri Mar 29 13:34:24 2019

@author: rodend
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:49:35 2019
@author: chrisw
"""

import numpy as np
import os
import itertools
from astropy.io import ascii
import matplotlib.pyplot as plt

print('------ FUNCTION STPARS ------')
print('-----------------------------')
print('PURPOSE:')
print('   Get stellar parameters from stellar evolutionary tracks')
print('CALLING SEQUENCE:')
print('   stpars(n_ms, n_rg, feh, afe, age, logg_cn, fig)')
print('INPUTS:')
print('   n_ms = number of desired main sequence stars')
print('   n_rg = number of desired red giant stars')
print('   feh  = iron abundance [Fe/H] (Available [Fe/H] range is -0.5 -> 0.5')
print('   afe  = [alpha/Fe] (available values -0.2, 0.0,'
                                                + ' +0.2, +0.4, +0.6, +0.8')
print('                    [alpha/Fe] = +0.4, +0.6, +0.8 -> available only'
                                                        + ' for [Fe/H] <=0)')
print('   age  = age of the population (in Gyr; available ages: see manual)')
print('logg_cn = stars with logg <= logg_cn --> RGBs with different [C/Fe]'
                                + ' and [N/Fe] ratios, default logg_cn = 3.0')
print('   fig  = logical; set to T to plot the isochrone')
print('OUTPUT:')
print('   text file containing the stellar parameters (Teff, logg, '
                                                            + 'Mass, logL)' )
print('     of (n_ms + n_rg + 1) stars' )
print('REQUIRED SCRIPTS:')
print('   ./atm_models/iso_interp_feh (fortran code)')
print('   ./atm_models/isolf_split (fortran code)')
print('EXAMPLE:')
print('   stpars(9, 6, 0.2, 0.2, 8)')

####################################################################
# FUNCTION STPARS
####################################################################

def stpars(n_ms, n_rg, feh, afe, age, logg_cn = 3, fig = False, iso = 'DARTMOUTH'):
    
    #---------------------------------
    # SET OUTPUT FILENAME
    #---------------------------------
    fileout = set_stpars_filename(n_ms, n_rg, feh, afe, age, logg_cn)
    
    #---------------------------------
    # ISOCHRONE [FE/H]-INTERPOLATION
    #---------------------------------
    isofile = gettracks(feh, afe, age, iso = iso)

    if iso.upper() == 'DARTMOUTH':
        print('Utilizing the DARTMOUTH isochrones')
    if iso.upper() == 'PADOVA':
        print('Utilizing the Padova isochrones')
    
    #---------------------------------
    # READ ISOCHRONE
    #---------------------------------
    # load in the isochrone file set by gettracks
    t = np.loadtxt(isofile)
    F0 = 1.021e-20 #Vega flux in erg cm-2 s-1 Hz-1
    pc = 3.086e18 # parsec in cm
    # checks what isochrone the user wants to use and retrieves the corrosponding coloumns from each
    if iso.upper() == 'DARTMOUTH':
        #isoteff in the file is in log form
        isoteff = 10**t[:,2]
        isologg = t[:,3]
        isomass = t[:,1]
        isologL = t[:,4]
    if iso.upper() == 'PADOVA':
        #isoteff in the file is in log form
        isoteff = 10**t[:,6]
        isologg = t[:,7]
        isomass = t[:,2]
        isologL = t[:,5]
        # have to calculate the H band luminosity by taking the H band absolute magnitudes and turning them to luminosities
        isoHFlux = 10**(-t[:,30]/2.5)*F0
        isologLH = np.log10(isoHFlux*4*np.pi*(10*pc)**2)
    
    #---------------------------------
    # GET STELLAR PARAMETERS
    #---------------------------------
    # logg_lim is the logg divider between ms stars and rg stars.
    # logg_lim is equal to the logg value at the maximum teff value on the turning point of the isochrone
    logg_lim = max(isologg[isoteff == max(isoteff)])

    # min teff for ms and rg stars depending on whether they are greater or smaller than this logg_lim
    min_teff_ms = min(isoteff[isologg>=logg_lim])
    min_teff_rg = min(isoteff[isologg<logg_lim])

    # delta_teff's for both ms and rg stars. This is the difference between the max teff and the min teff for either ms or rg, then divided by the number of stars of each type selected    
    delta_teff_ms = (max(isoteff) - min_teff_ms) / n_ms
    delta_teff_rg = (max(isoteff) - min_teff_rg) / n_rg
    
    # teff_grid is the final grid we will put our new steller teffs in
    teff_grid = np.linspace(1, n_ms + n_rg + 1, n_ms + n_rg + 1)
    phase = list(itertools.repeat("ms", n_ms + n_rg + 1))
    
    # a loop setting all the ms star teff values to the minimum teff + some delta teff. This is linear scale.
    for i in range(n_ms):
        teff_grid[i] = min_teff_ms + i * delta_teff_ms

    #this then sits the teff_grid at the changing point from ms to rg to the max isochrone teff
    teff_grid[n_ms] = max(isoteff)
    # using this additional variable j is important here as the for loop goes backwards from the end of teff_grid
    j = 0
    # for loop sets all the values linearly from the end of teff grid to the max value at teff_grid[n_ms]
    for i in range(n_ms+n_rg, n_ms, -1):
        teff_grid[i] = min_teff_rg + j*delta_teff_rg
        phase[i] = "rgb"
        j += 1
    # calculates the index where the teff_grid is the max value of the isoteff
    index_lim = np.where(teff_grid == max(isoteff))[0][0]
    # grid of logg, mass and luminosity values
    logg_grid = np.zeros(len(teff_grid))
    mass_grid = np.zeros(len(teff_grid))
    lumi_grid = np.zeros(len(teff_grid))
    # for each teff_grid value
    for i in range(len(teff_grid)):
        # if the index is less than or greater than index_lim and so is a rg star or a ms star respectively, to set xxxx to the appropriate condition
        if i <= index_lim:
            xxxx = isologg >= logg_lim
        else:
            xxxx = isologg < logg_lim


        # This next part of the code calculates the closest point on the isochrone to each of the current teff values
            
        # calculates the differences between the current teff_grid value and all the other ms or rg teff values on the isochrone
        temp = abs(isoteff[xxxx] - teff_grid[i])
        # temp teff_grid for the selected star type
        teff_grid_temp = isoteff[xxxx]
        # condition is the index of the value where temp == min(temp)
        condition = np.where(temp == min(temp))[0]
        #new teff_grid[i] is equal to the value where all of the prior conditions are correct
        teff_grid[i] = teff_grid_temp[condition][0]
        # This is then repeated for each of logg, mass and luminosity
        logg_grid_temp = isologg[xxxx]
        logg_grid[i] = logg_grid_temp[condition][0]
        mass_grid_temp = isomass[xxxx]
        mass_grid[i] = mass_grid_temp[condition][0]
        lumi_grid_temp = isologLH[xxxx]
        lumi_grid[i] = lumi_grid_temp[condition][0]

        # sets phase for rgb_cn stars to rgb_cn
        if logg_grid[i] <= logg_cn:
            phase[i] = "rgb_cn"
    
    # renaming of variables
    isoteffgrid = teff_grid
    isologggrid = logg_grid
    isomassgrid = mass_grid
    isologLgrid = lumi_grid
    
    #---------------------------------
    # WRITE OUTPUT FILE
    #---------------------------------
    exists = os.path.isfile('./Stellar_pars')

    if exists is False:
        os.system("mkdir Stellar_pars")
    
    ascii.write([isoteffgrid,isologggrid,isomassgrid,isologLgrid,phase],
                fileout, names = ['#Teff/k','logg',
                'Mass/Msun','logL/Lsun','phase'], overwrite = True)
    
    #---------------------------------
    # PLOT ISOCHRONE
    #---------------------------------
    if fig is True:
        isoplot = plt.plot(isoteff, isologg)
        isogridplot = plt.plot(isoteffgrid, isologggrid, 'o')
        plt.xlabel('Temperature (K)')
        plt.ylabel('log g (dex)')
        plt.xlim([2000, 8000])
        plt.ylim([-1, 6])
        plt.gca().invert_xaxis()
        plt.gca().invert_yaxis()
        plt.legend((('Isochrone [Fe/H] = ' + str(feh) + ', [a/Fe] = ' + str(afe) + ', Age = ' + str(age) + ' Gyr'),('Selected Stellar Parameters')),loc = 'upper left',fontsize = 'small')

        plt.show()
    
    #---------------------------------
    # PRINT SOME INFORMATION
    #---------------------------------
    
    print('Isochrone: [Fe/H] = ' + str(feh))
    print('           [a/Fe] = ' + str(afe))
    print('           Age = ' + str(age) + 'Gyr')
    print('STELLAR PARAMETERS OF: ' + str(n_ms + n_rg + 1) + 'STARS IN THE')
    print('OUTPUT FILE: ' + fileout)
    
##################################
# FUNCTION GET.TRACKS
##################################
def gettracks(feh, afe, age, iso = 'DARTMOUTH'):
    # gettracks prepares the DARTMOUTH isochrones and then finds the correct isochrones, either Padova or DARTMOUTH
    # The output is the file path to the isochrone files
    afes = np.array([[-0.2],[0.0],[0.2],[0.4],[0.6],[0.8]])
    afesr = abs(afes-afe)
    temp = np.where(afesr==min(afesr))[0] + 1
    
    scriptf = open("iso.sh","w+")
    scriptf.write("./iso_interp_feh << DONE\n")
    scriptf.write("13\n")
    scriptf.write("1\n")
    scriptf.write(str(temp[0]) + "\n")
    scriptf.write(str(feh) + "\n")
    scriptf.write("isotemp\n")
    scriptf.write("DONE\n")
    
    scriptf.write("#\n")
    
    scriptf.write("./isolf_split << DONE\n")
    scriptf.write("isotemp\n")
    scriptf.write("DONE")
    scriptf.close()
    
    os.system("chmod 777 iso.sh")
    os.system("mv iso.sh ./DARTMOUTH/iso.sh")
    os.chdir("./DARTMOUTH/")
    os.system("./iso.sh > temp.txt")
    os.chdir("../")
    if iso.upper() == 'DARTMOUTH':
        if age < 10:
            isofileout = "./DARTMOUTH/a0" + str(int(age)*1000) + "isotemp"
        else:
            isofileout = "./DARTMOUTH/a" + str(int(age)*1000) + "isotemp"
    
    if iso.upper() == 'PADOVA':
        if age < 10:
            isofileout = "./DARTMOUTH/isochrones/Padova/iso0" + str(int(age)) + ".dat"
        else:
            isofileout = "./DARTMOUTH/isochrones/Padova/iso" + str(int(age)) + ".dat"
        
    return(isofileout)

##################################
# FUNCTION SET.STPARS.FILENAME
##################################
def set_stpars_filename(n_ms, n_rg, feh, afe, age, logg_cn = 3):
    # sets up the stellar paramater file name and outputs it.

    if feh < 0:
        feh_s = "%1s%4.2f" % ('-', abs(feh))
    else:
        feh_s = "%1s%4.2f" % ('+', abs(feh))
        
    if afe < 0:
        afe_s = "%1s%4.2f" % ('-', abs(afe))
    else:
        afe_s = "%1s%4.2f" % ('+', abs(afe))
        
    if age < 10:
        age_s = "%1s%3.1f" % ('0', abs(age))
    else:
        age_s = "%4.1f" % abs(age)
        
    if n_ms < 10:
        n_ms_s = "%1s%1i" % ('0', abs(n_ms))
    else:
        n_ms_s = "%2i" % abs(n_ms)
        
    if n_rg < 10:
        n_rg_s = "%1s%1i" % ('0', abs(n_rg))
    else:
        n_rg_s = "%2i" % abs(n_rg)
    logg_cn_s = "%1.1f" % logg_cn
    fileout = "./Stellar_pars/stpars_fe" + feh_s + "a" + afe_s + "age" + age_s \
        + "ms" + n_ms_s + "rg" + n_rg_s + "loggCN" + logg_cn_s + ".dat"
        
    return(fileout)
