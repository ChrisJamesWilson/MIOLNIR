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

def stpars(n_ms, n_rg, feh, afe, age, logg_cn = 3, fig = False):
    
    #---------------------------------
    # SET OUTPUT FILENAME
    #---------------------------------
    fileout = set_stpars_filename(n_ms, n_rg, feh, afe, age, logg_cn)
    
    #---------------------------------
    # ISOCHRONE [FE/H]-INTERPOLATION
    #---------------------------------
    isofile = gettracks(feh, afe, age)
    
    #---------------------------------
    # READ ISOCHRONE
    #---------------------------------
    t = np.loadtxt(isofile)
    isoteff = 10**t[:,2]
    isologg = t[:,3]
    isomass = t[:,1]
    isologL = t[:,4]
    
    #---------------------------------
    # GET STELLAR PARAMETERS
    #---------------------------------
    logg_lim = max(isologg[np.where(isoteff == max(isoteff))])
    min_teff_ms = min(isoteff[np.where(isologg>=logg_lim)])
    min_teff_rg = min(isoteff[np.where(isologg<logg_lim)])
    
    delta_teff_ms = (max(isoteff) - min_teff_ms) / n_ms
    delta_teff_rg = (max(isoteff) - min_teff_rg) / n_rg
    
    teff_grid = np.linspace(1, n_ms + n_rg + 1, n_ms + n_rg + 1)
    phase = list(itertools.repeat("ms", n_ms + n_rg + 1))
    
    for i in range(n_ms):
        teff_grid[i] = min_teff_ms + i* delta_teff_ms

    teff_grid[n_ms] = max(isoteff)
    j = 0
    for i in range(n_ms+n_rg, n_ms, -1):
        teff_grid[i] = min_teff_rg + j*delta_teff_rg
        phase[i] = "rgb"
        j += 1
        
    index_lim = np.where(teff_grid == max(isoteff))[0]
    logg_grid = np.zeros(len(teff_grid))
    mass_grid = np.zeros(len(teff_grid))
    lumi_grid = np.zeros(len(teff_grid))
    
    for i in range(len(teff_grid)):
        
        if i <= index_lim:
            xxxx = np.where(isologg >= logg_lim)[0]
        else:
            xxxx = np.where(isologg < logg_lim)[0]
            
            
        temp = abs(isoteff[xxxx] - teff_grid[i])
        teff_grid_temp = isoteff[xxxx]
        condition = np.where(temp == min(temp))[0]
        teff_grid[i] = teff_grid_temp[condition][0]
        logg_grid_temp = isologg[xxxx]
        logg_grid[i] = logg_grid_temp[condition][0]
        mass_grid_temp = isomass[xxxx]
        mass_grid[i] = mass_grid_temp[condition][0]
        lumi_grid_temp = isologL[xxxx]
        lumi_grid[i] = lumi_grid_temp[condition][0]
#        
#        logg_grid[i] = isologg[(xxxx) and np.where(temp == min(temp))][0]
#        mass_grid[i] = isomass[(xxxx) and np.where(temp == min(temp))][0]
#        lumi_grid[i] = isologL[(xxxx) and np.where(temp == min(temp))][0]

        
        if logg_grid[i] <= logg_cn:
            phase[i] = "rgb_cn"
    
    isoteffgrid = np.round(teff_grid)
    isologggrid = np.round(logg_grid,2)
    isomassgrid = np.round(mass_grid,7)
    isologLgrid = np.round(lumi_grid,4)
    
    #---------------------------------
    # WRITE OUTPUT FILE
    #---------------------------------
    exists = os.path.isfile('./Stellar_pars')

    if exists == False:
        os.system("mkdir Stellar_pars")
    
    ascii.write([isoteffgrid,isologggrid,isomassgrid,isologLgrid,phase],
                fileout, names = ['#Teff/k','logg',
                'Mass/Msun','logL/Lsun','phase'], overwrite = True)
    
    #---------------------------------
    # PLOT ISOCHRONE
    #---------------------------------
#    plt.figure()
#    plt.plot()
#    lines(iso.teff, iso.logg, lwd = 2, col = 'red')
#    points(iso.teff.grid, iso.logg.grid, col = 'red', pch = 19, cex = 1.4)
    isoplot = plt.plot(isoteff, isologg)
    isogridplot = plt.plot(isoteffgrid, isologggrid, 'o')
    plt.xlabel('Temperature (K)')
    plt.ylabel('log g (dex)')
    plt.xlim([2000, 8000])
    plt.ylim([-1, 6])
    plt.gca().invert_xaxis()
    plt.gca().invert_yaxis()
    plt.legend((('Isochrone [Fe/H] = ' + str(feh) + ', [a/Fe] = ' + str(afe) + ', Age = ' + str(age) + ' Gyr'),('Selected Stellar Parameters')),loc = 'upper left',fontsize = 'small')
#    plt.text(4750,0.5,('Isochrone [Fe/H] = ' + str(feh) + '\n[a/Fe] = '+ str(afe) \
#    + '\nAge = ' + str(age) + ' Gyr'))
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
def gettracks(feh, afe, age):
    afes = np.array([[-0.2],[0.0],[0.2],[0.4],[0.6],[0.8]])
    afesr = abs(afes-afe)
    temp = np.where(afesr==min(afesr))[0] + 1
    
    scriptf = open("iso.sh","w+")
    scriptf.write("./iso_interp_feh << DONE\n")
    scriptf.write("11\n")
    scriptf.write("1\n")
    scriptf.write(str(temp[0]) + "\n")
#    scriptf.write('2\n')
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
    
    if age < 10:
        isofileout = "./DARTMOUTH/a0" + str(int(age*1000)) + "isotemp"
    else:
        isofileout = "./DARTMOUTH/a" + str(int(age*1000)) + "isotemp"
        
    return(isofileout)

##################################
# FUNCTION SET.STPARS.FILENAME
##################################
def set_stpars_filename(n_ms, n_rg, feh, afe, age, logg_cn = 3):

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
