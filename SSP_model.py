# -*- coding: utf-8 -*-
"""
Edited on Fri Mar 29 13:22:07 2019

@author: rodend
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:40:01 2019
@author: chris
"""
import numpy as np
import os
from datetime import datetime
from astropy.io import ascii
import matplotlib.pyplot as plt
import stpars
import itertools
import retrieve_irtf as ret
from scipy.interpolate import interp1d
from pyphot.phot import Filter
import pyphot

print('------ FUNCTION SSP.MODEL ------')
print('--------------------------------')
print('PURPOSE:')
print('   Create SSP spectra')
print('CALLING SEQUENCE:')
print('   ssp.model(feh, afe, age, imf, slope, fwhm, dl, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe,')
print('             n_ms, n_rg, parfile, lmin, lmax)')
print('INPUTS:')
print('   feh   = iron abundance [Fe/H]')
print('   afe   = [alpha/Fe]')
print('   age   = age of the isochrone (Gyr)')
print('   imf   = IMF (Kroupa, Salpeter or Unimodal)')
print('   slope = IMF slope (if IMF = Unimodal)')
print('   fwhm  = spectral resolution (default = 0.2 A)')
print('   dl    = sampling delta lambda (default = 0.1 A/pixel)')
print('   CFe   = [C/Fe] (set to 0.0 if omitted, i.e., default [C/Fe] = solar)')
print('CFe_rgb  = [C/Fe] of RGBs with logg <= logg_cn (set to CFe if omitted, i.e., default [C/Fe] = CFe)')
print('   NFe   = [N/Fe] (set to 0.0 if omitted, i.e., default [N/Fe] = solar)')
print('NFe_rgb  = [N/Fe] of RGBs with logg <= logg_cn (set to NFe if omitted, i.e., default [N/Fe] = NFe)')
print('   OFe   = [O/Fe] (set to [alpha/Fe] if omitted, i.e., default [O/Fe] = afe)')
print('   MgFe  = [Mg/Fe] (set to [alpha/Fe] if omitted, i.e., default [Mg/Fe] = afe)')
print('   SiFe  = [Si/Fe] (set to [alpha/Fe] if omitted, i.e., default [Si/Fe] = afe)')
print('   CaFe  = [Ca/Fe] (set to [alpha/Fe] if omitted, i.e., default [Ca/Fe] = afe)')
print('   TiFe  = [Ti/Fe] (set to a[alpha/Fe]fe if omitted, i.e., default [Ti/Fe] = afe)')
print('   NaFe  = [Na/Fe] (set to 0.0 if omitted, i.e., default [Na/Fe] = solar)')
print('   AlFe  = [Na/Fe] (set to 0.0 if omitted, i.e., default [Al/Fe] = solar)')
print('   BaFe  = [Na/Fe] (set to 0.0 if omitted, i.e., default [Ba/Fe] = solar)')
print('   EuFe  = [Na/Fe] (set to 0.0 if omitted, i.e., default [Eu/Fe] = solar)')
print('   n_ms  = number of main sequence stars')
print('   n_rg  = number of red giant stars')
print('logg_cn  = stars with logg <= logg_cn --> RGBs with different [C/Fe] and [N/Fe] ratios (CFe_rgb and NFe_rgb), default logg_cn = 3.0')
print(' parfile = file with list of stellar parameters (used only if n_ms/n_rg are not specified)')
print('   lmin  = lower lambda')
print('   lmax  = upper lambda')
print('OUTPUT:')
print('   SSP spectra in folder ./SSP_Spectra/' )
print('   logfile: SSP_model.log' )
print('REQUIRED FUNCTIONS:')
print('   nulbad (fortran code)')
print('   set.stpars.filename (R function, defined in stpars.R)')
print('   get.tracks          (R function, defined in stpars.R)')
print('   set.stspec.filename (R function, defined in pfant12.R)')
print('EXAMPLES:')
print("   ssp.model(0.2, 0.2, 8, 'Kroupa', n_ms = 9, n_rg = 6)")
print("   ssp.model(0.2, 0.2, 8, 'Unimodal', slope = 2.3, n_ms = 9, n_rg = 6)")

####################################################################
# FUNCTION SSP.MODEL
####################################################################
def ssp_model(Z, feh = None, afe = None, age = None, imf = None, slope = None,
             fwhm = 0.2, dl = 0.1, CFe = 0.0, CFe_rgb = None, NFe = 0.0, 
             NFe_rgb = None, OFe = None, MgFe = None, SiFe = None, CaFe = None,
             TiFe = None, NaFe = 0.0, AlFe = 0.0, BaFe = 0.0, EuFe = 0.0,
             n_ms = None, n_rg = None, logg_cn = 3, parfile = None, iso = 'DARTMOUTH'):
    
    #---------------------------------
    # CHECKING INPUTS
    #---------------------------------
    if feh is None:
        feh = input('Need to specify [Fe/H]: ')
    if afe is None:
        afe = input('Need to specify [alpha/Fe]: ')
    if age is None:
        age = input('Need to specify population age in Gyr: ')
    if imf is None:
        imf = input('Need to specify IMF (Salpeter, Kroupa, Unimodal): ')
    imf = imf.lower()
    
    if imf == 'unimodal' and slope is None:
        print('IMF = unimodal, you need to specify a slope')
    
    if n_ms is None and n_rg is None and parfile is None:
        print('FILE WITH STELLAR PARAMETERS --> ???')
        print(' Need to specify n_ms & n_rg OR parfile')
    elif n_ms is not None and n_rg is not None:
        parsfile = stpars.set_stpars_filename(n_ms, n_rg, feh, afe, age)
    else:
        parsfile = parfile
    
    #---------------------------------
    # SETTING DEFAULT VALUES
    #---------------------------------    
    if CFe_rgb is None:
        CFe_rgb = CFe
    if NFe_rgb is None:
        NFe_rgb = NFe
    if OFe is None:
        OFe = afe
    if MgFe is None:
        MgFe = afe
    if SiFe is None:
        SiFe = afe
    if CaFe is None:
        CaFe = afe
    if TiFe is None:
        TiFe = afe
        
    #---------------------------------
    # OUTPUT FILE NAMES
    #---------------------------------
    exists = os.path.isfile('SSP_Spectra/')
    if exists is False:
        os.system('mkdir SSP_Spectra')
        
    file_ssp = set_ssp_filename(feh, afe, age, imf, slope, CFe, CFe_rgb, NFe, NFe_rgb, OFe, MgFe,
                              SiFe, CaFe, TiFe, NaFe, AlFe, BaFe, EuFe)
    file_fig = file_ssp + '.jpg'
    logfile = file_ssp + '.log'
    
    #---------------------------------
    # WRITE INFO IN LOGFILE
    #---------------------------------
    log = open(logfile,'w+')
    
    log.write('----------------------- INPUT PARAMETERS -----------------------\n')
    log.write('----------------------------------------------------------------\n')
    
    log.write('System time:  ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S') 
                                                                    + '\n')
    
    temp = '%-26s%8.2f%4s' % ('Isochrone age: ', age, 'Gyr')
    log.write(temp + '\n')
    temp = '%-26s%8.2f%2s' % ('FWHM: ',fwhm, ' A')
    log.write(temp + '\n')

    if imf == 'unimodal':
        temp = '%-26s%8s%8s%8.2f' % ('IMF: ', imf, ' slope: ', slope)
    else:
        temp = '%-26s%8s' % ('IMF: ', imf)
    log.write(temp + '\n')
    
    log.write('Abundances:\n')
    temp = '%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s' % ('[Fe/H]', '[a/Fe]', 
                  '[C/Fe]', '[N/Fe]', '[O/Fe]', '[Mg/Fe]', '[Si/Fe]', 
                  '[Ca/Fe]', '[Ti/Fe]', '[Na/Fe]')
    log.write(temp + '\n')
    
    temp = '%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f' % (feh, afe,
                                CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe, NaFe)
    log.write(temp)
    
    
    #---------------------------------
    # CONSTANTS
    #---------------------------------
    Lsun = 3.839e33  # joules s-1 aka watts
    pc = 3.086e18 # parsec in cm
    G = 6.67300e-11  # m3 kg-1 s-2
    Msun = 1.9891e33 # g
    sb = 5.67e-5     # erg cm-2 s-1 K-4
    Rsun = np.sqrt(Lsun / (4 * np.pi * sb * 5777**4))*(10**(-2)) # = 6.955e10 cm
    aSun = 4 * np.pi * Rsun**2 # cm2
    
    #---------------------------------
    # READ INPUTS FILE
    #---------------------------------
    t = ascii.read(parsfile)

    nstars = len(t)
    
    Teffs = np.zeros(nstars)    
    loggs = np.zeros(nstars)
    masses = np.zeros(nstars)
    lumis = np.zeros(nstars)
    phase = list(itertools.repeat("ms", nstars))
    
    for i in range(nstars):
        Teffs[i] = t[i][0]
        loggs[i] = t[i][1]
        masses[i] = t[i][2]
        lumis[i] = 10**t[i][3]
        phase[i] = t[i][4]
    
    radii = np.sqrt(lumis * Lsun / (4 * np.pi * sb * Teffs**4)) # cm
    radii = radii*(10**(-2)) # m
    area = 4 * np.pi * radii**2 # m2
    
    
    log.write('\nReading stellar parameters from:  ' + parsfile + '\n')
    log.write('----------------------------------------------------------------\n\n')
    
    #---------------------------------
    # MASS BINS
    #---------------------------------
    logL = np.log10(lumis)
    logL_breaks = np.zeros(len(logL))
    nn = len(lumis)
    
    lum_breaks = np.zeros(nn+1)
    
    for i in range(nn+1):
        if i == 0:
            lum_breaks[i] = lumis[0]
        if i > 0 and i < nn:
            lum_breaks[i] = np.mean([lumis[i-1],lumis[i]])
        if i == nn:
            lum_breaks[i] = lumis[nn-1]

    
    isofile = stpars.gettracks(feh, afe, age, iso = iso)
    F0 = 1.021e-20 #Vega flux in erg cm-2 s-1 Hz-1
    t = np.loadtxt(isofile)
    if iso.upper() == 'DARTMOUTH':
	isoteff = 10**t[:,2]
        isologg = t[:,3]
        isomass = t[:,1]
        isologL = t[:,4]
    if iso.upper() == 'PADOVA':
        isoteff = 10**t[:,6]
        isologg = t[:,7]
        isomass = t[:,2]
        isologL = t[:,5]
	isoHFlux = 10**(-t[:,30]/2.5)*F0
        isologLH = np.log10(isoHFlux*4*np.pi*(10*pc)**2)
    print(iso)

    isoL = 10**isologLH
    mass_breaks = interp(isologLH, np.log10(lum_breaks), isomass)
    mass_breaks[0] = masses[0]
    mass_breaks[-1] = masses[-1]


    
    # Correct mass_breaks if interpolation fails for stars close to the turn-off
    for i in range(nn):
        if mass_breaks[i] > mass_breaks[i+1] or mass_breaks[i] > masses[i]:
            mass_breaks[i] = masses[i] - (masses[i] - masses[i-1])/2
            logL_breaks[i] = logL[i] - (logL[i] - logL[i-1])/2



    #---------------------------------
    # SSP SPECTRA
    #---------------------------------
    fraction_M = np.zeros(len(masses))
    fraction_L = np.zeros(len(masses))
    sum_flux = np.zeros(len(masses))
    L_corr = np.zeros(len(masses))
    #inter = interp1d(isomass, isoHmag, fill_value = "extrapolate")
    
    log.write('SSP divided in:  ' + str(nstars) + ' mass bins\n')
    log.write('Information about each bin: \n')
    log.write('----------------------\n')
    data = '%5s%10s%10s%10s%10s%7s%7s%10s%10s%10s%10s%10s%10s%10s' % ('i', 
        'logL_i', 'logL_f', 'Mass_i', 'Mass_f', 'Teff', 'logg', 'Mass',
        'logLstar', 'Radius', 'logLbin', 'Lcorr', 'frac_M', 'frac_L')
    log.write(data + '\n')
    log.write('    -----------------------------------------------------------------------------------------------------------------------------\n')
    
    for i in range(nstars):
        #if phase[i] != 'rgb_cn':
        #    file_flux = pfant12.set_stspec_filename(feh, afe, lmin, lmax, Teffs[i], 
        #                               loggs[i], CFe, NFe, OFe, MgFe, SiFe, 
        #                                CaFe, TiFe, NaFe, AlFe, BaFe, EuFe)
        #else:
        #    file_flux = pfant12.set_stspec_filename(feh, afe, lmin, lmax, Teffs[i], 
        #                               loggs[i], CFe, NFe, OFe, MgFe, SiFe,
        #                                CaFe, TiFe, NaFe, AlFe, BaFe, EuFe)
    
	file_flux = ret.set_spectra_name(Teffs[i], loggs[i], Z)
	print('-------------------------------------------------------------------------------')
	print('Implementing file: ' + file_flux)
        #---------------------------------
        # CONVOLUTION STELLAR SPECTRA
        #---------------------------------
	
	os.system('cp ./Stellar_Spectra/' + file_flux + ' ./')
	os.system('mv ./' + file_flux + ' st_spectra')
	
        #os.system('rm -f st_spectra')
        #print(file_flux)
        #os.system('cp ' + file_flux + ' temp.in')
        #nulbad_call = "nulbad --fn_flux temp.in --fn_cv st_spectra --flam T" \
        #    + " --pat %f --fn_progress progress.txt --fwhm %f" % (dl, fwhm)
        
        #os.system(nulbad_call)
        
        #---------------------------------
        # READ CONVOLVED STELLAR SPECTRA
        #---------------------------------
        tt = np.loadtxt('st_spectra')
        lamda = tt[:,0]
        hfilt = pyphot.get_library()['GROUND_BESSELL_H']
	units = 1
        #units = 1e-8 * 4 * np.pi  # erg / (s cm2 cm ster) --> erg / (s cm2 A)
        # Scale flux according to the star surface area and normalize by Lsun
        st_flux = tt[:,1]#*hfilt.get_flux(tt[:,0],tt[:,1])# * area[i] / Lsun # erg / (s cm2 A) --> Lsun / A
        
        if i == 0:
            ssp_flux = np.zeros(len(st_flux))
        
        # Number of stars formed between mass_breaks[i] and mass_breaks[i+1] (per unit of mass):
        #  sum(phi_m * dm) = phi_mdm
        dm_break = mass_breaks[i+1] - mass_breaks[i]
#        for k in range(len(mass_breaks)):
#            print(str(mass_breaks[k]))
        dm = dm_break/1000
        if dm > 0.01:
            dm = 0.001
        
#        ndm = round((mass_breaks[i+1] - mass_breaks[i]) / dm)
        masses_temp = np.arange(mass_breaks[i], mass_breaks[i+1], dm)
	
        iso_L_temp = 10**interp(isomass, masses_temp, isologLH)


#        print(str(iso_L_temp))
        phi_mdm = 0
        mass_bin = 0
        L_bin = 0
        L_star_bin = 0
        
        #---------------------------------
        # IMF: KROUPA
        #---------------------------------
        if imf == 'kroupa':
            for j in range(len(masses_temp)):
                cc = 1
                if masses_temp[j] >= 1.0 and masses_temp[j] < 100: x = 2.7
                if masses_temp[j] >= 0.5 and masses_temp[j] < 1:
                    x = 2.3
                    cc = ((1.0**(-2.7)) / 1.0**(-2.3))
                if masses_temp[j] >= 0.08 and masses_temp [j] < 0.5:
                    x = 1.3
                    cc = ((0.5**(-2.3)) / 0.5**(-1.3)) * \
                        ((1.0**(-2.7)) / 1.0**(-2.3))
                if masses_temp[j] >= 0.01 and masses_temp[j] < 0.08:
                    x = 0.3
                    cc = (0.08**(-1.3)) / 0.08**(-0.3) * ((0.5**(-2.3)) /
                            0.5**(-1.3)) * ((1.0**(-2.7)) / 1.0**(-2.3))
                #Hmag = inter(masses_temp[j])
                phi_m = masses_temp[j]**(-x) * cc
                phi_mdm = phi_mdm + phi_m * dm
                mass_bin = mass_bin + masses_temp[j] * phi_m * dm
                L_bin = L_bin + iso_L_temp[j] * phi_m * dm# * Hmag
                L_star_bin = L_star_bin + lumis[i] * phi_m * dm #* Hmag
        
        #---------------------------------
        # IMF: SALPETER
        #---------------------------------
        if imf == 'salpeter':
            x = 2.3
            for j in range(len(masses_temp)):
                phi_m = masses_temp[j]**(-x)
                phi_mdm += phi_m * dm
                mass_bin += masses_temp[j] * phi_m * dm
                L_bin += iso_L_temp[j] * phi_m * dm
                L_star_bin += lumis[i] * phi_m * dm
        
        #---------------------------------
        # IMF: UNIMODAL (SLOPE=x)
        #---------------------------------
        if imf == 'unimodal':
            x = slope
            for j in range(len(masses_temp)):
                phi_m = masses_temp[j]**(-x)
                phi_mdm += phi_m * dm
                mass_bin += masses_temp[j] * phi_m * dm
                L_bin += iso_L_temp[j] * phi_m * dm
                L_star_bin += lumis[i] * phi_m * dm
        
        #---------------------------------
        # SUM(S * PHI * DM)
        #---------------------------------
        L_corr[i] = L_bin / L_star_bin
        ssp_flux += st_flux * L_corr[i] * phi_mdm
        fraction_M[i] = mass_bin
        #sum_flux[i] = np.sum(st_flux) * (max(lamda) - min(lamda))
        fraction_L[i] = L_bin
        
    # normalization: sum(mass * phi_m * dm) = 1 Msun
    L_SSP = np.sum(fraction_L) 
    fraction_L = fraction_L / L_SSP
    ssp_flux = ssp_flux / np.sum(fraction_M)
    ssp_flux = ssp_flux  * 10**-7 * hfilt.get_flux(lamda,ssp_flux) # Hband unnormalisation and erg conversion
    fraction_M = fraction_M / np.sum(fraction_M)
    
    #---------------------------------
    # WRITE INFO IN LOGFILE
    #---------------------------------
    for i in range(nstars):
        data = '%5i%10.4f%10.4f%10.4f%10.4f%7.0f%7.2f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.2f' % (i+1, np.log10(lum_breaks[i]),\
                np.log10(lum_breaks[i+1]), mass_breaks[i], mass_breaks[i+1],\
                Teffs[i], loggs[i], masses[i], np.log10(lumis[i]),\
                radii[i] / Rsun, np.log10(L_SSP * fraction_L[i]), L_corr[i],\
                fraction_M[i]*100, fraction_L[i]*100)
        log.write(data + '\n')
    
    log.write('----------------------\n')
    
    data = ['logL_i, logL_f   -> lower and upper luminosity limits of the'
                                                    + ' bin (log[L/Lsun])', 
           'Mass_i, Mass_f   -> lower and upper mass limits of the bin (Msun)',
           'Teff, logg       -> stelar parameter of the synthetic stellar '
                           + 'spectrum representing the stars within bin',
           'Mass, logLstar   -> mass and luminosity of the star with'
                                                           + ' Teff, logg',
           'Radius           -> radius of a star with Teff, logg (Rsun,'
                       + ' Radius = sqrt[L / (4 * pi * sb * Teffs**4)])',
           'logLbin          -> total luminosity of the bin (log[L/Lsun])',
           'Lcorr            -> luminosity correction factor taking into '
                       + 'account the variation of stellar luminosities',
           '                      within bin ( = sum[L(m)*phi_m*dm]'
                                               + ' / sum[Lstar*phi_m*dm])',
           'frac_M           -> fraction of mass in bin i (%, frac_M = '
                                                       + '(Mbin/M_SSP)*100)',
           'frac_L           -> fraction of light in bin i (%, frac_L = '
                                                       + '(Lbin/L_SSP)*100)']
    for i in range(len(data)):
        log.write(data[i] + '\n')
    
    datum = open(file_ssp, 'w+')
    data = '%7.2f%13.5e' % (lamda[i], ssp_flux[i])
    datum.write(data + '\n')

    for i in range(1,len(lamda)):
        data = '%7.2f%13.5e' % (lamda[i], ssp_flux[i])
        datum.write(data + '\n')
    
    os.system('cp SSP_Spectra ' + file_ssp)
    
    log.write('\n    SSP spectra -->  ' + file_ssp)
    log.write('----------------------------------------------------------------\n')
    log.close()
    datum.close()

    print('-------------------------------------------------------------------------------')
    print('SSP succesfully created and saved in --> ' + file_ssp)
    
    #---------------------------------
    # PLOT SSP SPECTRA
    #---------------------------------
    
    t = np.loadtxt(file_ssp)
    t[:,1] = t[:,1]
    #t[:,1] = 3e-9*t[:,1]/(t[:,0]*10**-4)**2
#    t1 = np.zeros(len(t))
#    t2 = np.zeros(len(t))
#    for i in range(len(t)):
#        t1[i] = t[i][0]
#        t2[i] = t[i][1]
    plt.figure()
    
#    plt.yscale('log')
    plt.xlabel(r'$\lambda (\AA)$')
    plt.ylabel('Flux')
    #plt.savefig(file_fig)

    t_2 = np.loadtxt('./DATA/MarS/MARv_SAL_sed_NOCS_H_Z_0.029999999_Tg_1.0000000e+10')
    t_3 = np.loadtxt('./DATA/GirS/GIRv_SAL_sed_NOCS_H_Z_0.029999999_Tg_1.0000000e+10')
    t_4 = np.loadtxt('./DATA/BaSS/BASv_SAL_sed_NOCS_H_Z_0.029999999_Tg_1.0000000e+10')
    t_2[:,0] = t_2[:,0] *10000
    t_3[:,0] = t_3[:,0] *10000
    t_4[:,0] = t_4[:,0] *10000

    ax = plt.subplot(111)
    

    inter=interp1d(t_2[:,0],t_2[:,1])
    t_2norm = inter(12230)
    ax.plot(t_2[:,0],t_2[:,1],'r', linewidth = 0.25, label = 'MarS Model')
    inter=interp1d(t_3[:,0],t_3[:,1])
    ax.plot(t_3[:,0],t_3[:,1],'g', linewidth = 0.25, label = 'GirS Model')
    inter=interp1d(t_4[:,0],t_4[:,1])
    ax.plot(t_4[:,0],t_4[:,1],'m', linewidth = 0.25, label = 'BaSS Model')
    inter=interp1d(t[:,0],t[:,1])
    tnorm = inter(12230)
    ax.plot(t[:,0],t[:,1],'b', linewidth = 0.25, label = 'Our Model')
    ax.legend()
    print()

    plt.show(block=True)

    
##################################
# FUNCTION INTERP
##################################
def interp(xin, xout, yin):
    Ninterpol = len(xout)
    yout = np.zeros(Ninterpol)
#    xin = np.array(xin)    
    
    
    
    for k in range(Ninterpol):
#        t = xin[np.where(xin<xout[k])[0]]
#        tind = len(t)-1
        tind = 0
        for i in range(len(xin)):
            if xin[i]<xout[k]:
                tind += 1
        tind += -1
        
        if tind <= 0:
            tind = 1
        if tind >= len(xin) - 1:
            tind = len(xin) - 2
        t1 = xin[tind-1]
        t2 = xin[tind]
        t3 = xin[tind+1]
        tx = xout[k]
        
        A = (tx - t1) / (t3 - t1)
        E = (tx - t2) / (t3 - t1)
        B = (tx - t2) / (t3 - t2)
        D = (tx - t1) / (t3 - t2)
        C = (tx - t3) / (t2 - t1)
        
        yout[k] = yin[tind+1] * A * B - yin[tind] * D * C + yin[tind-1] * E * C
    return(yout)
    
##################################
# FUNCTION SET.SSP.FILENAME
##################################
def set_ssp_filename(feh, afe, age, imf, slope, CFe = 0.0, CFe_rgb = None, NFe = 0.0,
                   NFe_rgb = None, OFe = None, MgFe = None, SiFe = None,
                   CaFe = None, TiFe = None, NaFe = 0.0, AlFe = 0.0,
                   BaFe = 0.0, EuFe = 0.0):
    #---------------------------------
    # SETTING DEFAULT VALUES
    #---------------------------------    
    if CFe_rgb is None:
        CFe_rgb = CFe
    if NFe_rgb is None:
        NFe_rgb = NFe
    if OFe is None:
        OFe = afe
    if MgFe is None:
        MgFe = afe
    if SiFe is None:
        SiFe = afe
    if CaFe is None:
        CaFe = afe
    if TiFe is None:
        TiFe = afe
    imf = imf.lower()
    
    temp = np.array([feh, afe, CFe, NFe, OFe, MgFe, SiFe, CaFe, TiFe,
                     NaFe, AlFe, BaFe, EuFe])
    temp_s = list()
    for i in range(len(temp)):
        if temp[i] < 0:
            temp_s.append('-' + str(abs(temp[i])))
        else:
            temp_s.append('+' + str(abs(temp[i])))
            
    if age < 10:
        age_s = '%1s%3.1f' % ('0',age)
    else:
        age_s = '%4.1f' % age
    
    if CFe_rgb == CFe and NFe_rgb == NFe:
        CN_rgb = ''
    else:
        if CFe_rgb >= 0:
            tempC = '%1s%4.2f' % ('+', CFe_rgb)
        else:
            tempC = '%1s%4.2f' % ('-', CFe_rgb)
        if NFe_rgb >= 0:
            tempN = '%1s%4.2f' % ('+', NFe_rgb)
        else:
            tempN = '%1s%4.2f' % ('-', NFe_rgb)
        
        CN_rgb = 'CNrgb%s_%s' % (CFe_rgb, NFe_rgb)
    
    if imf == 'unimodal':
        slope_s = '%4.2f' % slope
        file_ssp = './SSP_Spectra/SSP_Fe' + temp_s[0] + '_a' + temp_s[1] \
        + '_C' + temp_s[2] + '_N' + temp_s[3] + '_O' + temp_s[4] + '_Mg' \
        + temp_s[5] + '_Si' + temp_s[6] + '_Ca' + temp_s[7] + '_Ti' + temp_s[8] \
        + '_Na' + temp_s[9] + '_Al' + temp_s[10] + '_Ba' + temp_s[11] + '_Eu' \
        + temp_s[12] + CN_rgb + '_age' + age_s + '_slope' + slope_s
    else:
        file_ssp = './SSP_Spectra/SSP_Fe' + temp_s[0] + '_a' + temp_s[1] + '_C' \
        + temp_s[2] + '_N' + temp_s[3] + '_O' + temp_s[4] + '_Mg' + temp_s[5] \
        + '_Si' + temp_s[6] + '_Ca' + temp_s[7] + '_Ti' + temp_s[8] + '_Na' \
        + temp_s[9] + '_Al' + temp_s[10] + '_Ba' + temp_s[11] + '_Eu' + \
        temp_s[12] + str(CN_rgb) + '_age' + age_s + '_' + imf
    
    return(file_ssp)

