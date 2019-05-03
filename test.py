# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 15:11:56 2019

@author: chrisw
"""

import stpars
import SSP_model
import interp
import numpy as np
from astropy.io import ascii

n_ms = 30
n_rg = 30
Z = np.log10(0.03/0.012)
feh = 0.281
afe = 0.0
age = 10
lmin = 9.353139996530000644e+03
lmax = 2.410741666080859795e+04
imf = 'kroupa'
dl = 0.1

stpars.stpars(n_ms, n_rg, feh, afe, age, fig = True)

parsfile = stpars.set_stpars_filename(n_ms, n_rg, feh, afe, age)

t = ascii.read(parsfile)

nstars = len(t)

for i in range(nstars):
	interp.interpolate(t[i][0],t[i][1],Z)


SSP_model.ssp_model(Z, feh = feh, afe = afe, age = age, imf = imf, fwhm = 2.5, n_ms = n_ms, n_rg = n_rg, dl = dl)



#file_ssp = SSP_model.set_ssp_filename(feh, afe, age = age, imf = imf)
#t = ascii.read(file_ssp)
