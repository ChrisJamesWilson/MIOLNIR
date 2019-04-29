import numpy as np
from astropy.io import ascii
from astropy.io import fits
import itertools
import os

# Code pertaining to reading and also determining file names in the irtf library

# Enter location of the irtf fits files relative to this directory here:


def param_retrieve():
	## This retrieves all the paramter data and outputs it as an array
	irtf_file = 'irtf_param.txt'
	# first it reads the param file
	t = ascii.read(irtf_file)


	# number of star parameters
	numstars = len(t)
	
	# assigns the arrays for all the coloumns
	ID = list(itertools.repeat('.', numstars))
	Teff = np.zeros(numstars)
	logg = np.zeros(numstars)
	Z = np.zeros(numstars)
	

	# writes the values to the array
	for i in range(numstars):
		ID[i] = t[i][0]
		Teff[i] = t[i][1]
		logg[i] = t[i][2]
		Z[i] = t[i][3]
	
	# creates an array of all the params and their ids
	t = [ID,Teff,logg,Z]
	# returns an array with their ids
	return(t)

def get_spectra(ID):
	## retrieves a spectra form IRTF with ID number, IRL***. Returns an array 
	## of the spectra on the second column and the wavelength in the first
	
	# First it retrieves the data from the corresponding fits file
	t = fits.getdata('irtf/' + ID + '.fits')
	l0 = fits.getheader('irtf/' + ID + '.fits')['CRVAL1']
	dl = fits.getheader('irtf/' + ID + '.fits')['CDELT1']

	# creates a new array to put all the values in
	tt = np.zeros([len(t),2])
	tt[:,1] = t
	tt[:,0] = np.arange(l0,l0 + dl*(len(t)-dl),dl)
	
	# returns the new array
	return(tt)

def set_spectra_name(Teff, logg, Z):

	exists = os.path.isdir('Stellar_Spectra/')

	if exists is False:
		os.system('mkdir Stellar_Spectra')

	Teff_s = '%4.0f' % Teff

	if logg < 0:
		logg_s = '%1s%4.2f' % ('-', abs(logg))
	else:
		logg_s = '%1s%4.2f' % ('+', abs(logg))
	if Z < 0:
		Z_s = '%1s%4.2f' % ('-', abs(Z))
	else:
		Z_s = '%1s%4.2f' % ('+', abs(Z))	

	temp_s = [Teff_s,logg_s,Z_s]
	
	sp_file = 'intspectra_Teff' + temp_s[0] + '_logg' + temp_s[1]\
		+ '_Z' +temp_s[2]

	return sp_file



	
    
