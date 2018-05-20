import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import linalg
import time
import random
import sys
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import signal
from scipy import io
import h5py
from scipy import stats

from electrodiffusion import Model
from electrodiffusion import makePSD # (Phi_of_t, N_t, delta_t): return(f, psd_max, location)
from electrodiffusion import PSD_of_LFP #() : # make the PSD of the LFP from Gratiy


file_name =  'mouse_1_lfp_trial_avg_3sec.h5'
h5=h5py.File(file_name,'r')
#print( h5.keys() )
lfp_on_flash = h5['lfp_on_flash'][...]
lfp_off_flash = h5['lfp_off_flash'][...]
z = h5['zdepth'][...]
time = h5['time'][...]

plt.plot(lfp_on_flash[10])
plt.show()

Dietzel1982_1 = Model(np.load('Dietzel1982_Phi_of_t.npy'), np.load('Dietzel1982_parameters.npy'), 'Dietzel1982', 'xkcd:teal') 
Cordingley1978    = Model(np.load('Cordingley1978_Phi_of_t.npy'), np.load('Cordingley1978_parameters.npy'), 'Cordingley1978', 'xkcd:indigo')
Halnes2016    = Model(np.load('Halnes2016_Phi_of_t.npy'), np.load('Halnes2016_parameters.npy'), 'Halnes2016', 'xkcd:olive')
Nicholson1987 = Model(np.load('Nicholson1987_Phi_of_t.npy'), np.load('Nicholson1987_parameters.npy'), 'Nicholson1987', 'xkcd:grey')
Herreras1993  = Model(np.load('Herreras1993_Phi_of_t.npy'), np.load('Herreras1993_parameters.npy'), 'Herreras1993', 'xkcd:teal')
#	EkstremeGradient = Model(np.load('EkstremeGradient_Phi_of_t.npy'), np.load('EkstremeGradient_parameters.npy'), 'EkstremeGradient','b')

Models = [ Halnes2016,  Dietzel1982_1,Nicholson1987, Cordingley1978]
sdModels =[Herreras1993]


for M in Models:
	M.f, M.psd, M.location = makePSD(M.phi, M.N_t, M.delta_t)
	plt.plot(np.log10(M.f[31:-1]), np.log10(M.psd[31:-1]), M.color,label = M.name)
	print(M.name, np.log10(M.f[100]), np.log10(M.psd[100]))
#	plt.plot(np.log10(f[31:-1]),np.log10(psd[31:-1]), M.color, label = M.name)

PSD_of_LFP()

plt.title('PSDs of a local field potential and diffusion potentials' )
plt.xlabel('$log_{10}(frequency)$ in Hz')
plt.ylabel('$log_{10}(PSD)$ in (mV)²/Hz')
#plt.plot(np.log10(f_Dietzel1982[31:-1]),np.log10(psd_Dietzel1982[31:-1]), 'xkcd:plum', label = 'diffusion')
#plt.plot(np.log10(f_Herreras1993[31:-1]), np.log10(psd_Herreras1993[31:-1]), 'xkcd:teal', label = 'SD diffusion')
plt.legend()

plt.savefig('PSD_Gratiy',dpi=500)
plt.show()

sys.exit()

