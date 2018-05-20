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

def makePSD(Phi_of_t, N_t, delta_t):
	fs = 1/delta_t # sampling frequency
	psd_max = np.zeros(int(N_t/2 +1))

	for i in range(int(Phi_of_t.shape[0])-1):
	    f, psd_new = signal.periodogram(Phi_of_t[i,:], fs)
	    if np.amax(psd_new[1:-1]) > np.amax(psd_max[1:-1]):
	   		psd_max = psd_new
	   		location = i

	return(f, psd_max, location)

file_name =  'mouse_1_lfp_trial_avg_3sec.h5'
h5=h5py.File(file_name,'r')
#print( h5.keys() )
lfp_on_flash = h5['lfp_on_flash'][...]
lfp_off_flash = h5['lfp_off_flash'][...]
z = h5['zdepth'][...]
time = h5['time'][...]

plt.plot(lfp_on_flash[10])
plt.show()

Dietzel1982_1_Phi_of_t = np.load('Dietzel1982_1_Phi_of_t.npy')
[N_t, delta_t, N_x, delta_x] = np.load('Dietzel1982_1_parameters.npy')
f_Dietzel1982, psd_Dietzel1982, location = makePSD(Dietzel1982_1_Phi_of_t, N_t, delta_t)

Herreras1993_Phi_of_t = np.load('Herreras1993_Phi_of_t.npy')
f_Herreras1993, psd_Herreras1993, loc_Herreras = makePSD(Herreras1993_Phi_of_t, N_t, delta_t)



fs=2500 # sampling rate
x=time
N=len(x)
fmax=fs//2

psd = np.zeros((len(z),int(N/2+1)))

f_log = np.logspace(-1,2,100)
psd_on = np.zeros((len(z),N//2+1))
psd_off = np.zeros((len(z),N//2+1))

for ich,data_ch in enumerate(lfp_on_flash):
    f, Pxx_den = signal.periodogram(data_ch, fs)
    psd_on[ich,:]=Pxx_den
#    psd_on_log[ich,:] =np.interp(f_log, f, Pxx_den)
    
psd_off_log = np.zeros((len(z),len(f_log)))
    
for ich,data_ch in enumerate(lfp_off_flash):
    f, Pxx_den = signal.periodogram(data_ch, fs)
    psd_off[ich,:]=Pxx_den
#    psd_off_log[ich,:] =np.interp(f_log, f, Pxx_den)


psd_off_mean = np.mean(psd_off,axis=0)
psd_on_mean = np.mean(psd_on,axis=0)
sem_on = stats.sem(psd_on,axis=0)
sem_off = stats.sem(psd_off,axis=0)

plt.plot(np.log10(f), np.log10(psd_on_mean),'xkcd:gold', label = 'LFP (flash on)')
plt.plot(np.log10(f), np.log10(psd_off_mean),'xkcd:coral', label = 'LFP (flash off)')
#plt.plot(np.log10(f), np.log10(psd_off_mean), label = 'flash off')
plt.title('PSDs of a local field potential and diffusion potentials' )
plt.xlabel('$log_{10}(frequency)$ in Hz')
plt.ylabel('$log_{10}(PSD)$ in (mV)²/Hz')
plt.plot(np.log10(f_Dietzel1982[31:-1]),np.log10(psd_Dietzel1982[31:-1]), 'xkcd:plum', label = 'diffusion')
plt.plot(np.log10(f_Herreras1993[31:-1]), np.log10(psd_Herreras1993[31:-1]), 'xkcd:teal', label = 'SD diffusion')
plt.legend()

plt.savefig('PSD_Gratiy',dpi=500)
plt.show()

sys.exit()


psd_on_mean_ls =np.interp(f_log, f, psd_on_mean)
psd_off_mean_ls =np.interp(f_log, f, psd_off_mean)

sem_on_ls =np.interp(f_log, f, sem_on)
sem_off_ls =np.interp(f_log, f, sem_off)