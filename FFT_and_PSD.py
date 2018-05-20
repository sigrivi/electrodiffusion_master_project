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

# I used this program to understand the difference betwen FFT and PSD
# I used this program to make the plots in figure 2.5

Fs = 150
Ts = 1./Fs
#x = np.linspace(0,10,100)
t = np.linspace(0,1,Fs)

y = np.sin(2*np.pi*1*t) + .4*np.sin(2*np.pi*4*t) + .1*np.sin(2*np.pi*32*t)
#y=np.exp(-t)
plt.plot(t,y)
plt.xlabel('t')
plt.ylabel('v(t)')
plt.title('$v(t) = \sin (2\pi t) + 0.4 \sin (8\pi t)+ 0.1 sin (64\pi t)$')
plt.savefig('sine_signal', dpi=500)
plt.show()
#y = np.sin(x)+np.sin(2*x)+np.sin(3*x)
n = len(y)
#plt.plot(t,y)
#plt.show()
fs = 100
fm, FFT = signal.periodogram(y, Fs)
plt.plot(fm,FFT)
plt.xlabel('frequency')
plt.ylabel('power')
plt.title('The power spectrum density of \n $v(t) = \sin (2\pi t) + 0.4 \sin (8\pi t)+ 0.1 sin (64\pi t)$')
plt.savefig('PSD_of_sine', dpi = 500)
#plt.show()
Y = np.fft.fft(y)/n #compute fft and norm the transformed signal

Y = np.abs(Y) # absolute values, use no imaginary frequencies
Y = Y**2
Y = Y[:(n//2)] # use only positive frequencies
Y = Y*2
plt.plot(Y) 
plt.legend()
plt.show()

plt.plot(np.log10(fm),np.log10(FFT))
plt.show()