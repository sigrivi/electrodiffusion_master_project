import numpy as np
import matplotlib.pyplot as plt
from scipy import io
import sys
from scipy import signal

from eldiff_ions import solveEquation #(Ions,lambda_n, N_t, delta_t, N, delta_x): return(Ions, Phi_of_t)
from eldiff_ions import Ion # (self, c_init, D, z, name ):
from eldiff_ions import integrate
from eldiff_ions import ConcentrationProfile # (values, c_values, delta_x, Na_base, K_base, name)
from eldiff_ions import plotIons #(Ions, x, filename):
from eldiff_ions import makePSD #(Phi_of_t, N_t, delta_t): return(f, psd_max, location)
from eldiff_ions import electroneutrality #(Ions, N, plot = 'false' ):
from eldiff_ions import makeAkses # (parameters):return(t,x)

if __name__=="__main__":



	N_t = 100000          # t_final = N_t * delta_t
	delta_t = 1/1000      # delta_t in seconds
	delta_x = 1/10000     # delta_x i meters

# valence: 
	zNa = 1                                         
	zCl = -1
	zK = 1

# diffusion constants: 
	DNa = 1.33e-9  
	DCl = 2.03e-9
	DK = 1.96e-9

# tortuosity
	lambda_n = 1.6 

# Phi is dimensionless, needs to be multiplied with Psi = RT/F = 0.0267 V 
	Psi = 8.3144598*310/(96485.333)

# 0 Gratiy2017 (baseline: Na:150, K:3)
	c_values = .001*np.asarray(\
		      [0,0.1350844277673544, 0.3827392120075046, 0.48780487804878053, 0.74296435272045, 1.1782363977485926,\
	           1.4859287054409005, 1.5759849906191372, 1.4934333958724202, 1.891181988742964, 1.590994371482176, \
	           1.0281425891181983, 0.6979362101313321, 0.6378986866791744, 0.6378986866791747, 0.6829268292682926,\
	           0.8405253283302062, 0.893058161350844, 0.8930581613508446, 0.8780487804878045, 0.8930581613508448, \
	           0.8255159474671664, 0.908067542213884, 0.8180112570356467, 0.8405253283302068, 0.7804878048780483,0])
	x_values = np.linspace(0,26,num = 27)
	N_x = len(c_values)

	Ions1 = [Ion(np.ones(N_x)*.15 - c_values, DNa, zNa, '$Na^+$'), Ion(np.ones(N_x)*.003 + c_values, DK, zK, '$K^+$'), \
     		        Ion(np.ones(N_x)*.153, DCl, zCl, '$Cl^-$' )]
	Ions2 = [Ion(np.ones(N_x)*.15 - c_values/2, DNa, zNa, '$Na^+$'), Ion(np.ones(N_x)*.003 + c_values, DK, zK, '$K^+$'), \
     		        Ion(np.ones(N_x)*.153 + c_values/2, DCl, zCl, '$Cl^-$' )]
	Ions3 = [Ion(np.ones(N_x)*.15 , DNa, zNa, '$Na^+$'), Ion(np.ones(N_x)*.003 + c_values, DK, zK, '$K^+$'), \
     		        Ion(np.ones(N_x)*.153 + c_values, DCl, zCl, '$Cl^-$' )]    		        
	el_sum = electroneutrality(Ions2, N_x)


#	plotIons(Ions2, x_values, 'ions')
	Ions11, Phi_of_t_1 = solveEquation(Ions1, lambda_n, N_t, delta_t, N_x, delta_x) #(Ions,lambda_n, N_t, delta_t, N, delta_x): return(Ions, Phi_of_t)
	Ions21, Phi_of_t_2 = solveEquation(Ions2, lambda_n, N_t, delta_t, N_x, delta_x)
	Ions31, Phi_of_t_3 = solveEquation(Ions3, lambda_n, N_t, delta_t, N_x, delta_x)
	f1, psd1, loc1 = makePSD(Phi_of_t_1, N_t, delta_t)
	f2, psd2, loc2 = makePSD(Phi_of_t_2, N_t, delta_t)
	f3, psd3, loc3 = makePSD(Phi_of_t_3, N_t, delta_t)

	diff_1_2 = np.mean(np.log10(psd1) - np.log10(psd2))
	diff_1_3 = np.mean(np.log10(psd1) - np.log10(psd3))
	plt.figure()
	plt.title('PSD of the diffusion potential' )
	plt.xlabel('$log_{10}(frequency)$ in Hz')
	plt.ylabel('$log_{10}(PSD)$ in (mV)²/Hz')
	plt.plot(np.log10(f1[1:-1]), np.log10(psd1[1:-1]), label = '$\Delta K = -\Delta Na$')
	plt.plot(np.log10(f1[1:-1]), np.log10(psd2[1:-1]), label = '$\Delta K = -0.5 \Delta Na$')
	plt.plot(np.log10(f1[1:-1]), np.log10(psd3[1:-1]), label = '$\Delta K = \Delta Cl$')

	plt.legend()
	plt.savefig('PSD', dpi = 500)
	sys.exit()
