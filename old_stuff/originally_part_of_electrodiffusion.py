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

# -----------------------------------------------------------------------------
# concentration profiles                          



# 0 Gratiy2017
	c_values = np.asarray(\
		      [0,0.1350844277673544, 0.3827392120075046, 0.48780487804878053, 0.74296435272045, 1.1782363977485926,\
	           1.4859287054409005, 1.5759849906191372, 1.4934333958724202, 1.891181988742964, 1.590994371482176, \
	           1.0281425891181983, 0.6979362101313321, 0.6378986866791744, 0.6378986866791747, 0.6829268292682926,\
	           0.8405253283302062, 0.893058161350844, 0.8930581613508446, 0.8780487804878045, 0.8930581613508448, \
	           0.8255159474671664, 0.908067542213884, 0.8180112570356467, 0.8405253283302068, 0.7804878048780483,0])
	x_values = np.linspace(0,26,num = 27)

	Gratiy2017 = ConcentrationProfile(x_values, c_values, delta_x, 150, 3, 'Gratiy2017')

# 1 Halnes2016
	halnes_delta_c = np.load('halnes_delta_c.npy')
	halnes_x_values = np.linspace(0,14,num = 15)
	Halnes2016 = ConcentrationProfile(halnes_x_values, halnes_delta_c, delta_x, 150, 3, 'Halnes2016')

# 2 Dietzel1982
	Dietzel1982_1 = ConcentrationProfile([0,2,4,6,8,10,16,17], [0,7.5,4.5,3.5,5,3,0.5,0], delta_x, 148,3, 'Dietzel1982_1')

# 3 Nicholson1987
	Nicholson1987 = ConcentrationProfile([0,1,2,3,4,5,6,7], [0, 4.4, 2.7, 1.6, 1., 0.8, 0.7, 0], delta_x, 150, 3, 'Nicholson1987')

# 4 EkstremeGradient

# List of all profiles
	Profiles = [Gratiy2017, Halnes2016, Dietzel1982_1, Nicholson1987]

# choose a profile from the list of profiles
	choose_profile = 1

# save the parameters used
	N_x = Profiles[choose_profile].N_x
	parameters = [N_t, delta_t, N_x, delta_x]
	np.save(Profiles[choose_profile].name +"_parameters.npy", parameters)

# vectors for the axes
	t,x = makeAkses(parameters)
# -----------------------------------------------------------------------------

# initialize ions
	Ions = [Ion(Profiles[choose_profile].c_Na,DNa,zNa,'$Na^+$'),Ion(Profiles[choose_profile].c_K, DK, zK,'$K^+$' ),Ion(Profiles[choose_profile].c_Cl, DCl, zCl,'$Cl^-$' )]


# check electroneutrality
	el_sum = electroneutrality(Ions, N_x, plot = 'true') # plot = 'true' if you want to plot 
	assert np.amax(el_sum) < 1.e-14       # unit test

# plot initial ion concentration
	plotIons(Ions, x, Profiles[choose_profile].name)

# solve the equation
	[sodium, chloride, potassium], Phi_of_t = solveEquation(Ions, lambda_n, N_t, delta_t, N_x, delta_x)

# Phi_of_t is dimensionless, needs to be multiplied with Psi = RT/F = 0.0267 V
# to get Phi_of_t in mV: *1000
	Phi_of_t = Phi_of_t*Psi*1000

	f, psd, location = makePSD(Phi_of_t, N_t, delta_t)
	plt.plot(np.log10(f[1:-1]),np.log10(psd[1:-1]), label=Profiles[choose_profile].name)
	plt.legend()
	plt.savefig('PSD',dpi=225)
	plt.show()

	

# Check electroneutrality
	el_sum = electroneutrality(Ions, N_x, plot = 'true') # true = 'true' if you want to plot
#	assert np.amax(el_sum) < 1.e-14       # unit test

# plot final ion concentration
	plotIons(Ions, x, Profiles[choose_profile].name + '_final')



# save Phi(x,t) 
	np.save( Profiles[choose_profile].name + "_Phi_of_t.npy" , Phi_of_t)

# contour plot of Phi
	X,Y = np.meshgrid(t,x[1:])
	plt.figure()
	cp = plt.contourf(X,Y,Phi_of_t)
	plt.colorbar(cp)
	plt.xlabel('time (s)')
	plt.ylabel('cortical depth (mm)')
	plt.title('$\Phi$ (mV)')
	plt.savefig(Profiles[choose_profile].name + '_Phi_X_T', dpi =225)
	plt.show()

#-----------------------------------------------------------------------------
	sys.exit() 


	plt.plot(x[1:], Phi_of_t[:,1])  # Phi is calculated at half-points => Phi is shorter than x
	plt.ylabel('Phi (mV)')
	plt.xlabel('cortical depth (mm)')
	plt.title('Phi(t=%.1f)'%(delta_t*N_t))
	plt.savefig('phi', dpi=225)
	plt.show()



	sys.exit()
# save to textfile
	#np.ndarray.tofile(Phi_of_t,'phi.txt')

	plt.plot(t,Phi_of_t[N_x//2,:])
	plt.ylabel('Phi (mV)')
	plt.xlabel('time (s) ')
	plt.title('Phi(x=%.1f)'%(N_x/2))
	plt.show()




#	sys.exit()

	plt.plot(x[1:], Phi_of_t[:,N_t-1])  # Phi is calculated at half-points => Phi is shorter than x
	plt.ylabel('Phi (mV)')
	plt.xlabel('cortical depth (mm)')
	plt.title('Phi(t=%.1f)'%(delta_t*N_t))
	plt.savefig('phi', dpi=225)
	plt.show()



	X,Y = np.meshgrid(t,x[1:])
	plt.figure()
	cp = plt.contourf(X,Y,Phi_of_t)
	plt.colorbar(cp)
	plt.xlabel('time (s)')
	plt.ylabel('cortical depth (mm)')
	plt.title('Phi (mV)')
	plt.savefig('Phi_X_T', dpi =225)
	plt.show()



	sys.exit()
	plt.plot(t,Phi_of_t[N_x//2,:])
	plt.ylabel('Phi (mV)')
	plt.xlabel('time (s) ')
	plt.title('Phi(x=%.1f)'%(N_x/2))
	
	plt.show()

#phi_average is the average of phi(t) over the whole cortical depth
#	phi_average = np.zeros(N_t)
#	for i in range(N_t-1):
#		phi_average[i] = np.sum(Phi_of_t[:,i])
#	phi_average = phi_average/(N_x-1)
