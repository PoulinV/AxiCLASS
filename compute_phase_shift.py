from classy import Class
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from classy import CosmoComputationError
from time import time
from matplotlib.backends.backend_pdf import PdfPages
from sys import exit
import matplotlib
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


start = time()

# ###############################################################################################################################################
# # REQUIRED / USEFUL FUNCTIONS
# ###############################################################################################################################################

# def is_number(s):
# 	# ---------------------------------- This func checks whether a thing is a number. Found online
# 	try:
# 		float(s)
# 		return True
# 	except ValueError:
# 		pass

# 	try:
# 		import unicodedata
# 		unicodedata.numeric(s)
# 		return True
# 	except (TypeError, ValueError):
# 		pass

# 	return False

# def read_ini_file(inifile, loc = ''):
# 	'''
# 	Function to read ini file and save it in a dictionary that can be passed to classy
# 	Takes the required argument inifile = filename with extension
# 	Takes the optional argument loc = location of your ini file, ending in a '/'
# 	Returns dictionary of everything contained in your ini file
# 	'''
# 	inivals = {}

# 	with open(loc + inifile) as f: # opening initialisation file as f
# 		content = f.readlines() # reading the initialisation file and turning it into a list

# 	q = {} # initialise q as an empty dictionary
# 	for s in content: # iterates over lines in .ini file and saves information, turning numbers into floats from strings
# 		if is_number(s[s.find('=')+2:]):
# 			q[s[:s.find(' =')]] = float(s[s.find('=')+2:])
# 		else:
# 			q[s[:s.find(' =')]] = s[s.find('=')+2:-1]

# 	return q # inivals dict has dict of initial values at key given by 'original'

###############################################################################################################################################
# SET PARAMETER VALUES WHETHER BY READING AN INI FILE OR BY TYPING IN
###############################################################################################################################################

'''
	This is likely the only section you will need to modify
	All parameters defined here
'''

# Functionality for reading in .ini file isn't built in yet.
# Will build that if needed.


params = { 'compute_phase_shift':'yes',
			'output':'tCl,mPk,dTk',
			'start_small_k_at_tau_c_over_tau_h':1e-07,
			'start_large_k_at_tau_h_over_tau_k':1e-07,
			'P_k_max_h/Mpc':5}
			# 'input_verbose':1,
			# 'background_verbose':1,
			# 'thermodynamics_verbose':1,
			# 'perturbations_verbose':1,
			# 'transfer_verbose':1,
			# 'primordial_verbose':1,
			# 'spectra_verbose':1,
			# 'nonlinear_verbose':1,
			# 'lensing_verbose':1,
			# 'output_verbose':1}
# params['100*theta_s'] = 1.042142
params['h'] = 0.67
# params['back_integration_stepsize'] = 5e-4
# params['reio_parametrization'] = 'reio_none'

output_file = 'omcdm_phase_shift_kmax5' # Just name, no extension!!
# file to either store a_c, Om_fld and H0 values computed below, or to load a_c, Om_fld and H0 values from

load_from_file = False # are we loading a_c, Omega_fld and H0 values from a file or computing them right now?

make_plot = True # do you want it to also make a plot or just save the arrays ?
# plot arguments will need to be toggled below under the plot section

N = 10 # Number of bins / resolution

# set sampling boundaries here
# Eps_min = 0.01 # min a_c
# Eps_max = 0.8 # max a_c

om_cdm_min = 0.01
om_cdm_max = 0.3

A = 7./8.*(4./11.)**(4./3.)

# Contours_at = (71.5,73.24,74.98) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
# # 					 # Can define just one value as (5. , )
# # 					 # Leave as None if don't want contours
# # Contours_at = (7,) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
# # 					 # Can define just one value as (5. , )
# # 					 # Leave as None if don't want contours
# Contours_at = (71.5,74.98) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
# Contours_at = (10**3,10**4) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
# 					 # Can define just one value as (5. , )
# 					 # Leave as None if don't want contours
# Contours_at = (7,) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
# 					 # Can define just one value as (5. , )
# 					 # Leave as None if don't want contours

# Contours_zc = (3,4,5) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
Contours_mu = (5,6,7) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma
Contours_alpha = (-1.5,-1,-0.5) # This must be a tuple with a at least one value and a comma. That is, must have at least one comma

# T_b_redshift = 20 # redshift at which we want gas temperature to be plot

###############################################################################################################################################
# RUN CLASS / FETCH OUTPUT
###############################################################################################################################################

if (output_file != None and load_from_file == True):
	exit('Error: Only have functionality for T_b and H_0 built in. Pick one. Type it correctly.')

elif load_from_file == False:
	# if we're calculating things now and not loading from a file

	# Eps_x = np.linspace(Eps_min, Eps_max, N, endpoint = True)
	# print Eps_x

	om_cdm = np.linspace(om_cdm_min, om_cdm_max, N, endpoint = True)

	cosmo = Class() # Create an instance of the CLASS wrapper


	phase_shift = np.zeros((N)) # initialise H0 as an array of zeroes
	amplitude = np.zeros((N)) # initialise H0 as an array of zeroes

	for i in range(N):
		# params['N_ur'] = Eps_x[i]/(A*(1-Eps_x[i]))
		params['omega_cdm'] = om_cdm[i]
		print  om_cdm[i]
		# print Eps_x[i]/(A*(1-Eps_x[i]))
		# params['m_axion'] = 10**mu[i]
		cosmo.empty()
		cosmo.struct_cleanup()
		# ensuring memory isn't being eaten up

		# try to solve with a certain cosmology, no worries if it cannot
		try:
			cosmo.set(params) # Set the parameters to the cosmological code
			cosmo.compute() # solve physics

			phase_shift[i] = (cosmo.phase_shift())
			amplitude[i] =(cosmo.amplitude())
			# zc[i][j] = np.log10(cosmo.zc())
			# print cosmo.phase_shift()

		except CosmoComputationError: # this happens when CLASS fails
			pass # eh, don't do anything


		# print('Nur = %e \t EpsX = %e \t phase_shift = %e \t amplitude = %e\n' %(Eps_x[i]/(A*(1-Eps_x[i])), Eps_x[i],  phase_shift[i], amplitude[i]))
		print('om_cdm = %e \t phase_shift = %e \t amplitude = %e\n' %(om_cdm[i],  phase_shift[i], amplitude[i]))
		# print('fEDE = %e \t mu = %e \t alpha = %.5f \t zc = %.5f\n' %(fEDE[i], mu[j], alpha[i][j], zc[i][j]))

		# # test that stuff is working by plotting the fluid energy density
		# bg = cosmo.get_background()
		# plt.loglog( bg['z'], bg['(.)rho_fld[0]'])
		# plt.show()


###############################################################################################################################################
# PLOT THINGS
###############################################################################################################################################

end = time()
print('\n\nTime taken for everything but plot = %.2f seconds' %(end - start))


if make_plot == True:



	# fig, (ax1,ax2) = plt.subplots() # new plot
	fig, (ax1,ax2) = plt.subplots(1,2, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(12,12))
	# fig_Pk, ax_Pk = plt.subplots()
	# fig.subplots_adjust(hspace=0)
	ax1.tick_params(axis='x',labelsize=23)
	ax1.tick_params(axis='y',labelsize=23)
	ax2.tick_params(axis='x',labelsize=23)
	ax2.tick_params(axis='y',labelsize=23)


	# The following code based on imshow doesn't altogether work
	# Seems like it depends on your python and matplotlib distribution and version
	# # label this object as cax to get colourbar to correspond to it
	# cax = ax.imshow( np.transpose(to_plot), # for some reason, need to transpose H0 to get it to plot correctly
	# 				interpolation='nearest', # can play with various interpolation techniques. May / may not make a difference
	# 				cmap=cm.plasma, # choose colour map
	# 				extent = ( min(z_c_p_1), max(z_c_p_1), min(Omega_fld), max(Omega_fld) ), # place colour map at correct a_c and Om_fld values
	# 				origin = 'lower', # same as transpose, to get pixels to go in the correct places
	# 				aspect = 'auto'
	# 				)
	# ax1.plot(Eps_x,np.arcsin(phase_shift),'r')
	# ax1.plot(Eps_x,0.191*3.14157*Eps_x,'b--')
	# ax2.plot(Eps_x,amplitude,'r')
	ax1.plot(om_cdm,np.arcsin(phase_shift),'r')
	# ax1.plot(Eps_x,0.191*3.14157*Eps_x,'b--')
	# ax2.plot(om_cdm,amplitude,'r')

	ax1.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	ax1.set_ylabel(r'$\Delta\phi$', fontsize = 23)
	# ax2.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	# ax2.set_ylabel(r'Amplitude', fontsize = 23)
	plt.savefig(output_file + '.png',bbox_inches='tight')
	plt.show()


else:
	print('You did not request a plot. Plot is skipped. ')
