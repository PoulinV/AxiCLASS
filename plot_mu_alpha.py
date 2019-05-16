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

# Parameters we won't be changing
params = {'scf_potential' : 'axion',
			'output':'tCl',
             'n_axion' : 3,
            'scf_parameters' : '3,0',
			'scf_tuning_index':0,
			'scf_evolve_like_axionCAMB':'no',
			'do_shooting':'yes',
			'do_shooting_scf':'yes',
			'use_big_theta_scf':'yes',
			'scf_has_perturbations':'yes',
			'attractor_ic_scf':'no',
			'adptative_stepsize':100,
			'scf_evolve_as_fluid':'no'}
			# 'threshold_scf_fluid_m_over_H': 1e-3,
			# 'back_integration_stepsize':'1e-3'
			#, 'Omega_many_fld' : 1e-10,
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

params['100*theta_s'] = 1.042142
# params['h'] = 0.67
# params['back_integration_stepsize'] = 5e-4
params['reio_parametrization'] = 'reio_none'

output_file = 'n3_Theta3' # Just name, no extension!!
# file to either store a_c, Om_fld and H0 values computed below, or to load a_c, Om_fld and H0 values from

load_from_file = False # are we loading a_c, Omega_fld and H0 values from a file or computing them right now?

make_plot = False # do you want it to also make a plot or just save the arrays ?
# plot arguments will need to be toggled below under the plot section

N = 10 # Number of bins / resolution

# set sampling boundaries here
fEDE_min = -3 # min a_c
fEDE_max = -0.5 # max a_c
ac_min=-5
ac_max=-3
# mu_min=4
# mu_max=8


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

	fEDE = np.linspace(fEDE_min, fEDE_max, N, endpoint = True)
	# a_c is N log spaced values includig your min and max
	print fEDE

	ac = np.linspace(ac_min, ac_max, N, endpoint = True)
	# mu = np.linspace(mu_min, mu_max, N, endpoint = True)
	# Om_fld is N log spaced values includig your min and max
	print ac

	cosmo = Class() # Create an instance of the CLASS wrapper


	mu = np.zeros((N,N)) # initialise H0 as an array of zeroes
	# zc = np.zeros((N,N)) # initialise H0 as an array of zeroes
	alpha = np.zeros((N,N)) # initialise H0 as an array of zeroes

	for i in range(N):
		params['log10_axion_ac'] = ac[i]
		# params['m_axion'] = 10**mu[i]

		for j in range(N):
			print i,j
			# going over a_c and Om_fld values
			params['log10_fraction_axion_ac'] = fEDE[j]

			cosmo.empty()
			cosmo.struct_cleanup()
			# ensuring memory isn't being eaten up

			# try to solve with a certain cosmology, no worries if it cannot
			try:
				cosmo.set(params) # Set the parameters to the cosmological code
				cosmo.compute() # solve physics

				alpha[i][j] = (cosmo.log10_f_axion())
				mu[i][j] =(cosmo.log10_m_axion())
				# zc[i][j] = np.log10(cosmo.zc())


			except CosmoComputationError: # this happens when CLASS fails
				pass # eh, don't do anything


			print('fEDE = %e \t ac = %e \t alpha = %.5f \t mu = %.5f\n' %(fEDE[i], ac[j], alpha[i][j], mu[i][j]))
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
	fig, (ax1, ax2) = plt.subplots(1,2, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(12,12))
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

	# pcolormesh is fast, simple and accurate
	cax1 = ax1.pcolormesh( ac, fEDE, np.transpose(alpha), # for some reason, need to transpose H0 to get it to plot correctly
					cmap=cm.plasma # choose colour map
					)
	cax2 = ax2.pcolormesh( ac, fEDE, np.transpose(mu), # for some reason, need to transpose H0 to get it to plot correctly
					cmap=cm.plasma # choose colour map
					)
	# cax1 = ax1.pcolormesh( mu, fEDE, np.transpose(alpha), # for some reason, need to transpose H0 to get it to plot correctly
	# 				cmap=cm.plasma # choose colour map
	# 				)
	# cax2 = ax2.pcolormesh( mu, fEDE, np.transpose(zc), # for some reason, need to transpose H0 to get it to plot correctly
	# 				cmap=cm.plasma # choose colour map
	# 				)

	# cax = ax.contourf( z_c_p_1, Omega_fld, np.transpose(to_plot), # for some reason, need to transpose H0 to get it to plot correctly
	# 				cmap=cm.plasma, # choose colour map
	# 				corner_mask = True,
	# 				extent = ( min(z_c_p_1), max(z_c_p_1), min(Omega_fld), max(Omega_fld) ), # place colour map at correct a_c and Om_fld values
	# 				origin = 'lower', # same as transpose, to get pixels to go in the correct places
	# 				aspect = 'auto'
	# 				)


	# add plot of constraints
	# constraints = np.loadtxt('results_ALP_n' + str(params['n_pheno_axion']) + '.dat')
	# ax.loglog(1./constraints[:,0], constraints[:,1], 'c', lw = 2) # 1/stuff because a_c given in constraint .dat files, 'c' cyan is the best showing colour


	CS1 = ax1.contour(ac, fEDE, np.transpose(alpha), Contours_alpha,
			extent = ( min(ac), max(ac), min(fEDE), max(fEDE) ), # place colour map at correct a_c and Om_fld values
			origin = 'lower',
			corner_mask = True,
			colors=('w','w','w'))
	CS2 = ax2.contour(ac, fEDE, np.transpose(mu), Contours_mu,
			extent = ( min(ac), max(ac), min(fEDE), max(fEDE) ), # place colour map at correct a_c and Om_fld values
			origin = 'lower',
			corner_mask = True,
			colors=('w','w','w'))

	# CS1 = ax1.contour(mu, fEDE, np.transpose(alpha), Contours_alpha,
	# 		extent = ( min(mu), max(mu), min(fEDE), max(fEDE) ), # place colour map at correct a_c and Om_fld values
	# 		origin = 'lower',
	# 		corner_mask = True,
	# 		colors=('w','w','w'))
	# CS2 = ax2.contour(mu, fEDE, np.transpose(zc), Contours_zc,
	# 		extent = ( min(mu), max(mu), min(fEDE), max(fEDE) ), # place colour map at correct a_c and Om_fld values
	# 		origin = 'lower',
	# 		corner_mask = True,
	# 		colors=('w','w','w'))

			# alpha = 0.5)
			# cmap = cm.Reds) # Greens seem to work well for the contour lines to be visible
	ax1.clabel(CS1, inline=1, fontsize=10)
	ax2.clabel(CS2, inline=1, fontsize=10)

	ax1.set_title(r'$n = $ %d, $\theta_i = 3$, ${\rm Log}_{10}\alpha$' %params['n_axion'], fontsize = 22)
	ax2.set_title(r'$n = $ %d, $\theta_i = 3$, ${\rm Log}_{10}\mu$' %params['n_axion'], fontsize = 22)
	# ax2.set_title(r'$n = $ %d, $\theta_i = 3$, ${\rm Log}_{10}z_c$' %params['n_axion'], fontsize = 22)

	ax1.set_xlabel(r'${\rm Log}_{10}a_c$', fontsize = 23)
	# ax1.set_xlabel(r'${\rm Log}_{10}\mu$', fontsize = 23)
	ax1.set_ylabel(r'${\rm Log}_{10}f(a_c)$', fontsize = 23)
	ax2.set_xlabel(r'${\rm Log}_{10}a_c$', fontsize = 23)
	# ax2.set_xlabel(r'${\rm Log}_{10}\mu$', fontsize = 23)
	ax2.set_ylabel(r'${\rm Log}_{10}f(a_c)$', fontsize = 23)

	# ax.set_xscale("log")
	# ax.set_yscale("log")
	# plt.xscale('log')
	# plt.yscale('log')


	ax1.set_xlim((ac.min(), ac.max()))
	# ax1.set_xlim((mu.min(), mu.max()))
	ax1.set_ylim((fEDE.min(), fEDE.max()))
	ax2.set_xlim((ac.min(), ac.max()))
	# ax2.set_xlim((mu.min(), mu.max()))
	ax2.set_ylim((fEDE.min(), fEDE.max()))


	# Add colorbar, make sure to specify tick locations to match desired ticklabels
	cbar1= fig.colorbar(cax1) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
	cbar2 = fig.colorbar(cax2) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
	# cbar.ax.set_yticklabels(['%.2f' % H0.min(), '%.2f' % (0.5*(H0.min() + H0.max())), '%.2f' % H0.max()])  # vertically oriented colorbar
	ticklabs = cbar1.ax.get_yticklabels()
	cbar1.ax.set_yticklabels(ticklabs, fontsize=20)
	cbar1.set_label(r'$\alpha$', fontsize = 22)
	ticklabs = cbar2.ax.get_yticklabels()
	cbar2.ax.set_yticklabels(ticklabs, fontsize=20)
	cbar2.set_label(r'$\mu$', fontsize = 22)
	# cbar2.set_label(r'$z_c$', fontsize = 22)
	plt.savefig(output_file + '_mu_alpha.png',bbox_inches='tight')
	# plt.savefig(output_file + '_zc_alpha.png',bbox_inches='tight')
	plt.show()


else:
	print('You did not request a plot. Plot is skipped. ')
