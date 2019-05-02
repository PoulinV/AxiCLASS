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
from scipy.interpolate import interp2d
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
			'output':'tCl,mPk',
			# 'use_ppf':'no',
			'start_small_k_at_tau_c_over_tau_h':1e-07,
			'start_large_k_at_tau_h_over_tau_k':1e-07,
			'k_per_decade_for_bao':400,
			'k_per_decade_for_pk':50,
			'P_k_max_h/Mpc':10}
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
A = 7./8.*(4./11.)**(4./3.)*3.046
#
params['h'] = 0.67
params['omega_cdm'] = 0.12
# params['N_ur'] = 0.001
# params['Omega_fld'] = 7./8.*(4./11.)**(4./3.)*3.046*5.38e-5
# params['back_integration_stepsize'] = 5e-4
# params['reio_parametrization'] = 'reio_none'

output_file = 'compute_phase_shift_kmax10' # Just name, no extension!!
# file to either store a_c, Om_fld and H0 values computed below, or to load a_c, Om_fld and H0 values from

load_from_file = False # are we loading a_c, Omega_fld and H0 values from a file or computing them right now?

make_plot = True # do you want it to also make a plot or just save the arrays ?
# plot arguments will need to be toggled below under the plot section

N = 10 # Number of bins / resolution

# set sampling boundaries here
Eps_min = 0.01 # min a_c
Eps_max = 0.8 # max a_c

om_cdm_min = 0.01
om_cdm_max = 0.3
cs2_fld_min = 0.01
cs2_fld_max = 0.5
w_fld_min = 0.01
w_fld_max = 0.5

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
Contours_at = (0.05,0.1,0.2,0.237,0.3)
# T_b_redshift = 20 # redshift at which we want gas temperature to be plot

###############################################################################################################################################
# RUN CLASS / FETCH OUTPUT
###############################################################################################################################################

if (output_file != None and load_from_file == True):
	# if we have a file to load from
	# 1+z_c and Om_fld are in one file, this grabs their lists separately

	Eps_x, cs2_fld = np.loadtxt(output_file + '_epsX_cs2.txt', comments='#', unpack = True)
	N = len(Eps_x)
	phase_shift = np.loadtxt(output_file + '_phase_shift.txt', comments='#')
	# T_b is in another file. We have rows corresponding to 1+z_c values and columns corresponding to Om_fld
	for i in range(N):
		for j in range(N):
			print('Eps_x = %e \t cs2_fld = %e \t phase_shift = %.2f' %(Eps_x[i], cs2_fld[j],  phase_shift[i][j]))


elif load_from_file == False:
	# if we're calculating things now and not loading from a file

	Eps_x = np.linspace(Eps_min, Eps_max, N, endpoint = True)
	# # print Eps_x
	#
	# om_cdm = np.linspace(om_cdm_min, om_cdm_max, N, endpoint = True)
	cs2_fld = np.linspace(cs2_fld_min, cs2_fld_max, N, endpoint = True)
	# print Eps_x

	# w_fld = np.linspace(w_fld_min, w_fld_max, N, endpoint = True)

	if output_file != None:
		# if you want this writted to an output file
		write_in = output_file + '_epsX_cs2.txt'
		with open(write_in, 'w') as f: # creates a new file to write in / overwrites
			# f.write('# 1 + z_c \t\t Omega_fld \n') # column titles
			f.write('# epsX \t\t cs2 \n') # column titles
			for i in range(N):
				f.write(str(Eps_x[i]) + '\t\t' + str(cs2_fld[i]) + '\n') # writes the array to file

		write_in = output_file + '_phase_shift.txt'

	cosmo = Class() # Create an instance of the CLASS wrapper


	phase_shift = np.zeros((N,N)) # initialise H0 as an array of zeroes
	amplitude = np.zeros((N,N)) # initialise H0 as an array of zeroes

	if output_file != None:
		with open(write_in, 'w') as f2: # creates a new file to write in / overwrites
			f2.write('# rows correspond to Eps_x values, columns correspond to cs2 \n') # info on format
			for i in range(N):
				# params['Omega_fld']
				params['N_ur'] = Eps_x[i]/(A*(1-Eps_x[i]))

				# print  om_cdm[i]
				for j in range(N):
					# params['use_ppf'] = 'no'
					# params['w0_fld'] = -0.6
					# params['wa_fld'] = 0

					params['ceff2_ur'] = cs2_fld[j]
					# print Eps_x[i]/(A*(1-Eps_x[i]))
					# params['m_axion'] = 10**mu[i]


					# try to solve with a certain cosmology, no worries if it cannot
					try:
						cosmo.set(params) # Set the parameters to the cosmological code
						cosmo.compute() # solve physics

						phase_shift[i][j] = (cosmo.phase_shift())
						amplitude[i][j] =(cosmo.amplitude())
						# zc[i][j] = np.log10(cosmo.zc())
						# print cosmo.phase_shift()

					except CosmoComputationError: # this happens when CLASS fails
						pass # eh, don't do anything

					cosmo.empty()
					cosmo.struct_cleanup()
					# ensuring memory isn't being eaten up

					print('Eps_x = %e \t cs2_fld = %e \t phase_shift = %.2f' %(Eps_x[i], cs2_fld[j],  phase_shift[i][j]))
					if output_file != None:
						f2.write(str(phase_shift[i][j]) + '\t\t') # write this value of H0 to the file

				if output_file != None:
					f2.write('\n') # new line after all columns for this row are done

			# print('Nur = %e \t EpsX = %e \t phase_shift = %e \t amplitude = %e\n' %(Eps_x[i]/(A*(1-Eps_x[i])), Eps_x[i],  phase_shift[i], amplitude[i]))
			# print('om_cdm = %e \t Neff = %e \t phase_shift = %e \t amplitude = %e\n' %(om_cdm[i],Eps_x[j]/(A*(1-Eps_x[j])),  phase_shift[i][j], amplitude[i][j]))
			# print('cs2_fld = %e \t Eps_x = %e \t phase_shift = %e \t amplitude = %e\n' %(cs2_fld[j], Eps_x[i],  phase_shift[i][j], amplitude[i][j]))
			# print('fEDE = %e \t mu = %e \t alpha = %.5f \t zc = %.5f\n' %(fEDE[i], mu[j], alpha[i][j], zc[i][j]))
	else:
		"please enter a file name"
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
	# fig, (ax1, ax2) = plt.subplots(1,2, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(12,12))
	fig, (ax1) = plt.subplots(1,1, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(12,12))
	# fig_Pk, ax_Pk = plt.subplots()
	# fig.subplots_adjust(hspace=0)
	ax1.tick_params(axis='x',labelsize=23)
	ax1.tick_params(axis='y',labelsize=23)
	# ax2.tick_params(axis='x',labelsize=23)
	# ax2.tick_params(axis='y',labelsize=23)

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
	# cax1 = ax1.pcolormesh( Eps_x, om_cdm, np.transpose(phase_shift), # for some reason, need to transpose H0 to get it to plot correctly
	# 				cmap=cm.plasma # choose colour map
	# 				)
	# cax2 = ax2.pcolormesh( Eps_x, om_cdm, np.transpose(amplitude), # for some reason, need to transpose H0 to get it to plot correctly
	# 				cmap=cm.plasma # choose colour map
	# 				)
	#
	# ax1.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	# # ax1.set_xlabel(r'${\rm Log}_{10}\mu$', fontsize = 23)
	# ax1.set_ylabel(r'$\omega_{\rm cdm}$', fontsize = 23)
	# ax2.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	# # ax2.set_xlabel(r'${\rm Log}_{10}\mu$', fontsize = 23)
	# ax2.set_ylabel(r'$\omega_{\rm cdm}$', fontsize = 23)

	cax1 = ax1.pcolormesh( Eps_x, cs2_fld, np.transpose(phase_shift), # for some reason, need to transpose H0 to get it to plot correctly
					cmap=cm.plasma # choose colour map
					)
	# cax2 = ax2.pcolormesh( Eps_x, cs2_fld, np.transpose(amplitude), # for some reason, need to transpose H0 to get it to plot correctly
	# 				cmap=cm.plasma # choose colour map
	# 				)

	if Contours_at != None:
		CS = ax1.contour(Eps_x, cs2_fld, np.transpose(phase_shift), Contours_at,
				extent = ( min(Eps_x), max(Eps_x), min(cs2_fld), max(cs2_fld) ), # place colour map at correct a_c and Om_fld values
				origin = 'lower',
				corner_mask = True,
				colors=('w','w','w'))
		ax1.clabel(CS, inline=1, fontsize=15)
	# X,Y = np.meshgrid(Eps_x, 0.33*cs2_fld/cs2_fld)
	# CS1 = ax1.contour(X, Y, np.sin(0.191*3.14157*X)*Y,colors=('w'))
   	# ax1.plot3D(Eps_x,0.33,np.arcsin(phase_shift),'w')
   	# ax1.plot3D(Eps_x,0.33,np.sin(0.191*3.14157*Eps_x),'w--')
	ax1.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	# ax1.set_xlabel(r'${\rm Log}_{10}\mu$', fontsize = 23)
	ax1.set_ylabel(r'$c_s^2$', fontsize = 23)
	# ax2.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	# # ax2.set_xlabel(r'${\rm Log}_{10}\mu$', fontsize = 23)
	# ax2.set_ylabel(r'$c_s^2$', fontsize = 23)

	# ax.set_xscale("log")
	# ax.set_yscale("log")
	# plt.xscale('log')
	# plt.yscale('log')


	ax1.set_xlim((Eps_x.min(), Eps_x.max()))
	# ax1.set_xlim((mu.min(), mu.max()))
	ax1.set_ylim((cs2_fld.min(), cs2_fld.max()))
	# ax2.set_xlim((Eps_x.min(), Eps_x.max()))
	# # ax2.set_xlim((mu.min(), mu.max()))
	# ax2.set_ylim((cs2_fld.min(), cs2_fld.max()))


	# Add colorbar, make sure to specify tick locations to match desired ticklabels
	cbar1= fig.colorbar(cax1) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
	# cbar2 = fig.colorbar(cax2) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
	# cbar.ax.set_yticklabels(['%.2f' % H0.min(), '%.2f' % (0.5*(H0.min() + H0.max())), '%.2f' % H0.max()])  # vertically oriented colorbar
	ticklabs = cbar1.ax.get_yticklabels()
	cbar1.ax.set_yticklabels(ticklabs, fontsize=20)
	cbar1.set_label(r'sin($\Delta\phi$)', fontsize = 22)
	# ticklabs = cbar2.ax.get_yticklabels()
	# cbar2.ax.set_yticklabels(ticklabs, fontsize=20)
	# cbar2.set_label('Amplitude', fontsize = 22)
	phase_shift_at = interp2d(Eps_x,cs2_fld,phase_shift)
	ax1.text(0.408,0.335,'%.3f'%(np.arcsin(phase_shift_at(0.408,0.33))),color='w',fontsize=15)
	ax1.plot(0.408,0.33,'*',markerfacecolor='white', markersize=15)
	# ax1.plot(Eps_x,0.33,'*',markerfacecolor='white', markersize=15)
	plt.savefig(output_file + '.png',bbox_inches='tight')
	plt.show()


else:
	print('You did not request a plot. Plot is skipped. ')
