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

# params['h'] = 0.67
# params['back_integration_stepsize'] = 5e-4
# params['reio_parametrization'] = 'reio_none'

output_file = 'n2_check_highTheta' # Just name, no extension!!
# file to either store a_c, Om_fld and H0 values computed below, or to load a_c, Om_fld and H0 values from

load_from_file = True # are we loading a_c, Omega_fld and H0 values from a file or computing them right now?

make_plot = False # do you want it to also make a plot or just save the arrays ?
# plot arguments will need to be toggled below under the plot section

N = 10 # Number of bins / resolution

# set sampling boundaries here
fEDE_min = -1 # min a_c
fEDE_max = -0.5 # max a_c
ac_min=-4.5
ac_max=-3
# mu_min=4
# mu_max=8
n_axion = 2.

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
theta_bug = []
theta_good = []
ac_bug = []
ac_good = []
ratio_max = []
security_small_Omega_scf = [1e-3,1e-3]
###############################################################################################################################################
# RUN CLASS / FETCH OUTPUT
###############################################################################################################################################

if (output_file != None and load_from_file == True):
	# if we have a file to load from
	# 1+z_c and Om_fld are in one file, this grabs their lists separately

	ac, Theta_initial = np.loadtxt(output_file + '_ac_Theta_initial.txt', comments='#', unpack = True)
	N = len(ac)
	ratio_max = np.loadtxt(output_file + '_ratio_max.txt', comments='#')
	# T_b is in another file. We have rows corresponding to 1+z_c values and columns corresponding to Om_fld
	for i in range(N):
		for j in range(N):
			print('ac = %e \t Theta_initial = %e \t ratio_max = %.2f' %(ac[i], Theta_initial[j], ratio_max[i][j]))


elif load_from_file == False:
	# if we're calculating things now and not loading from a file
	write_in = output_file + '_ac_Theta_initial.txt'
	Theta_initial = np.linspace(0.1, 3, N, endpoint = True)
	# a_c is N log spaced values includig your min and max
	print Theta_initial

	ac = np.linspace(ac_min, ac_max, N, endpoint = True)
	print ac

	with open(write_in, 'w') as f: # creates a new file to write in / overwrites
		f.write('# ac \t\t Theta_initial \n') # column titles
		for i in range(N):
			f.write(str(ac[i]) + '\t\t' + str(Theta_initial[i]) + '\n') # writes the array to file


	# mu = np.linspace(mu_min, mu_max, N, endpoint = True)
	# Om_fld is N log spaced values includig your min and max
	ratio_max = np.zeros((N,N)) # initialise H0 as an array of zeroes
	write_in = output_file + '_ratio_max.txt'
	cosmo = Class() # Create an instance of the CLASS wrapper
	with open(write_in, 'w') as f: # creates a new file to write in / overwrites
		f.write('# rows correspond to ac values, columns correspond to Theta_initial \n') # info on format

		for i in range(N):
			print ac[i]
			f.write('\n') # new line after all columns for this row are done
			# params['m_axion'] = 10**mu[i]
			for j in range(N):
				ratio = []
				# print i,j
				# going over a_c and Om_fld values
				print Theta_initial[j]
				for k in range(2):
					if k == 0:
						Theta_initial_k = 2.
					else:
						Theta_initial_k = Theta_initial[j]
					params={'scf_potential': 'axion',
				    'n_axion': n_axion,
					'scf_parameters':'%.2f,0.0'%(Theta_initial_k),
					'log10_axion_ac':ac[i],
				    'log10_fraction_axion_ac': -1, # Must input log10(fraction_axion_ac)
				    'adptative_stepsize': 100,
				    'scf_tuning_index': 0,
				    'do_shooting': 'yes',
				    'do_shooting_scf': 'yes',
					'h':0.72,
					'omega_b':0.02225,
					'omega_cdm':0.1198,
					'ln10^{10}A_s':3.094,
					'n_s':0.9645,
					'tau_reio':0.079,
					'output':'tCl',
					'security_small_Omega_scf':security_small_Omega_scf[k],
					# 'output':'tCl,lCl,mPk',
				    'use_big_theta_scf': 'yes',
				    'scf_has_perturbations': 'yes',
				    'attractor_ic_scf': 'no',
					'write thermodynamics':'yes',
					'compute damping scale':'yes'}
				# ensuring memory isn't being eaten up

					# try to solve with a certain cosmology, no worries if it cannot
					try:
						cosmo.set(params) # Set the parameters to the cosmological code
						cosmo.compute() # solve physics
						clM = cosmo.raw_cl(2500)


						if k == 0:
							ll_1 = clM['ell'][2:]
							clTT_1 = clM['tt'][2:]
						if k == 1:
							ll_2 = clM['ell'][2:]
							clTT_2 = clM['tt'][2:]
							for l in range(len(ll_1)):
								# print ll_1[l]
								sigma_CL = np.sqrt(2.0/(2.0*ll_1[l]+1.0))
								ratio.append(abs(clTT_2[l]/clTT_1[l]-1.0)/sigma_CL)
								# print l, ratio[l]
							ratio_max[i][j] = min(max(ratio),1)
							f.write(str(ratio_max[i][j]) + '\t\t') # write this value of H0 to the file
							print "ac %f theta_i %f ratio %f"%(ac[i],Theta_initial[j],max(ratio))
					except CosmoComputationError: # this happens when CLASS fails
						pass # eh, don't do anything


				# print('fEDE = %e \t ac = %e \t alpha = %.5f \t mu = %.5f\n' %(fEDE[i], ac[j], alpha[i][j], mu[i][j]))
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

fig, (ax1) = plt.subplots() # new plot
# pcolormesh is fast, simple and accurate
print ratio_max
cax1 = ax1.pcolormesh( ac, Theta_initial, np.transpose(np.log10(ratio_max)), # for some reason, need to transpose H0 to get it to plot correctly
				cmap=cm.plasma # choose colour map
				)
cbar1= fig.colorbar(cax1) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above

ax1.set_ylabel(r"$\Theta_i$", fontsize=16)
ax1.set_xlabel(r"$a_c$", fontsize=16)
ax1.set_title(r"n=%.1f"%(n_axion), fontsize=16)
# plt.savefig(output_file + '_zc_alpha.png',bbox_inches='tight')
plt.show()
