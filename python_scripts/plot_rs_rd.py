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
from scipy.interpolate import interp1d

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

#


# fig, (ax1,ax2) = plt.subplots() # new plot
fig, (ax1,ax2,ax3) = plt.subplots(1,3, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(20,15))
# fig, (ax1,ax2) = plt.subplots(1,2, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(30,30))
# fig, (ax1,ax2) = plt.subplots(1,2, sharex=False,gridspec_kw=dict(height_ratios=[1]),figsize=(15,15))
# fig_Pk, ax_Pk = plt.subplots()
# fig.subplots_adjust(hspace=0)
ax1.tick_params(axis='x',labelsize=23)
ax1.tick_params(axis='y',labelsize=23)
ax2.tick_params(axis='x',labelsize=23)
ax2.tick_params(axis='y',labelsize=23)
ax3.tick_params(axis='x',labelsize=23)
ax3.tick_params(axis='y',labelsize=23)
# ax3.tick_params(axis='y',labelsize=23)
# ax1.set_xlim(0,100)
# ax2.set_xlim(0,100)
# ax3.set_xlim(0,100)
ax1.set_xlabel(r'$\Theta_i$', fontsize = 23)
ax1.set_ylabel(r'$r_s^{\rm EDE}/r_s^{\rm LCDM}-1$', fontsize = 23)
ax2.set_xlabel(r'$\Theta_i$', fontsize = 23)
ax2.set_ylabel(r'$(r_d/r_s^{\rm EDE})/(r_d/r_s^{\rm LCDM})-1$', fontsize = 23)
ax3.set_xlabel(r'$\Theta_i$', fontsize = 23)
ax3.set_ylabel(r'$(d_A^{\rm EDE})/(d_A^{\rm LCDM})-1$', fontsize = 23)
# ax3.set_xlabel(r'$z$', fontsize = 23)
# ax3.set_ylabel(r'$\rho$', fontsize = 23)
# ax3.set_xscale('log')
# ax3.set_yscale('log')
# ax2.set_xlabel(r'$\epsilon_X$', fontsize = 23)
###############################################################################################################################################
# RUN CLASS / FETCH OUTPUT
###############################################################################################################################################
log10_axion_ac = -3.55314
n_axion = 3
N = 20
log10_fraction_axion_ac = -0.8788
Theta_initial = np.linspace(1.5, 3, N, endpoint = True)
rs_rec=[]
rd_rec=[]
da_rec=[]
# color=['r','b']
cosmo = Class()
cosmo.set({'h':0.72,
'write thermodynamics':'yes',
'compute damping scale':'yes'})
cosmo.compute() # solve physics
derived = cosmo.get_current_derived_parameters(['z_rec','tau_rec','conformal_age','rs_rec','rd_rec'])
rs_LCDM = int(1000.*derived['rs_rec'])/1000.
rd_LCDM = int(1000.*derived['rd_rec'])/1000.
# label=[r'$\Delta N_{\rm eff} = 0.5$',r'$f_{\rm EDE}(a_c) = 10\%,~a_c = 0.0002$']
for i in range(N):
	print Theta_initial[i]
	# params['Omega_fld']
	# params['N_ur'] = Eps_x[i]/(A*(1-Eps_x[i]))
	cosmo.set({'scf_potential': 'axion',
    'n_axion': n_axion,
    'log10_axion_ac': log10_axion_ac, # Must input log10(axion_ac)
    # log10_fraction_axion_ac': -1.922767 # Must input log10(fraction_axion_ac)
    'log10_fraction_axion_ac': log10_fraction_axion_ac, # Must input log10(fraction_axion_ac)
    # m_axion': 1.811412e+06
    # f_axion': 1
    'scf_parameters':'%.2f,0.0'%(Theta_initial[i]), #phi_i,phi_dot_i //dummy: phi_i will be updated.
    'adptative_stepsize': 500,
    'scf_tuning_index': 0,
    'do_shooting': 'yes',
    'do_shooting_scf': 'yes',
	'h':0.72,
	# 'output':'tCl,lCl,mPk',
    'back_integration_stepsize':1e-3,
    'use_big_theta_scf': 'yes',
    'scf_has_perturbations': 'yes',
    'attractor_ic_scf': 'no',
	'write thermodynamics':'yes',
	'compute damping scale':'yes'})
	cosmo.compute() # solve physics
	derived = cosmo.get_current_derived_parameters(['z_rec','tau_rec','conformal_age','rs_rec','rd_rec','da_rec'])
	#print derived.viewkeys()
	rs_rec.append(derived['rs_rec']/rs_LCDM-1) # round down at 4 digits after coma
	# rd_rec.append(derived['rs_rec']/derived['rd_rec']/(rs_LCDM/rd_LCDM)-1) # round down at 4 digits after coma
	rd_rec.append(derived['rd_rec']/rd_LCDM-1) # round down at 4 digits after coma
	da_rec.append(derived['da_rec']/rd_LCDM-1) # round down at 4 digits after coma

	cosmo.empty()
	cosmo.struct_cleanup()


ax1.plot(Theta_initial,rs_rec)
ax2.plot(Theta_initial,rd_rec)
ax3.plot(Theta_initial,da_rec)
ax1.axvline(2.58,color='red')
ax2.axvline(2.58,color='red')
ax3.axvline(2.58,color='red')
ax1.axvline(2.74,color='magenta')
ax2.axvline(2.74,color='magenta')
ax3.axvline(2.74,color='magenta')
# ax2.plot(Theta_initial,rs_LCDM/rd_LCDM*Theta_initial/Theta_initial)
		# ax3.plot(background['z'],background_density,color=color[i],label=label[i])
# ax1.legend(prop={'size':18},loc='upper right',numpoints=1,frameon=False,handlelength=1.5)

###############################################################################################################################################
# PLOT THINGS
###############################################################################################################################################


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

	# ax1.plot(Eps_x,0.191*3.14157*Eps_x,'b--')
	# ax2.plot(Eps_x,amplitude,'r')
	# ax1.plot(om_cdm,np.arcsin(phase_shift),'r')
	# ax1.plot(Eps_x,0.191*3.14157*Eps_x,'b--')
	# ax2.plot(om_cdm,amplitude,'r')

	# ax1.set_xlabel(r'$\epsilon_X$', fontsize = 23)
	# ax1.set_ylabel(r'$\Delta\phi$', fontsize = 23)

	# ax2.set_ylabel(r'Amplitude', fontsize = 23)
# plt.savefig(output_file + '.png',bbox_inches='tight')
plt.show()
