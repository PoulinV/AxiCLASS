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
from mpl_toolkits.mplot3d import Axes3D

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


start = time()

'''
	This is likely the only section you will need to modify
	All parameters defined here
'''

# Functionality for reading in .ini file isn't built in yet.
# Will build that if needed.

# Parameters we won't be changing
params = {
		  'output':'tCl,pCl,lCl,mPk',
          'lensing':'yes'
          }
params['do_shooting'] = 'no'
params['do_shooting_scf'] = 'no'
# params['back_integration_stepsize'] =1e-4
params['use_big_theta_scf'] = 'yes'
params['scf_has_perturbations'] = 'yes'
params['attractor_ic_scf'] = 'no'
params['scf_potential'] = 'axion'
params['n_axion'] = 2
params['omega_b'] = 0.022383
params['omega_cdm'] = 0.12011
params['tau_reio'] = 0.0543
params['ln10^{10}A_s'] = 3.0448
params['n_s'] = 0.96605
params['h'] = 0.7

#
# params['input_verbose'] =1
# params['background_verbose'] = 3
# params['thermodynamics_verbose'] = 1
# params['perturbations_verbose'] = 1
# params['transfer_verbose'] = 1
# params['primordial_verbose'] = 1
# params['spectra_verbose'] = 1
# params['nonlinear_verbose'] = 1
# params['lensing_verbose'] = 1
# params['output_verbose'] = 1
#


N=5

params['scf_tuning_index'] = 0

params['scf_evolve_as_fluid'] = 'no' ##set to yes for fluid approximation
params['scf_evolve_like_axionCAMB'] = 'no'
cosmo = Class() # Create an instance of the CLASS wrapper
maxion_table = np.logspace(1, 10, N, endpoint = True)
faxion_table = np.logspace(-2, 2, N, endpoint = True)
theta_initial_table = np.linspace(0.0001, 3.14, N, endpoint = True)

# good_maxion = np.zeros((N,N,N)) # initialise H0 as an array of zeroes
# good_faxion = np.zeros((N,N,N)) # initialise H0 as an array of zeroes
# good_theta = np.zeros((N,N,N)) # initialise H0 as an array of zeroes
good_maxion=[]
good_faxion=[]
good_theta=[]
bad_maxion=[]
bad_faxion=[]
bad_theta=[]
good_fede = []
good_zc = []
for i in range(N):
	params['m_axion'] = maxion_table[i]
	for j in range(N):
		params['f_axion'] = faxion_table[j]
		for k in range(N):
			params['scf_parameters'] = '%.5f,0.0'%theta_initial_table[k]
			print 'try:', maxion_table[i],faxion_table[j],theta_initial_table[k]
			# try to solve with a certain cosmology, no worries if it cannot
			try:
				cosmo.set(params) # Set the parameters to the cosmological code
				cosmo.compute() # solve physics
				print 'fac:', cosmo.fEDE(),'log10(zc):', cosmo.zc()
				print 'good!', maxion_table[i],faxion_table[j],theta_initial_table[k]
				# good_maxion[i][j][k]=maxion_table[i]
				# good_faxion[i][j][k]=faxion_table[j]
				# good_theta[i][j][k]=theta_initial_table[k]
				good_maxion.append(np.log10(maxion_table[i]))
				good_faxion.append(np.log10(faxion_table[j]))
				good_theta.append(theta_initial_table[k])
				good_fede.append(np.log10(cosmo.fEDE()))
				good_zc.append(cosmo.zc())
			except CosmoComputationError: # this happens when CLASS fails
				# print CosmoComputationError
				print 'this one was bad...'
				bad_maxion.append(np.log10(maxion_table[i]))
				bad_faxion.append(np.log10(faxion_table[j]))
				bad_theta.append(theta_initial_table[k])
				# good_maxion[i][j][k]=-1
				# good_faxion[i][j][k]=-1
				# good_theta[i][j][k]=-1
				pass # eh, don't do anything

			cosmo.empty()
			cosmo.struct_cleanup()
# good_maxion = np.ma.masked_where( (good_maxion < 0), good_maxion) # add a mask to exclude all values <= 0
# good_faxion = np.ma.masked_where( (good_faxion < 0), good_faxion) # add a mask to exclude all values <= 0
# good_theta = np.ma.masked_where( (good_theta < 0), good_theta) # add a mask to exclude all values <= 0

fig1 = plt.figure()
fig2 = plt.figure()
ax = fig1.add_subplot(111, projection='3d')
ax2 = fig2.add_subplot(111)
ax.tick_params(axis='x',labelsize=23)
ax.tick_params(axis='y',labelsize=23)
ax.tick_params(axis='z',labelsize=23)

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
print good_maxion
print good_faxion
print good_theta

cax=ax.scatter(good_maxion, good_faxion, good_theta,c='b')
cax=ax.scatter(bad_maxion, bad_faxion, bad_theta,c='r')
fax=ax2.scatter(good_zc,good_fede)
# cax=ax.scatter(maxion_table, faxion_table, theta_initial_table)
# cax=ax.scatter(maxion_table, faxion_table)


# cax = ax.pcolormesh( maxion_table, faxion_table, good_theta, # for some reason, need to transpose H0 to get it to plot correctly
# 				cmap=cm.plasma # choose colour map
# 				)
# plt.xscale('log')
# plt.yscale('log')

ax.set_xlabel(r'Log10($\mu$)', fontsize = 23)
ax.set_ylabel(r'Log10($\alpha$)', fontsize = 23)
ax.set_zlabel(r'$\theta_i$', fontsize = 23)

ax2.tick_params(axis='x',labelsize=23)
ax2.tick_params(axis='y',labelsize=23)
ax2.set_xlabel(r'Log10($z_c$)', fontsize = 23)
ax2.set_ylabel(r'Log10($f_{\rm ede}$)', fontsize = 23)
# cbar.set_label(r'$\theta_i$', fontsize = 22)
plt.show()
