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

# Functionality for reading in .ini file isn't built in yet.
# Will build that if needed.

# Parameters we won't be changing

# params['h'] = 0.67
# params['back_integration_stepsize'] = 5e-4
# params['reio_parametrization'] = 'reio_none'

# file to either store a_c, Om_fld and H0 values computed below, or to load a_c, Om_fld and H0 values from

load_from_file = True  # are we loading a_c, Omega_fld and H0 values from a file or computing them right now?

make_plot = False # do you want it to also make a plot or just save the arrays ?
# plot arguments will need to be toggled below under the plot section

N1 = 10 # Number of bins / resolution
N2 = 10
# set sampling boundaries here
fEDE_min = -1 # min a_c
fEDE_max = -0.5 # max a_c
ac_min=-4.5
ac_max=-3
# mu_min=4
# mu_max=8
n_axion = 2.2
# output_file = 'n%.0f'%(n_axion) # Just name, no extension!!
output_file = 'n2p2'# Just name, no extension!!
wn = (n_axion-1)/(n_axion+1)

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
security_small_Omega_scf = [-1,1e-3]
###############################################################################################################################################
# RUN CLASS / FETCH OUTPUT
###############################################################################################################################################
h = 0.72
_c_ = 299792.458 #km/s
K_star = 0.05
n_s = 0.9645
A_s = 2.22343e-9
TEST_PLOT = False
if (output_file != None and load_from_file == True):
	# if we have a file to load from
	# 1+z_c and Om_fld are in one file, this grabs their lists separately

	ac = np.loadtxt(output_file + '_ac.txt', comments='#', unpack = True)
	Theta_initial = np.loadtxt(output_file + '_Theta_initial.txt', comments='#', unpack = True)
	N1 = len(ac)
	N2 = len(Theta_initial)
	Omega_nl_scf = np.loadtxt(output_file + '_Omega_znl.txt', comments='#')
	z_nl_scf = np.loadtxt(output_file + '_znl.txt', comments='#')
	k_nl_scf = np.loadtxt(output_file + '_knl.txt', comments='#')
	# T_b is in another file. We have rows corresponding to 1+z_c values and columns corresponding to Om_fld
	# for i in range(N1):
	# 	for j in range(N2):
	# 		print('ac = %e \t Theta_initial = %e \t Omega_nl_scf = %.5f znl = %.5f' %(ac[i], Theta_initial[j], Omega_nl_scf[i][j], z_nl_scf[i][j]))


elif load_from_file == False:
    # if we're calculating things now and not loading from a file
    write_in_1 = output_file + '_ac.txt'
    write_in_2 = output_file + '_Theta_initial.txt'
    Theta_initial = np.linspace(0.1, 3, N2, endpoint = True)
    # a_c is N log spaced values includig your min and max
    print Theta_initial

    ac = np.linspace(ac_min, ac_max, N1, endpoint = True)
    print ac

    with open(write_in_1, 'w') as f: # creates a new file to write in / overwrites
    	f.write('# ac \n') # column titles
    	for i in range(N1):
    		f.write(str(ac[i]) + '\n') # writes the array to file
    with open(write_in_2, 'w') as f: # creates a new file to write in / overwrites
    	f.write('# Theta_initial \n') # column titles
    	for i in range(N2):
    		f.write(str(Theta_initial[i]) + '\n') # writes the array to file



    # mu = np.linspace(mu_min, mu_max, N, endpoint = True)
    # Om_fld is N log spaced values includig your min and max
    k_nl_scf = np.zeros((N1,N2)) # initialise H0 as an array of zeroes
    z_nl_scf = np.zeros((N1,N2)) # initialise H0 as an array of zeroes
    Omega_nl_scf = np.zeros((N1,N2)) # initialise H0 as an array of zeroes
    write_in_1 = output_file + '_Omega_znl.txt'
    write_in_2 = output_file + '_znl.txt'
    write_in_3 = output_file + '_knl.txt'
    cosmo = Class() # Create an instance of the CLASS wrapper
    f1 = open(write_in_1, 'w')  # creates a new file to write in / overwrites
    f1.write('# rows correspond to ac values, columns correspond to Theta_initial \n') # info on format
    f2 = open(write_in_2, 'w')  # creates a new file to write in / overwrites
    f2.write('# rows correspond to ac values, columns correspond to Theta_initial \n') # info on format
    f3 = open(write_in_3, 'w')  # creates a new file to write in / overwrites
    f3.write('# rows correspond to ac values, columns correspond to Theta_initial \n') # info on format

    for i in range(N1):
        print ac[i]
        f1.write('\n') # new line after all columns for this row are done
        f2.write('\n') # new line after all columns for this row are done
        f3.write('\n') # new line after all columns for this row are done
        # params['m_axion'] = 10**mu[i]
        for j in range(N2):
            ratio = []
            # print i,j
            # going over a_c and Om_fld values
            print Theta_initial[j]
            #first iteration: extract k_res
            params={'scf_potential': 'axion',
            'n_axion': n_axion,
            'scf_parameters':'%.2f,0.0'%(Theta_initial[j]),
            'log10_axion_ac':ac[i],
            'log10_fraction_axion_ac': -1, # Must input log10(fraction_axion_ac)
            'adptative_stepsize': 100,
            'scf_tuning_index': 0,
            'do_shooting': 'yes',
            'do_shooting_scf': 'yes',
            'h':h,
            'omega_b':0.02225,
            'omega_cdm':0.1198,
            'A_s':A_s,
            'n_s':n_s,
            'tau_reio':0.079,
            'output':'tCl,lCl,mPk,dTk,vTk',
            'use_big_theta_scf': 'yes',
            'scf_has_perturbations': 'yes',
            'attractor_ic_scf': 'no',
            'write background':'yes',
            'compute damping scale':'yes'}
            try:
                cosmo.set(params) # Set the parameters to the cosmological code

                cosmo.compute() # solve physics
                m_axion = 10**cosmo.log10_m_axion()*h*100/_c_
                print "m_axion %f"%(m_axion)
                # k_res_analytical = 0.79*m_axion*Theta_initial[j]*(10**ac[i])
                k_res_analytical = 0.79*m_axion*Theta_initial[j]*(10**ac[i])
                one_time = cosmo.get_transfer(0) # transfer functions at each time tau
                k = one_time['k (h/Mpc)']
                delta_phi_scf = np.abs(one_time['delta_phi_scf']) ##evaluated today
                k_of_delta = interp1d(delta_phi_scf,k)
                background = cosmo.get_background() # load background table
                background_tau = background['conf. time [Mpc]'] # read conformal times in background table
                background_z = background['z'] # read redshift
                background_Om_scf = background['(.)Omega_scf'] # read redshift
                background_z_at_tau = interp1d(background_tau,background_z)
                phi_enveloppe = Theta_initial[j]*10**cosmo.log10_f_axion()*(2/(10**ac[i]))**(-3.*(1+wn)/2/n_axion) ##evaluated today
                delta_max = np.sqrt(delta_phi_scf*delta_phi_scf*A_s*(k/K_star)**(n_s-1))/(phi_enveloppe)
                print max(delta_max), np.where(delta_max==max(delta_max))
                k_res_numerical = k[np.where(delta_max==max(delta_max))]*h ##evaluated today
                print "k_res analytical %f Mpc-1, k_res numerical %f Mpc-1"%(k_res_analytical,k_res_numerical)
            except CosmoComputationError: # this happens when CLASS fails
                print "bug!"
                k_res_numerical = 0
                k_nl_scf[i][j] = 100 #arbitrary large number: bug
                z_nl_scf[i][j] = 0 #arbitrary large number: bug
                Omega_nl_scf[i][j] = 1e-30 #arbitrary small number

            ##2em iteration: find z_nl at which k_res becomes non linear
            cosmo.empty()
            cosmo.struct_cleanup()
            if k_res_numerical != 0:
                print k_res_numerical
                params={'scf_potential': 'axion',
                'n_axion': n_axion,
                'scf_parameters':'%.2f,0.0'%(Theta_initial[j]),
                'log10_axion_ac':ac[i],
                'log10_fraction_axion_ac': -1, # Must input log10(fraction_axion_ac)
                'adptative_stepsize': 100,
                'scf_tuning_index': 0,
                'do_shooting': 'yes',
                'do_shooting_scf': 'yes',
                'h':h,
                'omega_b':0.02225,
                'omega_cdm':0.1198,
                'A_s':A_s,
                'n_s':n_s,
                'tau_reio':0.079,
                'k_output_values':k_res_numerical[0],
                'output':'tCl,lCl,mPk,dTk,vTk',
                'use_big_theta_scf': 'yes',
                'scf_has_perturbations': 'yes',
                'attractor_ic_scf': 'no',
                'write background':'yes',
                'compute damping scale':'yes'}
                try:
                    cosmo.set(params) # Set the parameters to the cosmological code
                    cosmo.compute() # solve physics

                    all_k = cosmo.get_perturbations()  # this potentially constains scalars/tensors and all k values
                    one_k = all_k['scalar'][0]     # this contains only the scalar perturbations for the requested k values
                    #
                    background = cosmo.get_background() # load background table
                    tau = one_k['tau [Mpc]']
                    delta_phi = one_k['delta_phi_scf']
                    background_tau = background['conf. time [Mpc]'] # read conformal times in background table
                    background_z = background['z'] # read redshift
                    background_Om_scf = background['(.)Omega_scf'] # read redshift
                    background_z_at_tau = interp1d(background_tau,background_z)
                    fudge = 1.6
                    phi_enveloppe = Theta_initial[j]*10**cosmo.log10_f_axion()*(fudge/(10**ac[i]*background_z_at_tau(tau)))**(-3.*(1+wn)/2/n_axion)
                    background_Om_scf_at_z = interp1d(background_z,background_Om_scf)
                    dimensionless_power_spectrum = np.sqrt(delta_phi*delta_phi*A_s*(k_res_numerical[0]/K_star)**(n_s-1))/(phi_enveloppe)
                    tau_at_delta_phi_over_phi = interp1d(dimensionless_power_spectrum,tau)
                    if max(1,max(abs(dimensionless_power_spectrum))) > 1:
                        k_nl_scf[i][j] = k_res_numerical
                        z_nl_scf[i][j] = background_z_at_tau(tau_at_delta_phi_over_phi(1))
                        Omega_nl_scf[i][j] = background_Om_scf_at_z(z_nl_scf[i][j])
                        ##FOR TESTING PURPOSE
                        if TEST_PLOT is True:
                            fig, (ax1) = plt.subplots() # new plot

                            ax1.plot(tau,dimensionless_power_spectrum)
                            ax1.axvline(tau_at_delta_phi_over_phi(1))
                            ax1.set_xscale('log')
                            ax1.set_yscale('log')
                            plt.show()
                    else:
                        k_nl_scf[i][j] = k_res_numerical #arbitrary large number: no visible nonlinearities
                        z_nl_scf[i][j] = -30 #arbitrary large number: no visible nonlinearities
                        Omega_nl_scf[i][j] = 1e-30 #arbitrary small number
                    print "z_nl %f fEDE %f"%(z_nl_scf[i][j],Omega_nl_scf[i][j])
                except CosmoComputationError: # this happens when CLASS fails
                    k_nl_scf[i][j] = 100 #arbitrary large number: bug
                    z_nl_scf[i][j] = -30 #arbitrary large number: bug
                    Omega_nl_scf[i][j] = 1e-30 #arbitrary small number
                f1.write(str(Omega_nl_scf[i][j]) + '\t\t') # write this value of H0 to the file
                f2.write(str(z_nl_scf[i][j]) + '\t\t') # write this value of H0 to the file
                f3.write(str(k_nl_scf[i][j]) + '\t\t') # write this value of H0 to the file
    			# print('fEDE = %e \t ac = %e \t alpha = %.5f \t mu = %.5f\n' %(fEDE[i], ac[j], alpha[i][j], mu[i][j]))
    			# print('fEDE = %e \t mu = %e \t alpha = %.5f \t zc = %.5f\n' %(fEDE[i], mu[j], alpha[i][j], zc[i][j]))

    			# # test that stuff is working by plotting the fluid energy density
    			# bg = cosmo.get_background()
    			# plt.loglog( bg['z'], bg['(.)rho_fld[0]'])
    			# plt.show()

    f1.close()
    f2.close()
    f3.close()
###############################################################################################################################################
# PLOT THINGS
###############################################################################################################################################

end = time()
print('\n\nTime taken for everything but plot = %.2f seconds' %(end - start))

fig = plt.figure(figsize=(15,10))

ax1 = fig.add_subplot(221)
# pcolormesh is fast, simple and accurate
# print ratio_max
print "z nl:",z_nl_scf
cax1 = ax1.pcolormesh( ac, Theta_initial, np.transpose(z_nl_scf), # for some reason, need to transpose H0 to get it to plot correctly
				cmap='coolwarm',vmin=0 # choose colour map
				)
cbar = fig.colorbar(cax1) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
cbar.set_label(r"${\rm Log}_{10}(z_{\rm nl})$", fontsize = 16)
cbar.ax.tick_params(labelsize=16)

ax1.set_ylabel(r"$\Theta_i$", fontsize=16)
ax1.set_xlabel(r"$a_c$", fontsize=16)
# ax1.set_title(r"$z_{\rm nl}$, n=%.1f"%(n_axion), fontsize=16)
ax1.tick_params('x',labelsize=16)
ax1.tick_params('y',labelsize=16)


ax2 = fig.add_subplot(222)
print "Omega_nl_scf:",Omega_nl_scf
cax2 = ax2.pcolormesh( ac, Theta_initial, np.transpose(np.log10(Omega_nl_scf)), # for some reason, need to transpose H0 to get it to plot correctly
				cmap='coolwarm',vmin=-6 # choose colour map
				)
cbar = fig.colorbar(cax2) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
cbar.set_label(r"${\rm Log}_{10}(f_{\rm EDE}(z_{\rm nl}))$", fontsize = 16)
cbar.ax.tick_params(labelsize=16)

ax2.set_ylabel(r"$\Theta_i$", fontsize=16)
ax2.set_xlabel(r"$a_c$", fontsize=16)
ax2.tick_params('x',labelsize=16)
ax2.tick_params('y',labelsize=16)
# ax2.set_title(r"$f_{\rm EDE}(z_{\rm nl})$, n=%.1f"%(n_axion), fontsize=16)

ax3 = fig.add_subplot(223)

print "k nl:",k_nl_scf


cax3 = ax3.pcolormesh( ac, Theta_initial, np.transpose(np.log10(k_nl_scf)), # for some reason, need to transpose H0 to get it to plot correctly
				cmap='coolwarm',vmin=-2,vmax=0# choose colour map
				)
cbar = fig.colorbar(cax3) #, ticks=[H0.min(), 0.5*(H0.min() + H0.max()), H0.max()]) # cax is the imshow object above
cbar.set_label(r"${\rm Log}_{10}(k_{\rm nl})$", fontsize = 16)

cbar.ax.tick_params(labelsize=16)

ax3.set_ylabel(r"$\Theta_i$", fontsize=16)
ax3.set_xlabel(r"$a_c$", fontsize=16)
ax3.tick_params('x',labelsize=16)
ax3.tick_params('y',labelsize=16)
# ax3.set_title(r"$k_{\rm nl}$, n=%.1f"%(n_axion), fontsize=16)

ax4 = fig.add_subplot(224)
print "znl:", z_nl_scf[np.where(ac==-4.0)][0]
ax4.plot(Theta_initial/np.pi,z_nl_scf[np.where(ac==-3.5)][0],'b')
ax4.set_xlabel(r"$\Theta_i/\pi$", fontsize=16)
# ax4.set_yticks([0.5,1,5])
ax4.set_ylabel(r"$z_{\rm nl}$", fontsize=16, color='b')
ax4.tick_params('x', labelsize=16)
ax4.tick_params('y', labelcolor='b',labelsize=16)
# ax4.set_yscale('log')

#
ax4_bis = ax4.twinx()
ax4_bis.plot(Theta_initial/np.pi,1000*Omega_nl_scf[np.where(ac==-4)][0],'r')
ax4_bis.set_ylabel(r"$1000\times f_{\rm EDE}(z_{\rm nl})$", fontsize=16, color='r')
ax4_bis.tick_params('y', labelcolor='r',labelsize=16)
# ax4_bis.set_yscale('log')
ax4.text(0.0,max(z_nl_scf[np.where(ac==-3.5)][0])*0.9,r"$a_c = 10^{-4}$, n=%.1f"%(n_axion), fontsize=16)
# plt.savefig(output_file + '_zc_alpha.png',bbox_inches='tight')
fig.tight_layout()
plt.show()
