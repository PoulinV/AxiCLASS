# coding: utf-8

# In[ ]:

# import necessary modules
#get_ipython().magic(u'matplotlib notebook')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# In[ ]:

# esthetic definitions for the plots
font = {'size'   : 16, 'family':'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


# In[ ]:

#############################################
#
# User settings controlling the figure aspect
#
z_max_pk = 46000       # highest redshift involved
k_per_decade = 100     # number of k values, controls final resolution
k_min_tau0 = 40.       # this value controls the minimum k value in the figure (it is k_min * tau0)
P_k_max_inv_Mpc =1.0   # this value is directly the maximum k value in the figure in Mpc
tau_num_early = 2000   # number of conformal time values before recombination, controls final resolution
tau_num_late = 2000     # number of conformal time values after recombination, controls final resolution
tau_ini = 10.          # first value of conformal time in Mpc
tau_label_Hubble = 20. # value of time at which we want to place the label on Hubble crossing
tau_label_ks = 40.     # value of time at which we want to place the label on sound horizon crossing
tau_label_kd = 230.    # value of time at which we want to place the label on damping scale crossing
#
# Cosmological parameters and other CLASS parameters
#
K_star = 0.05
log10_axion_ac = -3.7
n_axion = 3
wn = (n_axion-1)/(n_axion+1)
A_s = 2.215e-9
n_s = 0.9619
common_settings = {# which output? transfer functions only
                   'output':'mTk',
                   # LambdaCDM parameters
                   'h':0.67556,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':A_s,
                   'n_s':n_s,
                   'tau_reio':0.0925,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe':0.246,
                   # other output and precision parameters
                   'z_max_pk':z_max_pk,
                   'recfast_z_initial':z_max_pk,
                   'start_small_k_at_tau_c_over_tau_h':5e-08,
                   'start_large_k_at_tau_h_over_tau_k':5e-08,
                   #'k_step_sub':'0.01',
                   'k_per_decade_for_pk':k_per_decade,
                   'k_per_decade_for_bao':k_per_decade,
                   'k_min_tau0':k_min_tau0, # this value controls the minimum k value in the figure
                   # 'perturb_sampling_stepsize':'0.05',
                   'P_k_max_1/Mpc':P_k_max_inv_Mpc,
                   'compute damping scale':'yes', # needed to output and plot Silk damping scale
                   'gauge':'synchronous',
                   'scf_potential': 'axion',
                   'n_axion': n_axion,
                   'log10_axion_ac': log10_axion_ac, # Must input log10(axion_ac)
                   # log10_fraction_axion_ac': -1.922767 # Must input log10(fraction_axion_ac)
                   'log10_fraction_axion_ac': -1, # Must input log10(fraction_axion_ac)
                   # m_axion': 1.811412e+06
                   # f_axion': 1
                   'scf_parameters':'3.1,0.0', #phi_i,phi_dot_i //dummy: phi_i will be updated.
                   'adptative_stepsize': 100,
                   'scf_tuning_index': 0,
                   'do_shooting': 'yes',
                   'do_shooting_scf': 'yes',
                   # back_integration_stepsize':1'e-4
                   'use_big_theta_scf': 'yes',
                   'scf_has_perturbations': 'yes',
                   'attractor_ic_scf': 'no'}

###############
#
# call CLASS
#
###############
M = Class()
M.set(common_settings)
M.compute()
#
# define conformal time sampling array
#
times = M.get_current_derived_parameters(['tau_rec','conformal_age'])
tau_rec=times['tau_rec']
tau_0 = times['conformal_age']
tau1 = np.logspace(math.log10(tau_ini),math.log10(tau_rec),tau_num_early)
tau2 = np.logspace(math.log10(tau_rec),math.log10(tau_0),tau_num_late)[1:]
tau2[-1] *= 0.999 # this tiny shift avoids interpolation errors
tau = np.concatenate((tau1,tau2))
tau_num = len(tau)
#
# use table of background and thermodynamics quantitites to define some functions
# returning some characteristic scales
# (of Hubble crossing, sound horizon crossing, etc.) at different time
#
background = M.get_background() # load background table
#print background.viewkeys()
thermodynamics = M.get_thermodynamics() # load thermodynamics table
#print thermodynamics.viewkeys()
#
background_tau = background['conf. time [Mpc]'] # read conformal times in background table
background_z = background['z'] # read redshift
background_aH = 2.*math.pi*background['H [1/Mpc]']/(1.+background['z'])/M.h() # read 2pi * aH in [h/Mpc]
background_ks = 2.*math.pi/background['comov.snd.hrz.']/M.h() # read 2pi/(comoving sound horizon) in [h/Mpc]
background_rho_m_over_r =    (background['(.)rho_b']+background['(.)rho_cdm'])    /(background['(.)rho_g']+background['(.)rho_ur']) # read rho_r / rho_m (to find time of equality)
background_rho_l_over_m =    background['(.)rho_lambda']    /(background['(.)rho_b']+background['(.)rho_cdm']) # read rho_m / rho_lambda (to find time of equality)
background_rho_scf =    background['(.)rho_scf']    # read rho_m / rho_lambda (to find time of equality)
background_rho_g =    background['(.)rho_g']    # read rho_m / rho_lambda (to find time of equality)
background_rho_tot =  background['(.)rho_lambda'] +  background['(.)rho_scf']+background['(.)rho_b']+background['(.)rho_cdm']+ background['(.)rho_g']+background['(.)rho_ur']   # read rho_m / rho_lambda (to find time of equality)
thermodynamics_tau = thermodynamics['conf. time [Mpc]'] # read confromal times in thermodynamics table
thermodynamics_kd = 2.*math.pi/thermodynamics['r_d']/M.h() # read 2pi(comoving diffusion scale) in [h/Mpc]
#
# define a bunch of interpolation functions based on previous quantities
#
background_z_at_tau = interp1d(background_tau,background_z)
background_tau_at_z = interp1d(background_z,background_tau)
background_rho_scf_at_tau = interp1d(background_tau,background_rho_scf)
background_rho_tot_at_tau = interp1d(background_tau,background_rho_tot)
background_rho_g_at_tau = interp1d(background_tau,background_rho_g)
background_aH_at_tau = interp1d(background_tau,background_aH)
background_ks_at_tau = interp1d(background_tau,background_ks)
background_tau_at_mr = interp1d(background_rho_m_over_r,background_tau)
background_tau_at_lm = interp1d(background_rho_l_over_m,background_tau)
thermodynamics_kd_at_tau = interp1d(thermodynamics_tau, thermodynamics_kd)
#
# infer arrays of characteristic quantitites calculated at values of conformal time in tau array
#
aH = background_aH_at_tau(tau)
ks = background_ks_at_tau(tau)
kd = thermodynamics_kd_at_tau(tau)
#
# infer times of R/M and M/Lambda equalities
#
tau_eq = background_tau_at_mr(1.)
tau_lambda = background_tau_at_lm(1.)
#
# check and inform user whether intiial arbitrary choice of z_max_pk was OK
max_z_needed = background_z_at_tau(tau[0])
if max_z_needed > z_max_pk:
    print 'you must increase the value of z_max_pk to at least ',max_z_needed
    () + 1  # this strange line is just a trick to stop the script execution there
else:
    print 'in a next run with the same values of tau, you may decrease z_max_pk from ',z_max_pk,' to ',max_z_needed
#
# get transfer functions at each time and build arrays delta_phi_scf(tau,k) and phi(tau,k)
#
rho_scf = []
rho_tot = []
rho_g = []
one_time = M.get_transfer(10)
print one_time.viewkeys()
for i in range(tau_num):
    one_time = M.get_transfer(background_z_at_tau(tau[i])) # transfer functions at each time tau
    # rho_scf.append(background_rho_scf_at_tau(tau[i]))
    # rho_g.append(background_rho_g_at_tau(tau[i]))
    # rho_tot.append(background_rho_tot_at_tau(tau[i]))
    if i ==0:   # if this is the first time in the loop: create the arrays (k, delta_phi_scf, phi)
        k = one_time['k (h/Mpc)']
        k_num = len(k)
        delta_rho_phi_scf = np.zeros((tau_num,k_num))
        delta_phi_scf = np.zeros((tau_num,k_num))
        delta_tot = np.zeros((tau_num,k_num))
        phi = np.zeros((tau_num,k_num))
        delta_phi_prime_scf = np.zeros((tau_num,k_num))
        delta_cdm = np.zeros((tau_num,k_num))
        delta_g = np.zeros((tau_num,k_num))

    delta_rho_phi_scf[i,:] = np.abs(one_time['d_scf'][:])*background_rho_scf_at_tau(tau[i])
    delta_g[i,:] = np.abs(one_time['d_g'][:])*background_rho_g_at_tau(tau[i])
    delta_tot[i,:] = np.abs(one_time['d_tot'][:])*background_rho_tot_at_tau(tau[i])
    phi[i,:] = np.abs(one_time['phi'][:])
    delta_phi_prime_scf[i,:] = np.log10(np.abs(one_time['delta_phi_prime_scf'][:]))
    delta_phi_scf[i,:] = np.log10(np.abs(one_time['delta_phi_scf'][:]))
    delta_cdm[i,:] = np.abs(one_time['d_cdm'][:])
#
# find the global extra of delta_phi_scf(tau,k) and phi(tau,k), used to define color code later
#
delta_phi_scf_amp = max(delta_phi_scf.max(),-delta_phi_scf.min())
delta_phi_prime_scf_amp = max(delta_phi_prime_scf.max(),-delta_phi_prime_scf.min())
#
# reshaping of (k,tau) necessary to call the function 'pcolormesh'
#
K,T = np.meshgrid(k,tau)
#
# inform user of the size of the grids (related to the figure resolution)
#
print 'grid size:',len(k),len(tau),delta_phi_scf.shape
#
#################
#
# start plotting
#
#################
#
fig = plt.figure(figsize=(18,8))
#
# plot delta_phi_scf(k,tau)
#
background_phi_scf = fig.add_subplot(131)
background_phi_scf.plot(background_rho_scf/background_rho_tot,background_tau,linewidth=2)
background_phi_scf.set_ylim(tau[0],tau[-1])
background_phi_scf.set_yscale('log')
background_phi_scf.set_xlabel(r'$f_{\rm EDE}(\tau)$')
background_phi_scf.set_ylabel(r'$\tau   \,\,\, \mathrm{[Mpc]}$')
background_phi_scf.invert_yaxis()

ax_delta_phi_scf = fig.add_subplot(132)
print '> Plotting delta_phi_scf_0'
# fig_delta_phi_scf = ax_delta_phi_scf.pcolormesh(K,T,np.log10(delta_phi_scf*delta_phi_scf*K*K*K/2/np.pi/np.pi*A_s*(K/K_star)**(n_s-1)),cmap='coolwarm',vmin=np.log10(delta_phi_scf.min()*delta_phi_scf.min()*K.min()*K.min()*K.min()/2/np.pi*A_s*(K.min()/K_star)**(n_s-1)), vmax=np.log10(delta_phi_scf.max()*delta_phi_scf.max()*K.max()*K.max()*K.max()/2/np.pi*A_s*(K.max()/K_star)**(n_s-1))) #,shading='gouraud')
# fig_delta_phi_scf = ax_delta_phi_scf.pcolormesh(K,T,np.log10(delta_phi_scf*delta_phi_scf*A_s*(K/K_star)**(n_s-1)),cmap='coolwarm',vmin=np.log10(delta_phi_scf.min()*delta_phi_scf.min()*A_s*(K.min()/K_star)**(n_s-1)), vmax=np.log10(delta_phi_scf.max()*delta_phi_scf.max()*A_s*(K.max()/K_star)**(n_s-1))) #,shading='gouraud')
##delta_phi power spectrum
phi_enveloppe = 3.1*10**M.log10_f_axion()*(2/(10**log10_axion_ac*background_z_at_tau(T)))**(-3.*(1+wn)/2/n_axion)
background_phi_scf_v2 = background_phi_scf.twiny()
background_phi_scf_v2.set_ylim(tau[0],tau[-1])
background_phi_scf_v2.set_xlim(phi_enveloppe[0],2*phi_enveloppe[-1])
background_phi_scf_v2.plot(phi_enveloppe,tau,'r-',linewidth=2)
background_phi_scf_v2.set_xlabel(r'$\phi_{\rm env}$')
background_phi_scf_v2.invert_yaxis()

fig_delta_phi_scf = ax_delta_phi_scf.pcolormesh(K,T,np.log10(np.sqrt(delta_phi_scf*delta_phi_scf*A_s*(K/K_star)**(n_s-1))/(phi_enveloppe)),cmap='coolwarm') #,shading='gouraud')
##fractional contribution to delta_tot_squared
# fig_delta_phi_scf = ax_delta_phi_scf.pcolormesh(K,T,np.log10(delta_phi_scf/delta_tot),cmap='coolwarm',vmin=-14,vmax=0) #,shading='gouraud')

print '> Done'
#
# plot lines (characteristic times and scales)
#
# ax_delta_phi_scf.axhline(y=tau_rec,color='k',linestyle='-')
ax_delta_phi_scf.axhline(y=background_tau_at_z(1./(10**log10_axion_ac)),color='white',linestyle='-')
ax_delta_phi_scf.axhline(y=tau_eq,color='k',linestyle='-')
ax_delta_phi_scf.axhline(y=tau_lambda,color='k',linestyle='-')
ax_delta_phi_scf.plot(aH,tau,'r-',linewidth=2)
# ax_delta_phi_scf.plot(ks,tau,color='#FFFF33',linestyle='-',linewidth=2)
# ax_delta_phi_scf.plot(kd,tau,'b-',linewidth=2)
#
# dealing with labels
#
# ax_delta_phi_scf.set_title(r'${\rm Log}_{10}(A_s\left(\frac{k}{k_*}\right)^{n_s-1}\delta_\phi)$')
# ax_delta_phi_scf.set_title(r'${\rm Log}_{10}(\delta\rho_\phi/\delta\rho_{\rm tot})$')
ax_delta_phi_scf.set_title(r'${\rm Log}_{10}(\Delta_{\delta\phi}/\phi_{\rm env})$')
# ax_delta_phi_scf.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$')
ax_delta_phi_scf.text(1.5*k[0],0.9*background_tau_at_z(1./(10**log10_axion_ac)),r'$z_c$',color='white')
ax_delta_phi_scf.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$')
ax_delta_phi_scf.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$')
ax_delta_phi_scf.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
                  xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
                  xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
# ax_delta_phi_scf.annotate(r'$\mathrm{sound} \,\, \mathrm{horizon} \,\, \mathrm{cross.}$',
#                   xy=(background_ks_at_tau(tau_label_ks),tau_label_ks),
#                   xytext=(0.07*background_aH_at_tau(tau_label_ks),0.8*tau_label_ks),
#                   arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
# ax_delta_phi_scf.annotate(r'$\mathrm{damping} \,\, \mathrm{scale} \,\, \mathrm{cross.}$',
#                   xy=(thermodynamics_kd_at_tau(tau_label_kd),tau_label_kd),
#                   xytext=(0.2*thermodynamics_kd_at_tau(tau_label_kd),2.0*tau_label_kd),
#                   arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))

# dealing with axes
#
ax_delta_phi_scf.set_xlim(k[0],k[-1])
ax_delta_phi_scf.set_xscale('log')
ax_delta_phi_scf.set_yscale('log')
ax_delta_phi_scf.set_xlabel(r'$k  \,\,\, \mathrm{[h/Mpc]}$')
ax_delta_phi_scf.set_ylabel(r'$\tau   \,\,\, \mathrm{[Mpc]}$')
ax_delta_phi_scf.invert_yaxis()
#
# color legend
#
# axins = inset_axes(fig_delta_phi_scf,
#                    width="5%",  # width = 5% of parent_bbox width
#                    height="50%",  # height : 50%
#                    loc='lower right',
#                    bbox_to_anchor=(1.05, 0., 1, 1),
#                    # bbox_transform=fig_delta_phi_scf.transAxes,
#                    borderpad=0,
#                    )
fig.colorbar(fig_delta_phi_scf)
#
# plot phi(k,tau)
#
ax_delta_phi_scf_prime = fig.add_subplot(133)
ax_delta_phi_scf_prime.set_xlim(k[0],k[-1])
#ax_delta_phi_scf_prime.pcolor(K,T,phi,cmap='coolwarm')
print '> Plotting phi'
# fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,np.log10(delta_cdm*delta_cdm*A_s*(K/K_star)**(n_s-1)),cmap='coolwarm',vmin=np.log10(delta_cdm.min()*delta_cdm.min()*A_s*(K.min()/K_star)**(n_s-1)), vmax=np.log10(delta_cdm.max()*delta_cdm.max()*A_s*(K.max()/K_star)**(n_s-1))) #,shading='gouraud')
# fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,np.log10(delta_cdm*delta_cdm*A_s*(K/K_star)**(n_s-1)),cmap='coolwarm') #,shading='gouraud')
# fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,np.log10(delta_cdm*delta_cdm*A_s*(K/K_star)**(n_s-1)),cmap='coolwarm') #,shading='gouraud')
# fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,np.log10(delta_g/delta_tot),cmap='coolwarm',vmin=-14,vmax=0) #,shading='gouraud')
fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,np.log10(delta_rho_phi_scf/delta_tot),cmap='coolwarm',vmin=-14,vmax=0) #,shading='gouraud')
# fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,np.log10(delta_g/delta_tot),cmap='coolwarm',vmin=-14,vmax=0) #,shading='gouraud')
# fig_phi = ax_delta_phi_scf_prime.pcolormesh(K,T,delta_phi_prime_scf,cmap='coolwarm',vmin=delta_phi_prime_scf.min(), vmax=delta_phi_prime_scf.max())
print '> Done'
#
# plot lines (characteristic times and scales)
#
# ax_delta_phi_scf_prime.axhline(y=tau_rec,color='k',linestyle='-')
ax_delta_phi_scf_prime.axhline(y=tau_eq,color='k',linestyle='-')
ax_delta_phi_scf_prime.axhline(y=tau_lambda,color='k',linestyle='-')
ax_delta_phi_scf_prime.plot(aH,tau,'r-',linewidth=2)
# ax_delta_phi_scf_prime.plot(ks,tau,color='#FFFF33',linestyle='-',linewidth=2)
#
# dealing with labels
#
# ax_delta_phi_scf_prime.set_title(r'${\rm Log}_{10}(\delta\rho_g/\delta\rho_{\rm tot})$')
ax_delta_phi_scf_prime.set_title(r'${\rm Log}_{10}(\delta\rho_\phi/\delta\rho_{\rm tot})$')

# ax_delta_phi_scf_prime.set_title(r'${\rm Log}_{10}(\delta\dot{\phi})$')
# ax_delta_phi_scf_prime.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$')
ax_delta_phi_scf_prime.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$')
ax_delta_phi_scf_prime.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$')
ax_delta_phi_scf_prime.axhline(y=background_tau_at_z(1./(10**log10_axion_ac)),color='white',linestyle='-')
ax_delta_phi_scf_prime.text(1.5*k[0],0.9*background_tau_at_z(1./(10**log10_axion_ac)),r'$z_c$',color='white')

ax_delta_phi_scf_prime.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
                  xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
                  xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
                  arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
# ax_delta_phi_scf_prime.annotate(r'$\mathrm{sound} \,\, \mathrm{horizon} \,\, \mathrm{cross.}$',
#                   xy=(background_ks_at_tau(tau_label_ks),tau_label_ks),
#                   xytext=(0.07*background_aH_at_tau(tau_label_ks),0.8*tau_label_ks),
#                   arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
#
# dealing with axes

ax_delta_phi_scf_prime.set_xscale('log')
ax_delta_phi_scf_prime.set_yscale('log')
ax_delta_phi_scf_prime.set_xlabel(r'$k \,\,\, \mathrm{[h/Mpc]}$')
ax_delta_phi_scf_prime.set_ylabel(r'$\tau \,\,\, \mathrm{[Mpc]}$')
ax_delta_phi_scf_prime.invert_yaxis()
#
# color legend
#
fig.colorbar(fig_phi)
#
# produce the PDF
#
#plt.show()
plt.savefig('amin_criterion.png',dpi=300)
