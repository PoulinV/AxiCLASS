
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
from matplotlib import rc
from scipy.special import gamma, factorial
rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)

# In[ ]:

# esthetic definitions for the plots
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


# In[ ]:

#############################################
#
# User settings controlling the figure aspect
#
z_max_pk = 46000       # highest redshift involved
k_per_decade = 400     # number of k values, controls final resolution
k_min_tau0 = 40.       # this value controls the minimum k value in the figure (it is k_min * tau0)
P_k_max_inv_Mpc =1.0   # this value is directly the maximum k value in the figure in Mpc
tau_num_early = 2000   # number of conformal time values before recombination, controls final resolution
tau_num_late = 200     # number of conformal time values after recombination, controls final resolution
tau_ini = 10.          # first value of conformal time in Mpc
tau_label_Hubble = 70. # value of time at which we want to place the label on Hubble crossing
tau_label_ks = 40.     # value of time at which we want to place the label on sound horizon crossing
tau_label_kd = 230.    # value of time at which we want to place the label on damping scale crossing
#
# Cosmological parameters and other CLASS parameters
#
Contours_at = (0.51,0.6,0.9)
Theta_initial_fld=0.1
n_axion = 3
log10ac = -3.5
common_settings = {'scf_potential' : 'axion',
                    'output':'mTk',
                    'n_axion' : n_axion,
                    'scf_parameters' : '%.1f,0'%(Theta_initial_fld),
                    'scf_tuning_index':0,
                    'log10_fraction_axion_ac':-0.8769,
                    'log10_axion_ac':log10ac,
                    'omega_cdm':1.320351e-01,
                    'omega_b': 2.263127e-02,
                    # '100*theta_s':1.04262e+00,
                    # '100*theta_s':1.04262e+00,
                    'H0':7.238797e+01,
                    'tau_reio':6.402188e-02,
                    'A_s':2.178373e-09,
                    'n_s':9.917292e-01,
                    'scf_evolve_like_axionCAMB':'no',
                    'do_shooting':'yes',
                    'do_shooting_scf':'yes',
                    'use_big_theta_scf':'yes',
                    'scf_has_perturbations':'yes',
                    'attractor_ic_scf':'no',
                    'adptative_stepsize':100,
                    'scf_evolve_as_fluid':'no',
                    'z_max_pk':z_max_pk,
                    'recfast_z_initial':z_max_pk,
                    #'k_step_sub':'0.01',
                    'k_per_decade_for_pk':k_per_decade,
                    'k_per_decade_for_bao':k_per_decade,
                    'k_min_tau0':k_min_tau0, # this value controls the minimum k value in the figure
                    'perturb_sampling_stepsize':'0.05',
                    'P_k_max_1/Mpc':P_k_max_inv_Mpc,
                    'compute damping scale':'yes', # needed to output and plot Silk damping scale
                    'gauge':'synchronous'}

# common_settings['h'] = 0.67
# common_settings['back_integration_stepsize'] = 5e-4

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
# background_aH = 2.*math.pi*background['H [1/Mpc]']/(1.+background['z'])/M.h() # read 2pi * aH in [h/Mpc]
background_aH = background['H [1/Mpc]']/(1.+background['z'])/M.h() # read 2pi * aH in [h/Mpc]
background_ks = 2.*math.pi/background['comov.snd.hrz.']/M.h() # read 2pi/(comoving sound horizon) in [h/Mpc]
background_rho_m_over_r =    (background['(.)rho_b']+background['(.)rho_cdm'])    /(background['(.)rho_g']+background['(.)rho_ur']) # read rho_r / rho_m (to find time of equality)
background_rho_l_over_m =    background['(.)rho_lambda']    /(background['(.)rho_b']+background['(.)rho_cdm']) # read rho_m / rho_lambda (to find time of equality)
thermodynamics_tau = thermodynamics['conf. time [Mpc]'] # read confromal times in thermodynamics table
# thermodynamics_kd = 2.*math.pi/thermodynamics['r_d']/M.h() # read 2pi(comoving diffusion scale) in [h/Mpc]
#
# define a bunch of interpolation functions based on previous quantities
#
background_z_at_tau = interp1d(background_tau,background_z)
background_tau_at_z = interp1d(background_z,background_tau)
background_aH_at_tau = interp1d(background_tau,background_aH)
background_ks_at_tau = interp1d(background_tau,background_ks)
background_tau_at_mr = interp1d(background_rho_m_over_r,background_tau)
background_tau_at_lm = interp1d(background_rho_l_over_m,background_tau)
#
# infer arrays of characteristic quantitites calculated at values of conformal time in tau array
#
aH = background_aH_at_tau(tau)
ks = background_ks_at_tau(tau)
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
# get transfer functions at each time and build arrays Theta0(tau,k) and phi(tau,k)
#
i_at_ac =0
for i in range(tau_num):
    one_time = M.get_transfer(background_z_at_tau(tau[i])) # transfer functions at each time tau
    if i ==0:   # if this is the first time in the loop: create the arrays (k, Theta0, phi)
        k = one_time['k (h/Mpc)']
        k_num = len(k)
        Theta0 = np.zeros((tau_num,k_num))
        cs2_fld = np.zeros((tau_num,k_num))
        phi = np.zeros((tau_num,k_num))
    # Theta0[i,:] = 0.25*one_time['d_g'][:]
    a = 1/(1+background_z_at_tau(tau[i]))
    if a > 10**log10ac and i_at_ac == 0:
        i_at_ac=i
    a_osc=10**log10ac/1.7
    w_n = (n_axion-1)/(n_axion+1)
    theta_osc = Theta_initial_fld *(a/a_osc)**(-3*(1+w_n)/(2*n_axion))
    c = 3e5
    m = 10**M.log10_m_axion()*M.h()*100/c
    # print "m:",m
    omega = m*np.sqrt(np.pi)*gamma((1.0+n_axion)/(2.0*n_axion))/gamma(1.0+1.0/(2*n_axion))*2**(-(1.0+n_axion)/2.0)*theta_osc**(n_axion-1.0)
    # print "omega:",omega
    cs2_fld[i,:] = (2*a*a*(n_axion-1)*omega*omega+k*k)/(2*a*a*(n_axion+1)*omega*omega+k*k)
    phi[i,:] = one_time['phi'][:]
# Omega_scf = bg['(.)Omega_scf']
k_of_cs2_at_ac = interp1d(cs2_fld[i_at_ac,:],k)
print cs2_fld[i_at_ac,:],i_at_ac
print "k_of_cs2_at_ac=",k_of_cs2_at_ac(0.9)
H = background['H [1/Mpc]']
Da = background['ang.diam.dist.']
z = background['z']
# tau = background['conf. time [Mpc]']
# H_of_z = interp1d(z,H)
zc = 1/(10**log10ac)-1
tau_ac = background_tau_at_z(zc)
k_ac = background_aH_at_tau(tau_ac)
# print "k_ac=",k_ac
# for j in range(len(Omega_scf)):
#     if Omega_scf[j] > 0.01:
        # k = H[j]/(1+z[j])
#
# find the global extra of Theta0(tau,k) and phi(tau,k), used to define color code later
#
Theta_amp = max(Theta0.max(),-Theta0.min())
phi_amp = max(phi.max(),-phi.min())
#
# reshaping of (k,tau) necessary to call the function 'pcolormesh'
#
K,T = np.meshgrid(k,tau)
#
# inform user of the size of the grids (related to the figure resolution)
#
print 'grid size:',len(k),len(tau),Theta0.shape
#
#################
#
# start plotting
#
#################
#
if v1 == True:
    fig = plt.figure(figsize=(12,8))
    #
    # plot Theta0(k,tau)
    #
    ax_Theta = fig.add_subplot(111)
    print '> Plotting Theta_0'
    # fig_Theta = ax_Theta.pcolormesh(K,T,cs2_fld,cmap='coolwarm',vmin=cs2_fld.min(), vmax=cs2_fld.max()) #,shading='gouraud')
    fig_Theta = ax_Theta.pcolormesh(K,T,cs2_fld,cmap='coolwarm',vmin=cs2_fld.min(), vmax=cs2_fld.max()) #,shading='gouraud')
    print '> Done'
    #
    # plot lines (characteristic times and scales)
    #
    # ax_Theta.axhline(y=tau_rec,color='k',linestyle='-')
    # ax_Theta.axhline(y=tau_eq,color='k',linestyle='-')
    # ax_Theta.axhline(y=tau_lambda,color='k',linestyle='-')
    ax_Theta.axhline(y=tau_ac,color='white',linestyle='-',linewidth=2,zorder=1)
    # ax_Theta.axvline(x=k_ac,color='k',linestyle='-')
    ax_Theta.plot(aH,tau,'--',color='white',linewidth=2,zorder=1)
    # ax_Theta.plot(ks,tau,color='#FFFF33',linestyle='-',linewidth=2)
    # ax_Theta.plot(kd,tau,'b-',linewidth=2)
    #
    # dealing with labels
    #
    k_range = np.linspace(k_ac,k_of_cs2_at_ac(0.9),100)
    print "k_range",k_range
    ax_Theta.set_title(r'$n=%d,\theta_i=%.1f,{\rm log}_{10}(a_c)=%.1f$'%(n_axion,Theta_initial_fld,log10ac), fontsize = 22)
    # ax_Theta.text(1.5*k[0],0.9*tau_rec,r'$\mathrm{rec.}$', fontsize = 22)
    ax_Theta.text(0.6,0.9*tau_ac,r'$\tau(a_c)$', fontsize = 22,color='white')
    # ax_Theta.text(1.5*k[0],0.9*tau_eq,r'$\mathrm{R/M} \,\, \mathrm{eq.}$', fontsize = 22)
    # ax_Theta.text(1.5*k[0],0.9*tau_lambda,r'$\mathrm{M/L} \,\, \mathrm{eq.}$', fontsize = 22)
    # ax_Theta.text(1.1*k_ac,1.1*tau[0],r'$k=a_cH(a_c)$', fontsize = 22)
    # ax_Theta.annotate(r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',
    #                   xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
    #                   xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
    #                   arrowprops=dict(color='white', shrink=0.05, width=1, headlength=2, headwidth=2), fontsize = 22,color='white')
    # ax_Theta.text(0.3,17,r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',fontsize = 22,color='white',rotation=35)
    # ax_Theta.text(0.08,17,r'$\mathrm{Hubble}$',fontsize = 22,color='white',rotation=35)
    ax_Theta.text(0.04,17,r'$\mathrm{Hubble}$',fontsize = 22,color='white',rotation=35)
    # ax_Theta.text(0.3,17,r'$\mathrm{Hubble} \,\, \mathrm{cross.}$',fontsize = 22,color='white',rotation=35)
    #                   xy=(background_aH_at_tau(tau_label_Hubble),tau_label_Hubble),
    #                   xytext=(0.1*background_aH_at_tau(tau_label_Hubble),0.8*tau_label_Hubble),
    #                   arrowprops=dict(color='white', shrink=0.05, width=1, headlength=2, headwidth=2), fontsize = 22,color='white')
    # ax_Theta.annotate(r'$\mathrm{sound} \,\, \mathrm{horizon} \,\, \mathrm{cross.}$',
    #                   xy=(background_ks_at_tau(tau_label_ks),tau_label_ks),
    #                   xytext=(0.07*background_aH_at_tau(tau_label_ks),0.8*tau_label_ks),
    #                   arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
    # ax_Theta.annotate(r'$\mathrm{damping} \,\, \mathrm{scale} \,\, \mathrm{cross.}$',
    #                   xy=(thermodynamics_kd_at_tau(tau_label_kd),tau_label_kd),
    #                   xytext=(0.2*thermodynamics_kd_at_tau(tau_label_kd),2.0*tau_label_kd),
    #                   arrowprops=dict(facecolor='black', shrink=0.05, width=1, headlength=5, headwidth=5))
    #
    # dealing with axes
    #
    ax_Theta.set_xlim(k[0],k[-1])
    ax_Theta.set_xscale('log')
    ax_Theta.set_yscale('log')
    ax_Theta.set_xlabel(r'$k  \,\,\, \mathrm{[h/Mpc]}$', fontsize = 22)
    ax_Theta.set_ylabel(r'$\tau   \,\,\, \mathrm{[Mpc]}$', fontsize = 22)
    ax_Theta.invert_yaxis()
    #
    # color legend
    #
    cbar = fig.colorbar(fig_Theta)
    cbar.set_label(r'$c_s^2$', fontsize = 22)
    ax_Theta.tick_params(axis='x',labelsize=23)
    ax_Theta.tick_params(axis='y',labelsize=23)
    ticklabs = cbar.ax.get_yticklabels()

    cbar.ax.set_yticklabels(ticklabs, fontsize=20)



    contour = ax_Theta.contour(K, T, cs2_fld, Contours_at,
            extent = (min(K[0]), max(K[0]), min(T[0]), max(T[0]) ), # place colour map at correct a_c and Om_fld values
            origin = 'lower',
            corner_mask = True,
            colors=('w','w','w','w'),linewidths=('2','2','2','2'))
    # manual_locations = [(0.03, 300),(0.03, 40)]
    # manual_locations = [(0.03, 300),(0.01, 4000)]
    # #2p8: manual_locations = [(0.006, 200),(0.01, 400),(0.01, 1000)]
    #0p1:
    manual_locations = [(0.01, 400),(0.01, 30),(0.05, 50)]
    ax_Theta.clabel(contour, fmt='%.2f', colors='w', fontsize=22,ls='--',lw='2',manual=manual_locations,zorder=1)
    # ax_Theta.clabel(contour, fmt='%.1f', colors='w', fontsize=22,ls='--',lw='2')
    ax_Theta.fill_between(k_range, tau_ac-0.11*tau_ac, tau_ac+0.15*tau_ac, facecolor='blue', alpha=0.5,zorder=2)

    plt.savefig('cs2_n3_0p1_vAxiCLASSv2.png',dpi=300)
# plt.savefig('cs2_n3_2p8.png',dpi=300)
# plt.savefig('cs2_n3_2p8.pdf')
else:
    
