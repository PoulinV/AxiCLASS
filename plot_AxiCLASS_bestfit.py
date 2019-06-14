
# coding: utf-8

# In[ ]:

# import necessary modules
#get_ipython().magic(u'matplotlib inline')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
import math
from matplotlib import rc
import matplotlib.patches as patches

from scipy.interpolate import interp1d
from matplotlib.ticker import FixedLocator
from math import floor
from mpl_toolkits.axes_grid1 import make_axes_locatable

# In[ ]:

# esthetic definitions for the plots

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
# matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]


l_TT_low,Dl_TT_low,err_TT_low= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TT-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
l_TE_low,Dl_TE_low,err_TE_low= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TE-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
l_TT_high,Dl_TT_high,err_TT_high= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TT-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
l_TE_high,Dl_TE_high,err_TE_high= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-TE-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
l_EE_low,Dl_EE_low,err_EE_low= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-EE-loL-full_R2.02.txt",unpack=True,usecols=(0,1,2))
l_EE_high,Dl_EE_high,err_EE_high= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/COM_PowerSpect_CMB-EE-hiL-binned_R2.02.txt",unpack=True,usecols=(0,3,4))
lmin_phiphi,lmax_phiphi,cl_phiphi,err_phiphi= np.loadtxt("/Users/poulin/Dropbox/Labo/ProgrammeDarkAges/error_Planck/agressive_lensing.csv",unpack=True)

# In[ ]:

############################################
#
# Varying parameter (others fixed to default)
#
var_name1 = 'n_pheno_axion'
# var_name2 = 'Omega_fld_ac'
var_name2 = 'Omega_many_fld'
# var1_min = 1
# var1_max = 3
# var2_min = 1e14
n  = [2,3]
# n  = [3,3]
# n  = [1,1,1,1,1,1]
# n = [2,2]
# Omega_axion_ac = [1e14,5e14,1e15]
# Omega_axion_ac = [50,50,50]
# Omega_axion_ac = [3,4.5,7]
# Omega_many_fld = [3.1e-3,1.25e-5,0.4e-6]
# Omega_many_fld = [6.167e-3,11.95e-4,4.3e-4]
Omega_many_fld = [4.289e-06,15.8e-9,1,1]######95%CL ac=1em5
# Omega_many_fld = [0.02,4.289e-06,15.8e-9]######95%CL ac=1em5
# Omega_many_fld = [4.289e-06,4.289e-06]######95%CL ac=1em5
# Omega_many_fld = [3.1e-3,1.25e-5,0.3993e-06]######95%CL ac=1em3
# Omega_many_fld = [6.15e-3,1.2e-3,4.3e-4]#######95%CL ac=10
# Omega_many_fld = [6.15e-3,3.1e-3,0.02]
# Omega_many_fld = [0.266,0.266,0.266,0.266,0.266,0.266,0.266]
# Omega_axion_ac = [2.278e+13,4.278e+14,1.8e+15]
a_c = [1e-5,1e-5,1,1]
# color = ['r','b--','g-.','o','m--','k','c']
color = ['crimson','dodgerblue', 'darkorchid','orange']
# color = ['r--','b-.','o','m--','k','c']
var_num = 2
# var_legend = r'$n_{\rm a}$'
# legend_name= [r'$a_c =0.1$',r'$a_c =10^{-3}$',r'$a_c =10^{-5}$',r'$a_c =10^{-5}$',r'$a_c =10^{-5}$',r'$a_c =10^{-5}$']
# legend_name= [r'$a_c =0.1$',r'$a_c =10^{-3}$',r'$a_c =10^{-5}$',r'$a_c =10^{-5}$',r'$a_c =10^{-5}$',r'$a_c =10^{-5}$']
# legend_name= [r'$N_{\rm eff}$',r'$n=2$',r'$n=3$',r'$n=\infty$']
# legend_name= [r'$n=2$',r'$n=3$',r'$n=\infty$']
# legend_name= [r'$\theta_i=0.1$',r'$\theta_i$ bestfit',r'$n=\infty$']
legend_name= [r'$\theta_i$ bestfit',r'$\theta_i=3$',]
alpha = [0.5,0.3]
# legend_name= ['TT','TTTEEE']
var_figname = 'nalp'
dashes_array = [[10,0.0000000000001,10,0.0000000000001],[3,1,3,1],[5,2,1,2],[5,2,1,2],[5,2,1,2],[5,2,1,2]]
#
#############################################
#
# Fixed settings
#
common_settings = {'output':'tCl,pCl,lCl,mPk',
                   'lensing':'yes',
                   'format':'camb',
                   'N_ncdm': 1,
                   'deg_ncdm': 1,
                   'm_ncdm': 0.06,
                   'T_ncdm': 0.71611,
                   'N_ur':  2.0328,
                   # # LambdaCDM parameters
                   # # 'h':0.67556,
                   # '100*theta_s':1.04077,
                   # 'omega_b':0.02225,
                   # # 'omega_cdm':0.12038,
                   # # 'Omega_cdm':0.266,
                   # 'Omega_cdm':0.266,
                   # # 'A_s':2.215e-9,
                   # 'ln10^{10}A_s':3.094,
                   # 'n_s':0.9645,
                   # 'tau_reio':0.079,
                   # 'N_ur':3.046,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   #'YHe':0.246,
                   # other output and precision parameters
                   'P_k_max_1/Mpc':3.0,
                   'l_switch_limber':9,
                   'input_verbose':1,
                   'background_verbose':1,
                   'thermodynamics_verbose':1,
                   'perturbations_verbose':1,
                   'transfer_verbose':1,
                   'primordial_verbose':1,
                   'spectra_verbose':1,
                   'nonlinear_verbose':1,
                   'lensing_verbose':1,
                   'output_verbose':1}
                   #,'temperature contributions':'tsw'}
                   #'background_verbose':1}
#
# arrays for output
#
kvec = np.logspace(-4,np.log10(3),1000)
legarray = []
twopi = 2.*math.pi
#
# Create figures
#
# fig_Pk, ax_Pk = plt.subplots()
# fig_TT, ax_TT = plt.subplots()
# fig_EE, ax_EE = plt.subplots()
# fig_PP, ax_PP = plt.subplots()

# fig_TT, (ax_TT, ax_EE, ax_PP, ax_Pk) = plt.subplots(4,1, sharex=False,gridspec_kw=dict(height_ratios=[1,1,1,1]),figsize=(35,40))
# fig_TT, (ax_TT, ax_Pk) = plt.subplots(2,1, sharex=False,gridspec_kw=dict(height_ratios=[1,1]),figsize=(15,20))

####OLD#####
# fig_TT, (ax_TT, ax_EE) = plt.subplots(2,1, sharex=True,gridspec_kw=dict(height_ratios=[1,1]),figsize=(15,20))
# fig_TT.subplots_adjust(hspace=0)
# ax_TT.tick_params(axis='x',labelsize=30,length=20,which='both',width=2,direction='inout')
# ax_TT.tick_params(axis='y',labelsize=30,length=20,which='both',width=2,direction='inout')
# ax_EE.tick_params(axis='x',labelsize=30,length=20,which='both',width=2,direction='inout')
# ax_EE.tick_params(axis='y',labelsize=30,length=20,which='both',width=2,direction='inout')
# # ax_Pk.tick_params(axis='x',labelsize=30,length=20,which='both',width=2,direction='inout')
# # ax_Pk.tick_params(axis='y',labelsize=30,length=20,which='both',width=2,direction='inout')
# ax_TT.tick_params(axis='x',labelsize=30,length=10,which='minor',direction='inout')
# ax_TT.tick_params(axis='y',labelsize=30,length=10,which='minor',direction='inout')
# ax_EE.tick_params(axis='x',labelsize=30,length=10,which='minor',direction='inout')
# ax_EE.tick_params(axis='y',labelsize=30,length=10,which='minor',direction='inout')
# ######
# ax_Pk.tick_params(axis='x',labelsize=30,length=10,which='minor',direction='inout')
# ax_Pk.tick_params(axis='y',labelsize=30,length=10,which='minor',direction='inout')


plot_pk = False

if plot_pk == True:
    ax_Pk = plt.subplot(313)
    ax_EE_log = plt.subplot(312)
    ax_TT_log = plt.subplot(311)
    plt.setp(ax_Pk.get_xticklabels(), fontsize=15)
    plt.setp(ax_Pk.get_yticklabels(), fontsize=15)
    ax_Pk.axis([1.e-4,1.,-0.3,0.3])
    ax_Pk.set_xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$',fontsize=15)
    ax_Pk.set_ylabel(r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=15)
if plot_pk == False:
    ax_EE_log = plt.subplot(313)
    ax_TE_log = plt.subplot(312)
    ax_TT_log = plt.subplot(311)

plt.subplots_adjust(hspace=0.4)
plt.setp(ax_TT_log.get_xticklabels(), fontsize=15)
plt.setp(ax_TE_log.get_xticklabels(), fontsize=15)
plt.setp(ax_EE_log.get_xticklabels(), fontsize=15)
plt.setp(ax_TT_log.get_yticklabels(), fontsize=15)
plt.setp(ax_TE_log.get_yticklabels(), fontsize=15)
plt.setp(ax_EE_log.get_yticklabels(), fontsize=15)
# plt.subplots_adjust(hspace=0)
divider = make_axes_locatable(ax_TT_log)
ax_TT_lin = divider.append_axes("right", size=5, pad=0)

divider = make_axes_locatable(ax_TE_log)
ax_TE_lin = divider.append_axes("right", size=5, pad=0)

divider = make_axes_locatable(ax_EE_log)
ax_EE_lin = divider.append_axes("right", size=5, pad=0)


plt.setp(ax_TT_lin.get_yticklabels(), fontsize=15)
plt.setp(ax_TE_lin.get_yticklabels(), fontsize=15)
plt.setp(ax_EE_lin.get_yticklabels(), fontsize=15)
plt.setp(ax_TT_lin.get_xticklabels(), fontsize=15)
plt.setp(ax_TE_lin.get_xticklabels(), fontsize=15)
plt.setp(ax_EE_lin.get_xticklabels(), fontsize=15)

M = Class()
M.set(common_settings)


# ###LCDM BESTFIT ours####
# M.set({
# '100*theta_s':1.042129e+00,
# 'omega_b':2.239439e-02,
# 'omega_cdm':1.176791e-01,
# # 'ln10^{10}A_s':3.091635e+00,
# 'A_s':2.14e-9,
# 'n_s':9.687e-01,
# 'tau_reio':0.0684,
# })
##LCDM BESTFIT Planck 2018####
# M.set({
# '100*theta_s':1.040909,
# 'omega_b':0.022383,
# 'omega_cdm':0.12011,
# 'ln10^{10}A_s':3.0448,
# 'n_s':0.96605,
# 'tau_reio':0.0543,
# })
# #OLD LCDM
M.set({
'100*theta_s':1.04077,
'omega_b':0.02225,
'omega_cdm':0.1198,
'ln10^{10}A_s':3.094,
'n_s':0.9645,
'tau_reio':0.079,
'do_shooting':'yes'
})
#

M.compute()
all_k = M.get_perturbations()  # this potentially constains scalars/tensors and all k values
conversion = pow(2.7255*1.e6,2)
# one_k = all_k['scalar'][0]     # this contains only the scalar perturbations for the requested k values


clM = M.lensed_cl(2500)
ll_LCDM = clM['ell'][2:]
clTT_LCDM = clM['tt'][2:]
clEE_LCDM = clM['ee'][2:]
clTE_LCDM = clM['te'][2:]
clPP_LCDM = clM['pp'][2:]
# ax_TT.semilogx(ll_LCDM,(ll_LCDM)*(ll_LCDM+1)/(2*np.pi)*(clTT_LCDM),'k',lw=5)
# ax_EE.semilogx(ll_LCDM,(ll_LCDM)*(ll_LCDM+1)/(2*np.pi)*(clEE_LCDM),'k',lw=5)
# ax_PP.semilogx(ll_LCDM,(ll_LCDM)*(ll_LCDM+1)/(2*np.pi)*(clPP_LCDM),'k',lw=5)

pkM_LCDM = []

for k in kvec:
    pkM_LCDM.append(M.pk(k,0.))


#     i=i+1
fTT = interp1d(ll_LCDM,clTT_LCDM)
fEE = interp1d(ll_LCDM,clEE_LCDM)
fTE = interp1d(ll_LCDM,clTE_LCDM)

M.struct_cleanup()
for i in range(var_num):
    #
    # deal with varying parameters:
    #
    dashes=dashes_array[i]
    #
    l=[]
    print ' * Compute with %s=%d'%(var_name1,n[i])
    #
    # deal with colors and legends
    #

    # var_legend = r'$n_{\rm a}=%d$'%(var_n)
    var_legend = legend_name[i]
    var_alpha = alpha[i]
    var_color = color[i]
    # var_alpha = 1.*i/(var_num-1.)
    legarray.append(var_legend)
    # if i == var_num-1:
    #     legarray.append(var_legend)
    #
    # call CLASS
    #
    M = Class()
    M.set(common_settings)
    # if i == 0:
    #     M.set({
    #             'omega_cdm': 1.219020e-01,
    #             'omega_b': 2.263778e-02,
    #             '100*theta_s': 1.041373e+00,
    #             'tau_reio': 7.816028e-02,
    #             'ln10^{10}A_s': 3.093710e+00,
    #             # 'A_s': 2.11839e-9,
    #             # 'n_s':9.732537e-01,
    #             'n_s':9.816505e-01,
    #             'N_ncdm': 1,
    #             'deg_ncdm': 3,
    #             'm_ncdm': 0.02,
    #             'T_ncdm': 0.71611,
    #             'N_ur':  2.997612e-01})
    if i == 0:
        # M.set({'scf_potential' : 'axion',
        #         'n_axion' : 2,
        #         'scf_parameters' : '0.1,0',
    	# 		'scf_tuning_index':0,
        #         'log10_fraction_axion_ac':-1.080,
        #         'log10_axion_ac':-3.46492,
        #         'omega_cdm':0.129621,
        #         'omega_b': 2.26488e-2,
        #         '100*theta_s':1.0416,
        #         'tau_reio':0.067376,
        #         'A_s':2.19616e-9,
        #         'n_s':0.97937,
    	# 		'scf_evolve_like_axionCAMB':'no',
    	# 		'do_shooting':'yes',
    	# 		'do_shooting_scf':'yes',
    	# 		'use_big_theta_scf':'yes',
    	# 		'scf_has_perturbations':'yes',
    	# 		'attractor_ic_scf':'no',
    	# 		'adptative_stepsize':100,
    	# 		'scf_evolve_as_fluid':'no'})
        #Theta large
        M.set({'scf_potential' : 'axion',
                'n_axion' : 2,
                'scf_parameters' : '2.11954,0',
    			'scf_tuning_index':0,
                'log10_fraction_axion_ac':-0.9298,
                'log10_axion_ac':-3.50979,
                'omega_cdm':0.131354,
                'omega_b': 2.23759e-2,
                # '100*theta_s':1.04119,
                'H0':72.4265,
                'tau_reio':0.0795046,
                'A_s':2.24017e-9,
                'n_s':0.976422,
    			'scf_evolve_like_axionCAMB':'no',
    			'do_shooting':'yes',
    			'do_shooting_scf':'yes',
    			'use_big_theta_scf':'yes',
    			'scf_has_perturbations':'yes',
    			'attractor_ic_scf':'no',
    			'adptative_stepsize':100,
    			'scf_evolve_as_fluid':'no',
                'security_small_Omega_scf':1e-3})
        #small Theta
        # M.set({'scf_potential' : 'axion',
        #         'n_axion' : 3,
        #         'scf_parameters' : '0.1,0',
    	# 		'scf_tuning_index':0,
        #         'log10_fraction_axion_ac':-1.29877,
        #         'log10_axion_ac':-3.43756,
        #         'omega_cdm':0.12505,
        #         'omega_b': 2.24228e-2,
        #         '100*theta_s':1.0415,
        #         'tau_reio':0.063995,
        #         'A_s':2.16465e-9,
        #         'n_s':0.97584,
    	# 		'scf_evolve_like_axionCAMB':'no',
    	# 		'do_shooting':'yes',
    	# 		'do_shooting_scf':'yes',
    	# 		'use_big_theta_scf':'yes',
    	# 		'scf_has_perturbations':'yes',
    	# 		'attractor_ic_scf':'no',
    	# 		'adptative_stepsize':100,
    	# 		'scf_evolve_as_fluid':'no'})
        # M.set({'scf_potential' : 'axion',
        #         'n_axion' : 3,
        #         'scf_parameters' : '0.1,0',
        #         'scf_tuning_index':0,
        #         'log10_fraction_axion_ac':-0.8788,
        #         'log10_axion_ac':-3.55314,
        #         'omega_cdm':0.132745,
        #         'omega_b': 2.26145e-02,
        #         # '100*theta_s':1.04262e+00,
        #         'H0':73.4205,
        #         'tau_reio':0.075,
        #         'A_s':2.23266e-9,
        #         'n_s':0.988002,
        #         'scf_evolve_like_axionCAMB':'no',
        #         'do_shooting':'yes',
        #         'do_shooting_scf':'yes',
        #         'use_big_theta_scf':'yes',
        #         'scf_has_perturbations':'yes',
        #         'attractor_ic_scf':'no',
        #         'adptative_stepsize':500,
        #         'scf_evolve_as_fluid':'no'})
    elif i == 1:
        #small Theta
        # M.set({'scf_potential' : 'axion',
        #         'n_axion' : 3,
        #         'scf_parameters' : '0.1,0',
    	# 		'scf_tuning_index':0,
        #         'log10_fraction_axion_ac':-1.29877,
        #         'log10_axion_ac':-3.43756,
        #         'omega_cdm':0.12505,
        #         'omega_b': 2.24228e-2,
        #         '100*theta_s':1.0415,
        #         'tau_reio':0.063995,
        #         'A_s':2.16465e-9,
        #         'n_s':0.97584,
    	# 		'scf_evolve_like_axionCAMB':'no',
    	# 		'do_shooting':'yes',
    	# 		'do_shooting_scf':'yes',
    	# 		'use_big_theta_scf':'yes',
    	# 		'scf_has_perturbations':'yes',
    	# 		'attractor_ic_scf':'no',
    	# 		'adptative_stepsize':100,
    	# 		'scf_evolve_as_fluid':'no'})
        # #large Theta
        # M.set({'scf_potential' : 'axion',
        #         'n_axion' : 3,
        #         'scf_parameters' : '2.74031,0',
        #         'scf_tuning_index':0,
        #         'log10_fraction_axion_ac':-0.8788,
        #         'log10_axion_ac':-3.55314,
        #         'omega_cdm':0.132745,
        #         'omega_b': 2.26145e-02,
        #         # '100*theta_s':1.04262e+00,
        #         # '100*theta_s':1.04262e+00,
        #         'H0':72.2205,
        #         'tau_reio':0.075,
        #         'A_s':2.23266e-9,
        #         'n_s':0.988002,
        #         'scf_evolve_like_axionCAMB':'no',
        #         'do_shooting':'yes',
        #         'do_shooting_scf':'yes',
        #         'use_big_theta_scf':'yes',
        #         'scf_has_perturbations':'yes',
        #         'attractor_ic_scf':'no',
        #         'adptative_stepsize':500,
        #         'scf_evolve_as_fluid':'no'})
        M.set({'scf_potential' : 'axion',
                'n_axion' : 2,
                'scf_parameters' : '3,0',
    			'scf_tuning_index':0,
                'log10_fraction_axion_ac':-0.9298,
                'log10_axion_ac':-3.50979,
                'omega_cdm':0.131354,
                'omega_b': 2.23759e-2,
                # '100*theta_s':1.04119,
                'H0':72.4265,
                'tau_reio':0.0795046,
                'A_s':2.24017e-9,
                'n_s':0.976422,
    			'scf_evolve_like_axionCAMB':'no',
    			'do_shooting':'yes',
    			'do_shooting_scf':'yes',
    			'use_big_theta_scf':'yes',
    			'scf_has_perturbations':'yes',
    			'attractor_ic_scf':'no',
    			'adptative_stepsize':100,
    			'scf_evolve_as_fluid':'no',
                'security_small_Omega_scf':1e-3})
                # 2.258053e-02	 1.298650e-01	 9.880137e-01	 2.176655e-09	 6.833666e-02	 1.041246e+00	-3.736778e+00 4.717743e-09	 2.791761e+00


    M.compute()
    #
    # get Cls
    #
    # clM = M.raw_cl(2500
    clM = M.lensed_cl(2500)
    ll = clM['ell'][2:]
    clTT = clM['tt'][2:]
    clTE = clM['te'][2:]
    clEE = clM['ee'][2:]
    bg = M.get_background()
    Omega_scf = bg['(.)Omega_scf']
    H = bg['H [1/Mpc]']
    Da = bg['ang.diam.dist.']
    z = bg['z']
    tau = bg['conf. time [Mpc]']
    tau_interp = interp1d(z,tau)
    tau_0 =  tau_interp(0)
    derived = M.get_current_derived_parameters(['z_rec','tau_rec'])
    tau_rec =  derived['tau_rec']
    for j in range(len(Omega_scf)):
        if Omega_scf[j] > 0.01:
            k = H[j]/(1+z[j])
            # l.append(k*Da[j])
            l.append(k*(tau_0-tau_rec))
    print 450./(tau_0-tau_rec)
    #    store P(k) for common k values
    lmax = max(l)
    lmin = min(l)
    l=[]
    l.append(lmin)
    l.append(2*lmax)
    if plot_pk == True:
        pkM = []
        for k in kvec:
            pkM.append(M.pk(k,0.))
        ax_Pk.semilogx(kvec,(np.array(pkM)-np.array(pkM_LCDM))/np.array(pkM_LCDM),var_color,lw=2,dashes=dashes,label=var_legend)
    print l,(clTE)-(clTE_LCDM)
    # ax_TT_lin.fill_between(l,-2,2,alpha=var_alpha, facecolor=var_color,
    #                  linewidth=0)
    # ax_TE_lin.fill_between(l,-2,2,alpha=var_alpha, facecolor=var_color,
    #                  linewidth=0)
    # ax_EE_lin.fill_between(l,-2,2,alpha=var_alpha, facecolor=var_color,
    #                  linewidth=0)
    sigma_CV = np.sqrt(clTT_LCDM*clEE_LCDM+clTE_LCDM**2)
    print ((clTE)-(clTE_LCDM))/sigma_CV,np.sqrt(1/(2*ll+1)),np.sqrt(clTT_LCDM*clEE_LCDM+clTE_LCDM**2)
    ax_TT_lin.plot(ll,(clTT-clTT_LCDM)/clTT_LCDM,var_color,lw=2,dashes=dashes,label=var_legend)
    ax_TE_lin.plot(ll,((clTE)-(clTE_LCDM))/sigma_CV,var_color,lw=2,dashes=dashes,label=var_legend)
    ax_EE_lin.plot(ll,(clEE-clEE_LCDM)/clEE_LCDM,var_color,lw=2,dashes=dashes,label=var_legend)
    ax_TT_log.plot(ll,(clTT-clTT_LCDM)/clTT_LCDM,var_color,lw=2,dashes=dashes,label=var_legend)
    ax_TE_log.plot(ll,((clTE)-(clTE_LCDM))/sigma_CV,var_color,lw=2,dashes=dashes,label=var_legend)
    ax_EE_log.plot(ll,(clEE-clEE_LCDM)/clEE_LCDM,var_color,lw=2,dashes=dashes,label=var_legend)
    M.struct_cleanup()

#######TRIS SCRIPT####

#fig, axMain = plt.subplots(figsize=(10, 3))

# axMain.errorbar(ellTT, LCDM, yerr=Plerr, fmt='.')
###OLDSCALES###
# ax_TT_log.set_xscale('log')
# ax_TT_log.set_xlim((2,29))
# ax_TT_log.set_ylim((-0.2,0.2))
# ax_EE_log.set_xscale('log')
# ax_EE_log.set_xlim((2,29))
# ax_EE_log.set_ylim((-0.7,0.7))
#
# ax_TT_log.set_xticks([1., 10.])
# ax_EE_log.set_xticks([1., 10.])
# ax_TT_log.set_yticks([-0.2,-0.1,0.,0.1,0.2])
# ax_EE_log.set_yticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
#
# ax_TT_log.spines['right'].set_visible(False)
# ax_EE_log.yaxis.set_ticks_position('left')
# ax_TT_log.spines['right'].set_visible(False)
# ax_EE_log.yaxis.set_ticks_position('left')
#
# ax_TT_log.tick_params('both', length=10, width=1, which='major')
# ax_TT_log.tick_params('both', length=5, width=1, which='minor')
# ax_EE_log.tick_params('both', length=10, width=1, which='major')
# ax_EE_log.tick_params('both', length=5, width=1, which='minor')
#
# # axLin.plot(ellCIP, CIPTT,color="cyan", linestyle='-',lw=2)
# # axLin.plot(ellCIP, CIPTT2,color="red", linestyle='-',lw=2)
#
# ax_TT_lin.set_xscale('linear')
# ax_TT_lin.set_xlim((30, 2000))
# ax_TT_lin.set_ylim((-0.05,0.05))
# ax_TT_lin.tick_params('both', length=10, width=1, which='major')
# ax_TT_lin.tick_params('both', length=5, width=1, which='minor')
# ax_EE_lin.set_xscale('linear')
# ax_EE_lin.set_xlim((30, 2000))
# ax_EE_lin.set_ylim((-0.05,0.05))
# ax_EE_lin.tick_params('both', length=10, width=1, which='major')
# ax_EE_lin.tick_params('both', length=5, width=1, which='minor')
#
# major_ticks = [30, 500, 1000, 1500, 2000]
#
# ax_TT_lin.set_yticks([-0.04,-0.02,0,0.02,0.04])
# # Removes bottom axis line
# ax_TT_lin.spines['left'].set_visible(False)
# ax_TT_lin.yaxis.set_ticks_position('right')
#
# ax_TT_lin.axvline(30, color='black', linestyle='--',lw=2)
#
# ax_TT_lin.set_xticks(major_ticks)
# ax_EE_lin.set_xticks(major_ticks)
# ax_EE_lin.set_yticks([-0.04,-0.02,0,0.02,0.04])
#
# # Removes bottom axis line
# ax_EE_lin.spines['left'].set_visible(False)
# ax_EE_lin.yaxis.set_ticks_position('right')
#
# ax_EE_lin.axvline(30, color='black', linestyle='--',lw=2)



###NEW SCALES###
ax_TT_log.set_xscale('linear')
ax_TT_log.set_xlim((1,29))
ax_TT_log.set_ylim((-1.5,1.5))
ax_EE_log.set_xscale('linear')
ax_EE_log.set_xlim((1,29))
ax_EE_log.set_ylim((-10,10))
ax_TE_log.set_xscale('linear')
ax_TE_log.set_xlim((1,29))
# ax_TE_log.set_ylim((-10,10))

ax_TT_log.set_xticks([2., 10.])
ax_EE_log.set_xticks([2., 10.])
ax_TE_log.set_xticks([2., 10.])
ax_TT_log.set_yticks([-0.8,-0.4,0.,0.4,0.8])
ax_EE_log.set_yticks([-7.5,-5,-2.5,0.0,2.5,5,7.5])
ax_TE_log.set_yticks([-5,-2.5,0.0,2.5,5])

ax_TT_log.spines['right'].set_visible(False)
ax_EE_log.yaxis.set_ticks_position('left')
ax_TE_log.yaxis.set_ticks_position('left')
ax_TT_log.spines['right'].set_visible(False)
ax_EE_log.yaxis.set_ticks_position('left')
ax_TE_log.yaxis.set_ticks_position('left')
ax_TT_log.tick_params('both', length=10, width=1, which='major')
ax_TT_log.tick_params('both', length=5, width=1, which='minor')
ax_EE_log.tick_params('both', length=10, width=1, which='major')
ax_EE_log.tick_params('both', length=5, width=1, which='minor')
ax_TE_log.tick_params('both', length=10, width=1, which='major')
ax_TE_log.tick_params('both', length=5, width=1, which='minor')

# axLin.plot(ellCIP, CIPTT,color="cyan", linestyle='-',lw=2)
# axLin.plot(ellCIP, CIPTT2,color="red", linestyle='-',lw=2)

ax_TT_lin.set_xscale('linear')
ax_TT_lin.set_xlim((30, 2510))
ax_TT_lin.set_ylim((-0.12,0.12))
ax_TT_lin.tick_params('both', length=10, width=1, which='major',labelsize=15)
ax_TT_lin.tick_params('both', length=5, width=1, which='minor',labelsize=15)
ax_EE_lin.set_xscale('linear')
ax_EE_lin.set_xlim((30, 2010))
ax_EE_lin.set_ylim((-0.12,0.12))
ax_EE_lin.tick_params('both', length=10, width=1, which='major',labelsize=15)
ax_EE_lin.tick_params('both', length=5, width=1, which='minor',labelsize=15)
ax_TE_lin.set_xscale('linear')
# ax_TE_lin.set_yscale('log')
ax_TE_lin.set_xlim((30, 2010))
ax_TE_lin.set_ylim((-0.12,0.12))
ax_TE_lin.tick_params('both', length=10, width=1, which='major',labelsize=15)
ax_TE_lin.tick_params('both', length=5, width=1, which='minor',labelsize=15)

major_ticks = [30, 500, 1000, 1500, 2000,2500]

ax_TT_lin.set_yticks([-0.08,-0.04,0,0.04,0.08])

# Removes bottom axis line
ax_TT_lin.spines['left'].set_visible(False)
ax_TT_lin.yaxis.set_ticks_position('right')

ax_TT_lin.axvline(30, color='black', linestyle='--',lw=2)

ax_TT_lin.set_xticks(major_ticks)
ax_EE_lin.set_xticks(major_ticks)
ax_EE_lin.set_yticks([-0.08,-0.04,0,0.04,0.08])

# Removes bottom axis line
ax_EE_lin.spines['left'].set_visible(False)
ax_EE_lin.yaxis.set_ticks_position('right')

ax_EE_lin.axvline(30, color='black', linestyle='--',lw=2)
ax_TE_lin.set_xticks(major_ticks)
ax_TE_lin.set_yticks([-0.08,-0.04,0,0.04,0.08])

# Removes bottom axis line
ax_TE_lin.spines['left'].set_visible(False)
ax_TE_lin.yaxis.set_ticks_position('right')

ax_TE_lin.axvline(30, color='black', linestyle='--',lw=2)
#
# ###############
# ###############
# axMain = plt.subplot(312)
# #FigSize = (15,15)
# #fig, axMain = plt.subplots(figsize=(10, 3))
# axMain.plot(ellCIP, CIPEE ,color="red", linestyle='-',lw=2)
# axMain.errorbar(ellEE, LCDMEE, yerr=PlerrEE, fmt='.')
#
# axMain.set_xscale('log')
# axMain.set_xlim((2,29))
# axMain.set_ylim((-.5,.5))
#
# axMain.tick_params('both', length=10, width=1, which='major')
# axMain.tick_params('both', length=5, width=1, which='minor')
#
# axMain.set_xticks([1, 10])
#
# axMain.spines['right'].set_visible(False)
# axMain.yaxis.set_ticks_position('left')
# axMain.set_ylabel(r'$\Delta \mathcal{D}^{\rm EE}_{l}\ (\mu{\rm K})$', fontsize = 22)
#
# divider = make_axes_locatable(axMain)
# axLin = divider.append_axes("right", size=5, pad=0)
# axLin.plot(ellCIP, CIPEE ,color="red", linestyle='-',lw=2)
# axLin.errorbar(ellEE, LCDMEE, yerr=PlerrEE, fmt='.')
#
# axLin.tick_params('both', length=10, width=1, which='major')
# axLin.tick_params('both', length=5, width=1, which='minor')
#
# axLin.set_xscale('linear')
# axLin.set_xlim((30, 2000))
# axLin.set_ylim((-4,4))
#
# major_ticks = [30, 500, 1000, 1500, 2000]
#
# axLin.set_xticks(major_ticks)
#
# # Removes bottom axis line
# axLin.spines['left'].set_visible(False)
# axLin.yaxis.set_ticks_position('right')
#
# axLin.axvline(30, color='black', linestyle='--',lw=5)
#
# ###############
# ###############
# axMain = plt.subplot(313)
#
# #FigSize = (15,15)
# #fig, axMain = plt.subplots(figsize=(10, 3))
# axMain.plot(ellCIP, CIPTE ,color="red", linestyle='-',lw=2)
# axMain.errorbar(ellTE, LCDMTE, yerr=PlerrTE, fmt='.')
#
# axMain.set_xscale('log')
# axMain.set_xlim((2,29))
# axMain.set_ylim((-20,20))
#
# axMain.tick_params('both', length=10, width=1, which='major')
# axMain.tick_params('both', length=5, width=1, which='minor')
#
# axMain.set_xticks([1, 10])
#
# axMain.spines['right'].set_visible(False)
# axMain.yaxis.set_ticks_position('left')
# axMain.set_ylabel(r'$\Delta \mathcal{D}^{\rm TE}_{l} \ (\mu{\rm K})$', fontsize = 22)
# axMain.set_xlabel(r'$l$', fontsize = 22)
#
#
#
# divider = make_axes_locatable(axMain)
# axLin = divider.append_axes("right", size=5, pad=0)
# axLin.plot(ellCIP, CIPTE ,color="red", linestyle='-',lw=2)
# axLin.errorbar(ellTE, LCDMTE, yerr=PlerrTE, fmt='.')
#
# axLin.tick_params('both', length=10, width=1, which='major')
# axLin.tick_params('both', length=5, width=1, which='minor')
#
# axLin.set_xscale('linear')
# axLin.set_xlim((30, 2000))
# axLin.set_ylim((-20,20))
#
# major_ticks = [30, 500, 1000, 1500, 2000]
#
# axLin.set_xticks(major_ticks)
#
# # Removes bottom axis line
# axLin.spines['left'].set_visible(False)
# axLin.yaxis.set_ticks_position('right')
#
# axLin.axvline(30, color='black', linestyle='--',lw=5)
#
# plt.tight_layout()
# plt.show
# plt.savefig('Planck_diff.pdf')

######
# ax_Pk.semilogx(kvec,(np.array(pkM_LCDM)-np.array(pkM_LCDM))/np.array(pkM_LCDM),'k',lw=5,label=var_legend)

    #
#

# # ax_Pk.text(0.007,0.11,r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=65)###ac1em5
# ax_Pk.text(0.007,0.05,r'$\frac{\Delta P(k)}{P(k)^{\Lambda{\rm CDM}}}$',fontsize=65)
# ax_PP.set_xlabel(r'$\ell$',fontsize=65)
# ax_PP.text(50,0.03,r'$\frac{\Delta C_\ell^{\phi\phi}}{C_\ell^{\phi\phi}({\Lambda{\rm CDM}})}$',fontsize=65)
# ax_PP.text(50,0.06,r'$\frac{\Delta C_\ell^{\phi\phi}}{C_\ell^{\phi\phi}({\Lambda{\rm CDM}})}$',fontsize=65)



# ###NEW: GREY BANDS AROUND 0###
# ax_TT_lin.fill_between(l_TT_high,-err_TT_high/Dl_TT_high,err_TT_high/Dl_TT_high,facecolor='grey',alpha=0.2)
# ax_TT_log.fill_between(l_TT_low,-err_TT_low/Dl_TT_low,err_TT_low/Dl_TT_low,facecolor='grey',alpha=0.2)
# ax_EE_lin.fill_between(l_EE_high,-err_EE_high/Dl_EE_high,err_EE_high/Dl_EE_high,facecolor='grey',alpha=0.2)
# ax_EE_log.fill_between(l_EE_low,-err_EE_low/Dl_EE_low,err_EE_low/Dl_EE_low,facecolor='grey',alpha=0.2)


# ###ERROR BARS W/R LCDM###
T_cmb = 2.7225
factor1 = l_TT_high*(l_TT_high+1)/2./np.pi;
conversion1 = 1/(factor1*(T_cmb*1.e6)**2)
factor2 = l_TT_low*(l_TT_low+1)/2./np.pi;
conversion2 = 1/(factor2*(T_cmb*1.e6)**2)
factor3 = l_EE_high*(l_EE_high+1)/2./np.pi;
conversion3 = 1/(factor3*(T_cmb*1.e6)**2)
factor4 = l_EE_low*(l_EE_low+1)/2./np.pi;
conversion4 = 1/(factor4*(T_cmb*1.e6)**2)
factor5 = l_TE_high*(l_TE_high+1)/2./np.pi;
conversion5 = 1/(factor3*(T_cmb*1.e6)**2)
factor6 = l_TE_low*(l_TE_low+1)/2./np.pi;
conversion6 = 1/(factor4*(T_cmb*1.e6)**2)

print Dl_TT_high*conversion, fTT(l_TT_high)
ax_TT_lin.errorbar(l_TT_high, Dl_TT_high*conversion1/fTT(l_TT_high)-1, yerr=err_TT_high*conversion1/fTT(l_TT_high), fmt='.')
ax_TT_log.errorbar(l_TT_low, Dl_TT_low*conversion2/fTT(l_TT_low)-1, yerr=err_TT_low*conversion2/fTT(l_TT_low), fmt='.')
ax_EE_lin.errorbar(l_EE_high, Dl_EE_high*conversion3/fEE(l_EE_high)-1, yerr=err_EE_high*conversion3/fEE(l_EE_high), fmt='.')
ax_EE_log.errorbar(l_EE_low, Dl_EE_low*conversion4/fEE(l_EE_low)-1, yerr=err_EE_low*conversion4/fEE(l_EE_low), fmt='.')
ax_TE_lin.errorbar(l_TE_high, (Dl_TE_high*conversion5-fTE(l_TE_high))/np.sqrt(fTT(l_TE_high)*fEE(l_TE_high)+fTE(l_TE_high)**2), yerr=err_TE_high*conversion5/np.sqrt(fTT(l_TE_high)*fEE(l_TE_high)+fTE(l_TE_high)**2), fmt='.')
ax_TE_log.errorbar(l_TE_low, (Dl_TE_low*conversion6-fTE(l_TE_low))/np.sqrt(fEE(l_TE_low)*fTT(l_TE_low)+fTE(l_TE_low)**2), yerr=err_TE_low*conversion6/np.sqrt(fTT(l_TE_low)*fEE(l_TE_low)+fTE(l_TE_low)**2), fmt='.')



# ax_EE.fill_between(l_EE_low,-err_EE_low/Dl_EE_low,err_EE_low/Dl_EE_low,facecolor='grey',alpha=0.2)
#
#
# ax_TT.axis([2,2500,-0.06,0.06])
# ax_TT.set_xlabel(r'$\ell$',fontsize=35)
# ax_TT_lin.text(1050,0.02,r'$\frac{\Delta C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}(\Lambda{\rm CDM})}$',fontsize=20)
# ax_TT_log.set_ylabel(r'$\Delta C_\ell^\mathrm{TT}/C_\ell^\mathrm{TT}(\Lambda{\rm CDM})$',fontsize=17)
ax_TT_log.set_ylabel(r'$\frac{\Delta C_\ell^\mathrm{TT}}{C_\ell^\mathrm{TT}}$',fontsize=19)

# ax_Pk.legend(frameon=False,prop={'size':30},loc='upper left',borderaxespad=0.)
ax_TT_lin.set_xlabel(r'$\ell$',fontsize=20,labelpad=-20)
ax_EE_lin.legend(frameon=False,prop={'size':12},loc='upper right',borderaxespad=0.)

# ax_EE.axis([2,2500,-0.06,0.06])
ax_EE_lin.set_xlabel(r'$\ell$',fontsize=20,labelpad=-20)
# ax_EE_lin.text(200,-0.1,r'$\frac{\Delta C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}(\Lambda{\rm CDM})}$',fontsize=20)
# ax_EE_log.set_ylabel(r'$\Delta C_\ell^\mathrm{EE}/C_\ell^\mathrm{EE}(\Lambda{\rm CDM})$',fontsize=19)
ax_EE_log.set_ylabel(r'$\frac{\Delta C_\ell^\mathrm{EE}}{C_\ell^\mathrm{EE}}$',fontsize=19)

ax_TE_lin.legend(frameon=False,prop={'size':12},loc='upper right',borderaxespad=0.)

# ax_TE.axis([2,2500,-0.06,0.06])
ax_TE_lin.set_xlabel(r'$\ell$',fontsize=20,labelpad=-20)
# ax_TE_lin.text(200,-0.1,r'$\frac{\Delta C_\ell^\mathrm{TE}}{C_\ell^\mathrm{TE}(\Lambda{\rm CDM})}$',fontsize=20)
# ax_TE_log.set_ylabel(r'$\Delta C_\ell^\mathrm{TE}/C_\ell^\mathrm{TE}(\Lambda{\rm CDM})$',fontsize=19)
ax_TE_log.set_ylabel(r'$\frac{\Delta C_\ell^\mathrm{TE}}{\sqrt{C_\ell^\mathrm{EE}*C_\ell^\mathrm{TT}+(C_\ell^\mathrm{TE})^2}}$',fontsize=19)


# fig_TT.savefig('spectra_nalp_zc1e3_ErrorPlanck.pdf')
# fig_TT.savefig('spectra_nalp_zc1e5_ErrorPlanck.pdf', bbox_inches='tight')
# fig_TT.savefig('spectra_nalp_n1_ErrorPlanck.pdf', bbox_inches='tight')
# figxe.savefig('xe_nalp_zc1e3_ErrorPlanck.pdf')



def binned_cosmic_variance (result,l_ini,width):
    central_l = l_ini+width/2
    weight_total = 0
    result = 0
    Clb = 0
    # for i in range(0,int(width)):
    #     weight_total += (l_ini+i)*(l_ini+i+1)
    #     result += np.sqrt(2/(2*(l_ini+float(i))+1))*(l_ini+float(i))*(l_ini+float(i)+1)
    # print l_ini,l_ini+width,result, weight_total
    # return result/weight_total/np.sqrt(2)/np.pi
    for i in range(0,int(width)):
        weight_total += (l_ini+i)*(l_ini+i+1)
        result += 2/(2*(l_ini+float(i))+1)*(l_ini+float(i))*(l_ini+float(i)+1)*(l_ini+float(i))*(l_ini+float(i)+1)*fTT(l_ini+i)*fTT(l_ini+i)
        Clb += (l_ini+float(i))*(l_ini+float(i)+1)*fTT(l_ini+i)
    # print l_ini,l_ini+width,np.sqrt(result), Clb, np.sqrt(result)/Clb
    return np.sqrt(result)/Clb

#
#
#
#
l_min = 2.;
l_max = 2000;
n_step = 100.;
j=0.
step = l_min
width= 25
# while j < n_step:
while step < l_max:
        result = 0.0
        if step < 29:
            width = 1
            step = l_min+j*width
            j+=1
            if step == 29:
                j = 0
                l_min = 30
        else:
            width = 30
            step = l_min+j*width
            j+=1

        # step = l_min*(l_max/l_min)**(j/n_step)
        # step_plus_1 = l_min*(l_max/l_min)**((j+1)/n_step)
        # print step
        # print int(step), int(step_plus_1)
        # width = (step_plus_1) - (step)
        ax_TT_lin.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        ax_TT_log.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        ax_EE_lin.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        ax_EE_log.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        ax_TE_lin.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        ax_TE_log.add_patch(
            patches.Rectangle(
                (int(step), -1*binned_cosmic_variance(result,int(step),width)),   # (x,y)
                width,          # width
                2*binned_cosmic_variance(result,int(step),width),          # height
                color='r',
                alpha=0.1
            )
        )
        # print j, step
plt.savefig('AxiCLASS_bestfit_smallvsbestfitTheta.pdf', bbox_inches='tight')
# plt.savefig('spectra_nalp_TT_vs_TTTEEE.pdf', bbox_inches='tight')


# In[ ]:




# In[ ]:
