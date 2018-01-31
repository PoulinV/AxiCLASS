import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_28_no_shoot_fluidbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_28_no_shoot_kgbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_24_no_shoot_fluidbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_24_no_shoot_kgbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_20_no_shoot_fluidbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_20_no_shoot_kgbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_31_no_shoot_fluidbackground.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cubed_m_31_no_shoot_kgbackground.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['cubed_m_28_no_shoot_fluidbackground', 'cubed_m_28_no_shoot_kgbackground', 'cubed_m_24_no_shoot_fluidbackground', 'cubed_m_24_no_shoot_kgbackground', 'cubed_m_20_no_shoot_fluidbackground', 'cubed_m_20_no_shoot_kgbackground', 'cubed_m_31_no_shoot_fluidbackground', 'cubed_m_31_no_shoot_kgbackground']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 1, data[1]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 2, data[2]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 3, data[3]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 4, data[4]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 5, data[5]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 6, data[6]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

index, curve = 7, data[7]
y_axis = [u'rho_scf']
tex_names = [u'(8\\pi G/3)rho_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 14]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()