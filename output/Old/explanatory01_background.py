import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/charlotte/Dropbox/JHU/class_perso-ScalarField/class_perso-ScalarField/output/explanatory01_background.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['explanatory01_background']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'rho_cdm', u'rho_g', u'rho_b', u'rho_lambda', u'rho_ur']
tex_names = [u'(8\\pi G/3)rho_g', u'(8\\pi G/3)rho_b', u'(8\\pi G/3)rho_cdm', u'(8\\pi G/3)rho_lambda', u'(8\\pi G/3)rho_ur']
x_axis = 'z'
ylim = []
xlim = []
ax.loglog(curve[:, 0], abs(curve[:, 10]))
ax.loglog(curve[:, 0], abs(curve[:, 8]))
ax.loglog(curve[:, 0], abs(curve[:, 9]))
ax.loglog(curve[:, 0], abs(curve[:, 11]))
ax.loglog(curve[:, 0], abs(curve[:, 12]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()