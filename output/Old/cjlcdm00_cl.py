import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/cjlcdm00_cl.dat', '/home/charlotte/class_perso-ScalarField/class_perso-ScalarField/output/CJ_SF60_cl.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['cjlcdm00_cl', 'CJ_SF60_cl']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = [2.0, 2000.0]
ax.loglog(curve[:, 0], abs(curve[:, 1]))

index, curve = 1, data[1]
y_axis = [u'TT']
tex_names = ['TT']
x_axis = 'l'
ylim = []
xlim = [2.0, 2000.0]
ax.loglog(curve[:, 0], abs(curve[:, 1]))

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('$\ell$', fontsize=16)
ax.set_xlim(xlim)
ax.set_ylim()
plt.show()