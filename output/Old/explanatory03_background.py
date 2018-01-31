import matplotlib.pyplot as plt
import numpy as np
import itertools

files = ['/home/charlotte/Dropbox/JHU/class_perso-ScalarField/class_perso-ScalarField/output/explanatory03_background.dat']
data = []
for data_file in files:
    data.append(np.loadtxt(data_file))
roots = ['explanatory03_background']

fig, ax = plt.subplots()

index, curve = 0, data[0]
y_axis = [u'w_fld', u'w_scf']
tex_names = [u'(8\\pi G/3)w_fld', u'(8\\pi G/3)w_scf']
x_axis = 'z'
ylim = []
xlim = []
ax.semilogx(curve[:, 0], curve[:, 13])
ax.semilogx(curve[:, 0], curve[:, 18])

ax.legend([root+': '+elem for (root, elem) in
    itertools.product(roots, y_axis)], loc='best')

ax.set_xlabel('z', fontsize=16)
plt.show()