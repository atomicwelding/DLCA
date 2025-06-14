import numpy as np
import matplotlib.pyplot as plt

from rcparams import named_colors as nc
import rcparams


phi05, P05, _, _ = np.loadtxt('L05.dat', comments='@', unpack=True)
phi10, P10, _, _ = np.loadtxt('L10.dat', comments='@', unpack=True)
phi20, P20, _, _ = np.loadtxt('L20.dat', comments='@', unpack=True)
phi40, P40, _, _ = np.loadtxt('L40.dat', comments='@', unpack=True)

plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$P(\phi_0,L)$')


plt.plot(phi40, P40,
         color = nc['red'],
         linestyle='dotted',
         marker='o',
         markersize=9,
         label='L = 40')

plt.plot(phi20, P20,
          color = nc['blue'],
         linestyle='dotted',
         marker='h',
         markersize=9,
         label='L = 20')

plt.plot(phi10, P10,
         color = nc['green'],
         linestyle='dotted',
         marker='v',
         markersize=9,
         label='L = 10')

plt.plot(phi05, P05,
         color = nc['black'],
         linestyle='dotted',
         marker='^',
         markersize=9,
         label='L = 05')



plt.legend()
plt.show()


