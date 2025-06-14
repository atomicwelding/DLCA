import numpy as np
import matplotlib.pyplot as plt

import rcparams
from rcparams import named_colors as nc


phi20L,   P20L, _, _ = np.loadtxt('L20L.dat', comments="@", unpack=True)
phi20LL,  P20LL, _, _ = np.loadtxt('L20LL.dat', comments="@", unpack=True)
phi20LLL, P20LLL, _, _ = np.loadtxt('L20LLL.dat', comments="@", unpack=True)

plt.xlabel(r'$\phi_0$')
plt.ylabel(r'$P(\phi_0, L = 20)$')

plt.plot(phi20L, P20L, label='800',
         marker='o',
         linestyle='dotted',
         color=nc['blue'])

plt.plot(phi20LL, P20LL, label='4e3',
         marker='h',
         linestyle='dotted',
         color=nc['red'])

plt.plot(phi20LLL, P20LLL, label='8e3',
         marker='^',
         linestyle='dotted',
         color=nc['green'])


plt.legend()
plt.show()
