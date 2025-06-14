import numpy as np
import matplotlib.pyplot as plt

from rcparams import named_colors as nc
import rcparams


_, _, t40, _ = np.loadtxt('L40.dat', comments='@', unpack=True)
_, _, t50, _ = np.loadtxt('L50.dat', comments='@', unpack=True)
_, _, t60, _ = np.loadtxt('L60.dat', comments='@', unpack=True) 
_, _, t70, _ = np.loadtxt('L70.dat', comments='@', unpack=True) 

fig, axs = plt.subplots(2, 2, figsize=(10, 8))

axs[0,0].hist(t40,
            color=nc['black'],
            label=f"tgel = {np.mean(t40):.2f}")
axs[0,0].set_title("L = 40", fontsize=18)
axs[0,0].set_xlabel(r"$t$")
axs[0,0].legend()

axs[0,1].hist(t50,
              color=nc['black'],
              label=f"tgel = {np.mean(t50):.2f}")
axs[0,1].set_title("L = 50", fontsize=18)
axs[0,1].set_xlabel(r"$t$")
axs[0,1].legend()

axs[1,0].hist(t60,
              color=nc['black'],
              label=f"tgel = {np.mean(t60):.2f}")
axs[1,0].set_title("L = 60", fontsize=18)
axs[1,0].set_xlabel(r"$t$")
axs[1,0].legend()

axs[1,1].hist(t70,
              color=nc['black'],
              label=f"tgel = {np.mean(t70):.2f}")
axs[1,1].set_title("L = 70", fontsize=18)
axs[1,1].set_xlabel(r"$t$")
axs[1,1].legend()

fig.tight_layout()
plt.show()
