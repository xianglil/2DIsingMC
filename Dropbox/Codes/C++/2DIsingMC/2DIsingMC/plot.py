import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

with open("output.txt", 'r') as file:
	text = file.readlines()
data = np.array([i.split() for i in text][1:])
temperature = data[:,0]
energy = data[:,1]
magnetization = data[:,2]
specificHeat = data[:,3]
magneticSusceptibility = data[:,4]

x = temperature.astype(float)
fig, axes = plt.subplots(2, 2, sharex=True)

y = energy.astype(float)
ax = axes[0,0]
ax.plot(x, y, "ro")
ax.set_xlabel('T')
ax.set_ylabel('energy')

y = magnetization.astype(float)
ax = axes[0,1]
ax.plot(x, y, "bo")
ax.set_xlabel('T')
ax.set_ylabel('magnetization')

y = specificHeat.astype(float)
ax = axes[1,0]
ax.plot(x, y, "ro")
ax.set_xlabel('T')
ax.set_ylabel('specifit heat')

y = magneticSusceptibility.astype(float)
ax = axes[1,1]
ax.plot(x, y, "bo")
ax.set_xlabel('T')
ax.set_ylabel('magnetic susceptibility')

plt.tight_layout()
plt.show()
fig.savefig("Ising2DvsT.pdf", bbox_inches='tight')