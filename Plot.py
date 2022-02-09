import numpy as np
import matplotlib.pyplot as plt
from math import exp

infile = open("cosmology.txt", "r")

x = []
eta_of_x = []
Hp_of_x = []
dHpdx_of_x = []
Total_Matter = []
get_OmegaCDM = []
get_OmegaLambda = []
get_OmegaR = []
get_OmegaNu = []
get_OmegaK = []

for line in infile:
    values = line.split()
    x.append(float(values[0]))
    eta_of_x.append(float(values[1]))
    Hp_of_x.append(float(values[2]))
    dHpdx_of_x.append(float(values[3]))
    Total_Matter.append(float(values[4]))
    get_OmegaCDM.append(float(values[5]))
    get_OmegaLambda.append(float(values[6]))
    get_OmegaR.append(float(values[7]))
    get_OmegaNu.append(float(values[8]))
    get_OmegaK.append(float(values[9]))

a = np.array(x)

plt.plot(a,eta_of_x, label="$\eta (x)$")
plt.xlabel("Scalefactor a")
plt.ylabel("$\eta (x)$ [Mpc]")
plt.title("Evolution of Conformal Time")
# plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()


plt.plot(a,Total_Matter, label="$\Omega_B$")
plt.xlabel("Scalefactor a")
plt.ylabel("$\eta (x)$ [Mpc]")
plt.title("Evolution of Conformal Time")
# plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.show()
