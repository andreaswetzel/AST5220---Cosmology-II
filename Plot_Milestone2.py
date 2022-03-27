import numpy as np
import matplotlib.pyplot as plt

file = np.loadtxt('recombination.txt')



x     = file[:, 0]
Xe    = file[:, 1]
ne    = file[:, 2]
tau   = file[:, 3]
dtau  = file[:, 4]
ddtau = file[:, 5]
g     = file[:, 6]
dg    = file[:, 7]
ddg   = file[:, 8]
tau_at_1 = file[:, 9]

dec1 = np.linspace(tau_at_1[0],tau_at_1[0],100)
dec2 = np.linspace(0,1.5,100)

'''Plots of Xe'''
plt.plot (x, Xe, label = r"$X_e$", linewidth = 1)
plt.plot(dec1,dec2, '--', c='k', label=r'$\tau=1$')
plt.title (r"Evolution of the Electron Fraction")
plt.ylim (10**(-4), 10**1)
plt.xlim (-12, 0)
plt.yscale ("log")
plt.xlabel (r"x")
plt.ylabel(r"$X_e(x)$")
plt.legend()
plt.savefig('Electron_fraction.png')
plt.show ()
'''Plots of g and its derivatives'''
#print(np.argmin())
g_max = max(g)
dg_max = max(dg)
ddg_max = max(ddg)


plt.plot (x, g/g_max, label = r"$\tilde{g}(x)$", linewidth = 1)
plt.plot (x, dg/dg_max, label = r"$\frac{d\tilde{g}(x)}{dx}$", linewidth = 1)
plt.plot (x, ddg/ddg_max, label = r"$\frac{d^2\tilde{g}}{dx^2}(x)$", linewidth = 1)
plt.plot(dec1,dec2, '--', c='k', label=r'$\tau=1$')
plt.title (r"Evolution of the Visibility Function")
plt.ylim (-1.4, 1.1)
plt.xlim (-7.75, -5.7)
plt.xlabel (r"x")
plt.ylabel(r'$\tildeg(x)$')
plt.legend()
plt.savefig('Visibility_function')
plt.show ()


'''Plots of tau and its derivatives'''
plt.plot (x, tau, label = r"$\tau(x)$", linewidth = 1)
plt.plot (x, -1*dtau, label = r"$-\frac{d\tau(x)}{dx}$", linewidth = 1)
plt.plot (x, ddtau, label = r"$\frac{d^2\tau(x)}{dx^2}$", linewidth = 1)
plt.plot(dec1,dec2, '--', c='k', label=r'$\tau=1$')
plt.title (r"Evolution of the optical depth")
plt.ylim (1e-8,1e7)
plt.xlim (-12, 0)
plt.yscale ("log")
plt.xlabel (r"x")
plt.ylabel(r'$\tau$,$\frac{d\tau}{dx}$,$\frac{d^2\tau}{dx^2}$')
plt.legend()
plt.savefig('Optical Depth.png')
plt.show ()

# Plot of ne
'''plt.plot (x, ne, label = r"$n_e$", linewidth = 1)
plt.title (r"Evolution of Electron Number")
plt.xlabel (r"x")
plt.legend()
plt.savefig('Electron_number.png')
plt.show()'''
