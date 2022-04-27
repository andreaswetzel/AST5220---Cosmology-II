import numpy as np
import matplotlib.pyplot as plt
import os


OmegaB0      = 0.05
OmegaCDM0    = 0.267
OmegaR0      = 5.50896e-05
OmegaLambda0 = 0.682945
a_r          = OmegaR0/(OmegaB0 + OmegaCDM0)                # Matter-radiation equality
a_d          = ((OmegaB0 + OmegaCDM0)/OmegaLambda0)**(1/3)  # Matter-dark energy equality

kvals = ['0.1', '0.01', '0.001']
col   = {'0.1': 'tab:blue', '0.01': 'tab:orange', '0.001': 'tab:green'}

#for k in kvals:
    #file    = 'perturbations_k' + k + '.txt'
data      = np.loadtxt("perturbations_k0.01.txt")
x         = data[:, 0]
theta0    = data[:, 1]
theta1    = data[:, 2]
theta2    = data[:, 3]
phi       = data[:, 4]
psi       = data[:, 5]
delta_b   = data[:, 6]
delta_cdm = data[:, 7]
v_b       = data[:, 8]
v_cdm     = data[:, 9]
ST        = data[:, 10]
a         = np.exp(x)
z         = 1/a - 1


plt.plot(x,ST)
plt.show()


plt.figure(1)
plt.title(r'$\delta_\gamma$')
plt.plot(x, 4*theta0)#=r'$k = $' + k + '/Mpc')
'''plt.xscale('log')'''
plt.legend()
plt.xlabel(r'$z$')
plt.show()

'''plt.figure(2)
plt.title(r'$v_\gamma$')
plt.plot(x, -3*theta1)#=r'$k = $' + k + '/Mpc')
#plt.xscale('log')
plt.legend()
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()'''

'''plt.figure(3)
plt.title(r'$\Theta_2$')
plt.plot(x, theta2)#=r'$k = $' + k + '/Mpc')
#plt.xscale('log')
plt.legend()
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()'''

'''plt.figure(4)
plt.title(r'$\Phi$')
plt.plot(x, phi)#=r'$k = $' + k + '/Mpc')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()'''

'''plt.figure(5)
plt.title(r'$\Psi$')
plt.plot(x, psi)#=r'$k = $' + k + '/Mpc')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()'''

plt.figure(6)
plt.title(r'$\delta_b, \delta_{CDM}$')
plt.plot(x, delta_cdm)#=r'$k = $' + k + '/Mpc', color=col[k])
plt.plot(x, np.abs(delta_b))#, ls='--', color=col[k])
plt.legend()
plt.yscale("log")
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()

'''plt.figure(7)
plt.title(r'$v_b, v_{CDM}$')
plt.log(x, v_cdm)#=r'$k = $' + k + '/Mpc', color=col[k])
plt.log(x, np.abs(v_b))#, ls='--', color=col[k])
plt.legend()
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()'''

'''plt.figure(8)
plt.title(r'$\Psi$')
plt.log(x, psi)#=r'$k = $' + k + '/Mpc', color=col[k])
plt.legend()
plt.xlabel(r'$z$')
plt.gca().invert_xaxis()
plt.tight_layout()
plt.show()'''

'''os.chdir('./Plots')
plt.figure(1)
#plt.savefig('theta0.pdf')
plt.figure(2)
#plt.savefig('theta1.pdf')
plt.figure(4)
#plt.savefig('phi.pdf')
plt.figure(5)
#plt.savefig('psi.pdf')
plt.figure(6)
#plt.savefig('delta.pdf')
plt.figure(7)
#plt.savefig('v.pdf')

plt.figure(8)
plt.title(r'$\Theta_0 + \Psi$')
plt.plot(z, theta0 - phi)
plt.xlabel(r'$z$')
plt.xscale('log')
plt.gca().invert_xaxis()
plt.tight_layout()

plt.show()
'''
