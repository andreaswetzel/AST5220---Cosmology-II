import numpy as np
import matplotlib.pyplot as plt

data1       = np.loadtxt("perturbations_k0.1.txt")
x           = data1[:, 0]
theta_1     = data1[:, 1]
theta1_1    = data1[:, 2]
theta2_1    = data1[:, 3]
phi_1       = data1[:, 4]
psi_1       = data1[:, 5]
delta_b_1   = data1[:, 6]
delta_cdm_1 = data1[:, 7]
v_b_1       = data1[:, 8]
v_cdm_1     = data1[:, 9]
ST_1        = data1[:, 10]


data01       = np.loadtxt("perturbations_k0.01.txt")
theta_01     = data01[:, 1]
theta1_01    = data01[:, 2]
theta2_01    = data01[:, 3]
phi_01       = data01[:, 4]
psi_01       = data01[:, 5]
delta_b_01   = data01[:, 6]
delta_cdm_01 = data01[:, 7]
v_b_01       = data01[:, 8]
v_cdm_01     = data01[:, 9]
ST_01        = data01[:, 10]


data001       = np.loadtxt("perturbations_k0.001.txt")
theta_001     = data001[:, 1]
theta1_001    = data001[:, 2]
theta2_001    = data001[:, 3]
phi_001       = data001[:, 4]
psi_001       = data001[:, 5]
delta_b_001   = data001[:, 6]
delta_cdm_001 = data001[:, 7]
v_b_001       = data001[:, 8]
v_cdm_001     = data001[:, 9]
ST_001        = data001[:, 10]


plt.plot(x,ST_1)
plt.show()
'''plt.plot(x, delta_cdm_1,label=r'$\delta_{CDM}$ - k=0.1[Mpc]')
plt.plot(x, delta_cdm_01,label=r'$\delta_{CDM}$ - k=0.01[Mpc]')
plt.plot(x, delta_cdm_001,label=r'$\delta_{CDM}$ - k=0.001[Mpc]')
plt.plot(x, np.abs(delta_b_1),'--',label=r'$\delta_b$ - k=0.1[Mpc]')
plt.plot(x, np.abs(delta_b_01),'--',label=r'$\delta_b$ - k=0.01[Mpc]')
plt.plot(x, np.abs(delta_b_001),'--',label=r'$\delta_b$ - k=0.001[Mpc]')
plt.yscale('log')
plt.legend()
plt.title(r'$\delta_{CDM},\delta_b$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\delta_{CDM}(x,k),\delta_b(x,k)$')
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M3/delta_.png')
plt.show()


plt.plot(x, v_cdm_1,label=r'$v_{CDM}$ - k=0.1[Mpc]')
plt.plot(x, v_cdm_01,label=r'$v_{CDM}$ - k=0.01[Mpc]')
plt.plot(x, v_cdm_001,label=r'$v_{CDM}$ - k=0.001[Mpc]')
plt.plot(x, np.abs(v_b_1),'--',label=r'$v_b$ - k=0.1[Mpc]')
plt.plot(x, np.abs(v_b_01),'--',label=r'$v_b$ - k=0.01[Mpc]')
plt.plot(x, np.abs(v_b_001),'--',label=r'$v_b$ - k=0.001[Mpc]')
plt.yscale('log')
plt.legend()
plt.title(r'$v_b, v_{CDM}$')
plt.xlabel(r'$x$')
plt.ylabel(r'$v_b(x,k), v_{CDM}(x,k)$')
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M3/v_.png')
plt.show()


plt.plot(x, 4*theta_1,label='k=0.1[Mpc]')
plt.plot(x, 4*theta_01,label='k=0.01[Mpc]')
plt.plot(x, 4*theta_001,label='k=0.001[Mpc]')
plt.legend()
plt.title(r'$\delta_\gamma$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\delta_\gamma(x,k)$')
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M3/delta_gamma.png')
plt.show()


plt.plot(x, -3*theta1_1,label='k=0.1[Mpc]')
plt.plot(x, -3*theta1_01,label='k=0.01[Mpc]')
plt.plot(x, -3*theta1_001,label='k=0.001[Mpc]')
plt.legend()
plt.title(r'$v_\gamma$')
plt.xlabel(r'$x$')
plt.ylabel(r'$v_\gamma(x,k)$')
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M3/v_gamma.png')
plt.show()


plt.plot(x, phi_1,label='k=0.1[Mpc]')
plt.plot(x, phi_01,label='k=0.01[Mpc]')
plt.plot(x, phi_001,label='k=0.001[Mpc]')
plt.legend()
plt.title(r'$\Phi$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\Phi(x,k)$')
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M3/phi.png')
plt.show()


plt.plot(x, psi_1,label='k=0.1[Mpc]')
plt.plot(x, psi_01,label='k=0.01[Mpc]')
plt.plot(x, psi_001,label='k=0.001[Mpc]')
plt.legend()
plt.title(r'$\Psi$')
plt.xlabel(r'$x$')
plt.ylabel(r'$\Psi(x,k)$')
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M3/psi.png')
plt.show()
'''
