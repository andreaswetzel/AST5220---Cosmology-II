import numpy as np
import matplotlib.pyplot as plt


CMB_data = np.loadtxt('cells.txt')
CMB_ell  = CMB_data[:, 0]
CMB_cell = CMB_data[:, 1]


low_TT_data = np.loadtxt('low_ell_TT_data.txt',skiprows=1)
low_TT_ell      = low_TT_data[:,0]
low_TT_cell     = low_TT_data[:,1]
low_TT_err_up   = low_TT_data[:,2]
low_TT_err_down = low_TT_data[:,3]


high_TT_data = np.loadtxt('high_ell_TT_data.txt',skiprows=1)
high_TT_ell      = high_TT_data[:,0]
high_TT_cell     = high_TT_data[:,1]
high_TT_err_up   = high_TT_data[:,2]
high_TT_err_down = high_TT_data[:,3]



theta_data = np.loadtxt('theta_output.txt')
k_theta    = theta_data[:, 0]
kcH0       = theta_data[:, 1]
theta_6    = theta_data[:, 2]
theta_100  = theta_data[:, 3]
theta_200  = theta_data[:, 4]
theta_500  = theta_data[:, 5]
theta_1000 = theta_data[:, 6]
eta0       = 14253.3 # Mpc

matter_data = np.loadtxt("matter_output.txt")
kch         = matter_data[:,0]
MPS         = matter_data[:,1]
k_eq        = matter_data[:,2]

theta_values = [theta_6,theta_100,theta_200,theta_500,theta_1000]
theta_label  = [r"$\ell=6$",r"$\ell=100$",r"$\ell=200$",r"$\ell=500$",r"$\ell=1000$"]
ell_values   = [6,100,200,500,1000]

'''plt.plot(CMB_ell,CMB_cell)
plt.errorbar(low_TT_ell,low_TT_cell,np.array([low_TT_err_down,low_TT_err_up]),fmt='o',markersize=2,c='r')
plt.errorbar(high_TT_ell,high_TT_cell,np.array([high_TT_err_down,high_TT_err_up]),fmt='o',markersize=2,c='r')
plt.xscale('log')
plt.title("CMB Power Spectrum")
plt.xlabel(r'Multipole ${\ell}$')
plt.ylabel(r'$\ell(\ell+1)C_\ell$ $(\mu K)^2$')
plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M4/CMB_PS.png')
plt.show()
'''



for i in range(5):
    plt.plot(k_theta,theta_values[i],label=theta_label[i])

plt.axvline(x=6/eta0, ls='--', color='black',label=r'$k = \ell/\eta_0$')
plt.axvline(x=100/eta0, ls='--', color='black')
plt.axvline(x=200/eta0, ls='--', color='black')
plt.axvline(x=500/eta0, ls='--', color='black')
plt.axvline(x=1000/eta0, ls='--', color='black')
plt.xlabel(r'$k$' + ' [Mpc' + r'$^{-1}$' + ']')
plt.ylabel(r'$\Theta_\ell$')
plt.tight_layout()
plt.legend()
plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M4/theta_ell.png')
plt.show()

'''
for i in range(5):
    plt.plot(kcH0,theta_values[i]**2/k_theta,label=theta_label[i])

plt.legend()
plt.title('Spectrum integrand $\Theta_{\ell}(k)^2/k$')
plt.xlabel(r'$ck/H_0$')
plt.ylabel(r'$\Theta_\ell^2H_0/(kc)$')
plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M4/theta_l_2_k.png')
plt.show()

x = np.linspace(0.0202677,0.0202677,100)
y = np.linspace(300,32000,100)
plt.plot(kch,MPS)
plt.plot(x,y, c='k',label=r'$k_{eq}$')
plt.xscale('log')
plt.yscale('log')
plt.title('Total Matter Power-Spectrum')
plt.xlabel('Wavenumber $k[h/Mpc]$')
plt.ylabel(r'$P(k)[(Mpc/h)^3]$')
plt.legend()
plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M4/MPS.png')
plt.show()
'''
