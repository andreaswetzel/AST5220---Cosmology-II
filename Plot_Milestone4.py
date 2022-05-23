import numpy as np
import matplotlib.pyplot as plt


CMB_data = np.loadtxt('cells.txt')
CMB_ell  = CMB_data[:, 0]
CMB_cell = CMB_data[:, 1]


TT_data = np.loadtxt('low_ell_TT_data.txt',skiprows=1)
TT_ell      = TT_data[:,0]
TT_cell     = TT_data[:,1]
TT_err_up   = TT_data[:,2]
TT_err_down = TT_data[:,3]


theta_data = np.loadtxt('theta_output.txt')
k_theta    = theta_data[:, 0]
theta_6  = theta_data[:, 1]
theta_100  = theta_data[:, 2]
theta_200  = theta_data[:, 3]
theta_1000= theta_data[:, 4]


theta_values = [theta_6,theta_100,theta_200,theta_1000]
theta_label = [r"$\ell=6$",r"$\ell=100$",r"$\ell=200$",r"$\ell=1000$"]

plt.plot(CMB_ell,CMB_cell)
plt.errorbar(TT_ell,TT_cell,np.array([TT_err_down,TT_err_up]),fmt='o',markersize=2)
plt.xscale('log')
plt.show()


for i in range(4):
    plt.plot(k_theta,theta_values[i],label=theta_label[i])

plt.legend()
plt.show()


for i in range(4):
    plt.plot(k_theta,theta_values[i]**2/k_theta,label=theta_label[i])

plt.legend()
plt.show()
