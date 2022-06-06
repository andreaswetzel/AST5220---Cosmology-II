import numpy as np
import matplotlib.pyplot as plt
from math import exp

infile = open("cosmology.txt", "r")
infile_SN = open("supernovadata.txt", "r")




x = []
eta_of_x       = []
Hp_of_x         = []
dHpdx_of_x      = []
get_OmegaB      = []
get_OmegaCDM    = []
get_OmegaLambda = []
get_OmegaR      = []
get_OmegaNu     = []
get_OmegaK      = []
t_of_x        = []
ddHpddx_of_x = []
H_of_x = []
H0 = []
z = []
dA = []
dL = []
chi = []
t0 = []
time_today = []
eta_Hp_c = []

for line in infile: #Filling in the empty lists
    values = line.split()
    x.append(float(values[0]))
    eta_of_x.append(float(values[1]))
    Hp_of_x.append(float(values[2]))
    dHpdx_of_x.append(float(values[3]))
    get_OmegaB.append(float(values[4]))
    get_OmegaCDM.append(float(values[5]))
    get_OmegaLambda.append(float(values[6]))
    get_OmegaR.append(float(values[7]))
    get_OmegaNu.append(float(values[8]))
    get_OmegaK.append(float(values[9]))
    t_of_x.append(float(values[10]))
    ddHpddx_of_x.append(float(values[11]))
    H_of_x.append(float(values[12]))
    H0.append(float(values[13]))
    z.append(float(values[14]))
    dA.append(float(values[15]))
    dL.append(float(values[16]))
    chi.append(float(values[17]))
    t0.append(float(values[18]))
    time_today.append(float(values[19]))
    eta_Hp_c.append(float(values[20]))


z_SN = []
dL_SN = []
Error_SN = []


for line in infile_SN: #Filling lists from SN data
    values = line.split()
    z_SN.append(float(values[0]))
    dL_SN.append(float(values[1]))
    Error_SN.append(float(values[2]))


x = np.array(x) #Making x to an array
a = np.exp(x)

x = np.log(a)
idx = np.argmin(abs(x)) #Todays x

a_rm = (get_OmegaR[idx]+get_OmegaNu[idx])/(get_OmegaB[idx]+get_OmegaCDM[idx]) #Matter-radiation equality
a_mlam = ((get_OmegaB[idx]+get_OmegaCDM[idx]) / get_OmegaLambda[idx])**(1/3) #Matter-Dark energy equality
a_acc = ((get_OmegaB[idx]+get_OmegaCDM[idx]) / (2*get_OmegaLambda[idx]))**(1/3) #Accelerating universe


Hx_H0 = [] #H/H_0
for i in range(len(H0)):
    Hx_H0.append(float(H_of_x[i]/H0[i]))

Hp_H0 = [] #H/H_0
for i in range(len(H0)):
    Hp_H0.append(float(Hp_of_x[i]/H0[i]))

plt.plot(a,Hx_H0)
plt.ylabel(r"$\frac{H}{H_0}$")
plt.xlabel('Scalefactor a')
plt.title('Hubble Parameter')
plt.vlines(a_rm,-1e-1,1e10,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-1e-1,1e10,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-1e-1,1e10,linestyles ="dashed",label='Accelerating universe',color='r')
plt.legend()
plt.yscale("log")
plt.xscale("log")
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M1/Hubble Parameter.png')
plt.show()

plt.plot(a,Hp_H0)
plt.ylabel(r"$\frac{\mathcal{H}}{H_0}$")
plt.xlabel('Scalefactor a')
plt.title('Scaled Hubble Parameter')
plt.vlines(a_rm,1e0,1e5,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,1e-1,1e5,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,1e-1,1e5,linestyles ="dashed",label='Accelerating universe',color='r')
plt.legend()
plt.yscale("log")
plt.xscale("log")
#plt.savefig('/Users/andreaswetzel/Documents/Master_Astro/2_semester/AST5220/AST5220-Cosmology/Figures_M1/Scaled Hubble Parameter.png')
plt.show()





dHp_Hp = [] #Derivative of Hp
for i in range(len(Hp_of_x)):
    dHp_Hp.append(float(dHpdx_of_x[i]/Hp_of_x[i]))

ddHp_Hp = [] # Second Derivative of Hp
for i in range(len(Hp_of_x)):
    ddHp_Hp.append(float(ddHpddx_of_x[i]/Hp_of_x[i]))




plt.plot(a,eta_of_x)
plt.xlabel("Scalefactor a")
plt.ylabel("$\eta $ [Mpc]")
plt.title("Evolution of Conformal Time")
plt.xscale("log")
plt.yscale("log")
plt.vlines(a_rm,-1.5,2e4,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-1.5,2e4,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-1.5,2e4,linestyles ="dashed",label='Accelerating universe',color='r')
plt.legend()
plt.savefig('Evolution of Conformal Time.png')
plt.show()



plt.plot(a,dHp_Hp, label=r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$")
plt.plot(a,ddHp_Hp, label=r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")
plt.xlabel("Scalefactor a")
plt.ylabel(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx},\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$")
plt.title(r"Evolution of Hubble Factors (derivatives of $\mathcal{H}$)")
plt.xscale("log")

plt.vlines(a_rm,-1.5,1.5,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-1.5,1.5,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-1.5,1.5,linestyles ="dashed",label='Accelerating universe',color='r')
plt.legend()
plt.savefig('der_Hp.png')
plt.show()

OB = np.array(get_OmegaB) #From list to array
OCDM = np.array(get_OmegaCDM) #From list to array

plt.plot(a,OB+OCDM, label="Total matter") #Feil  - BAryonic
plt.plot(a,get_OmegaLambda, label="Dark energy") #Dark energy
plt.plot(a,get_OmegaR, label="Radiation") #radiation
plt.xlabel("Scalefactor a")
plt.ylabel(r'$\Omega$')
plt.title('Evolution of the Density Parameters')
plt.vlines(a_rm,-0.1,1.1,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-0.1,1.1,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-0.1,1.1,linestyles ="dashed",label='Accelerating universe',color='r')
plt.xscale("log")
plt.legend()
plt.savefig('Omega')
plt.show()


#Hx/H0 plot if you want to see
'''plt.plot(a,Hx_H0)
plt.yscale("log")
plt.xscale("log")
plt.show()'''


plt.plot(a,eta_Hp_c)
plt.yscale("log")
plt.xscale("log")
plt.ylabel(r'$\frac{\eta(x)\mathcal{H}(x)}{c}$')
plt.xlabel('Scalefactor $a$')
plt.vlines(a_rm,-1.5,2e1,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-1.5,2e1,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-1.5,2e1,linestyles ="dashed",label='Accelerating universe',color='r')
plt.legend()
plt.savefig('eta_Hp_c.png')
plt.show()


t_of_x = np.array(t_of_x)
plt.plot(a,t_of_x/(1e9*365*24*60**2))
plt.yscale("log")
plt.xscale("log")
plt.title('Age of the Universe')
plt.ylabel('$t(x)$')
plt.xlabel('Scalarfactor $a$')
plt.vlines(a_rm,-1.5,2e4,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-1.5,2e4,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-1.5,2e4,linestyles ="dashed",label='Accelerating universe',color='r')
plt.legend()
plt.savefig('t_of_x')
plt.show()


#Plot for distance Measure, was not a part of the assignment, more like a test
'''H0_z = []

for i in range(len(H0)):
    H0_z.append(float(H0[i]*z[i]))

H0_z_ = []#[element * (3e8*H0) for element in H0_z]

for i in range(len(H0_z)):
    H0_z_.append(float(H0_z[i]*3e8/H0[i]))



H0_z__ = [element * (1e-5) for element in H0_z_]



dL_ = []
dA_ = []
H0_z_ = []
chi_ = []
for i in range(len(dL)): #Luminosity-, Angular distance, chi, and H0*z
    dL_.append(float(dL[i]*1000))
    dA_.append(float(dA[i]*1000))
    H0_z_.append(float(H0_z[i]*1000))
    chi_.append(float(chi[i]*1000))


plt.plot(z,dL, label='Luminosity distance')
plt.plot(z,dA, label='Angular diameter distance')
plt.plot(z,H0_z__, label="Naive Hubble Parameter")
plt.plot(z,chi, label='Comoving distance')
plt.vlines(a_rm,-1.5,2e8,linestyles ="dashed",label='Matter-Radiation equality',color='b')
plt.vlines(a_mlam,-1.5,2e8,linestyles ="dashed",label='Matter-Dark energy equality',color='g')
plt.vlines(a_acc,-1.5,2e8,linestyles ="dashed",label='Accelerating universe',color='r')
plt.xlabel('Redshift z')
plt.ylabel('Distance [Mpc]')
plt.title('Evolution of Cosmological Distance Measure')
plt.yscale("log")
plt.xscale("log")
plt.legend()
plt.savefig('Cos Dis.png')
plt.show()'''


TP = [] #Theoretical predictions

for i in range(len(dL)):
    TP.append(float(dL[i]/1000))
plt.plot(z,TP,label="Therotical predictions")
plt.plot(z_SN,dL_SN, "o", label="Supernova data")
plt.legend()
plt.yscale("log")
plt.xscale("log")
plt.xlabel('Redshift z')
plt.ylabel('Distance [Gpc]')
plt.title('Luminosity distance vs Supernova Data')
plt.xlim(0.5e-2,1.5e0)
plt.ylim(1e-2,1.5e1)
plt.savefig('LD_vs_SN')
plt.show()
