
# coding: utf-8

# # Load Libs

# In[1]:

### Load Python Lib##### 
########
import numpy as np
from math import *
import matplotlib.pyplot as plt
from random import gauss, randint
from matplotlib import rc
from matplotlib.ticker import FixedLocator, MultipleLocator, FormatStrFormatter

#very import to enable this so the plots can be showed in the page
get_ipython().magic(u'matplotlib inline')

# Use LaTeX font.
plt.rc('text', usetex=True)
plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman'],'size':20})

import matplotlib.font_manager as font_manager
font_prop = font_manager.FontProperties(size=12)


# # Load MC Data

# ## e-p Data

# In[3]:

### Loading the theoretical calculations
filename = 'xz_bins/prot_output_Q2_1.dat'
lines = open(filename,'r').readlines()

N = len(lines)-1
W = np.zeros(N, dtype=float)
Q2 = np.zeros(N, dtype=float)
x = np.zeros(N, dtype=float)
z = np.zeros(N, dtype=float)
pt = np.zeros(N, dtype=float)
asy = np.zeros(N, dtype=float)     
xs_inc = np.zeros(N, dtype=float)  
xs_pip = np.zeros(N, dtype=float)  
xs_pim = np.zeros(N, dtype=float) 
N_inc = np.zeros(N, dtype=float)   
N_pip = np.zeros(N, dtype=float)   
N_pim = np.zeros(N, dtype=float)   
mulp_pip = np.zeros(N, dtype=float)   
mulp_pim = np.zeros(N, dtype=float)   
Asym = np.zeros(N, dtype=float)   
Astat = np.zeros(N, dtype=float)   
Ax = np.zeros(N, dtype=float)   
mulp_err_pip = np.zeros(N, dtype=float)  
mulp_err_pim = np.zeros(N, dtype=float)   
Q2_avg = 0.0
#bin Q2          x          W          z         pt     xs_inc     xs_pip     xs_pim
#N_inc        N_pip        N_pim     mulp_pip     mulp_pim         Asym        Astat

## Read-In 
for i in range(0, N):
    values=lines[i+1].split()
    values = np.array(values,dtype=float)
    Q2[i]=(values[1])
    x[i]=(values[2])
    W[i]=(values[3])
    z[i]=(values[4])
    pt[i]=(values[5])
    xs_inc[i]=(values[6])
    xs_pip[i]=(values[7])
    xs_pim[i]=(values[8])
    N_inc[i]=(values[9])
    N_pip[i]=(values[10])
    N_pim[i]=(values[11])
    mulp_pip[i]=(values[12])
    mulp_pim[i]=(values[13])
    Asym[i]=(values[14])
    #Asym[i] = (mulp_pip[i]-mulp_pim[i]) #/(mulp_pip[i]+mulp_pim[i])
    Astat[i]=(values[15])
    Ax[i]=(values[16])

    Q2_avg += Q2[i]
Q2_avg /= (N-1)


# ## e-C12 Data

# In[4]:

### Loading the theoretical calculations
filename = 'xz_bins/c12_output_Q2_1.dat'
lines = open(filename,'r').readlines()

N = len(lines)-1
W_c12 = np.zeros(N, dtype=float)
Q2_c12 = np.zeros(N, dtype=float)
x_c12 = np.zeros(N, dtype=float)
z_c12 = np.zeros(N, dtype=float)
pt_c12 = np.zeros(N, dtype=float)
asy_c12 = np.zeros(N, dtype=float)     
xs_inc_c12 = np.zeros(N, dtype=float)  
xs_pip_c12 = np.zeros(N, dtype=float)  
xs_pim_c12 = np.zeros(N, dtype=float) 
N_inc_c12 = np.zeros(N, dtype=float)   
N_pip_c12 = np.zeros(N, dtype=float)   
N_pim_c12 = np.zeros(N, dtype=float)   
mulp_pip_c12 = np.zeros(N, dtype=float)   
mulp_pim_c12 = np.zeros(N, dtype=float)   
Asym_c12 = np.zeros(N, dtype=float)   
Astat_c12 = np.zeros(N, dtype=float)   
Ax_c12 = np.zeros(N, dtype=float)   
mulp_err_pip_c12 = np.zeros(N, dtype=float)  
mulp_err_pim_c12 = np.zeros(N, dtype=float)   
Q2_avg_c12 = 0.0
#bin Q2          x          W          z         pt     xs_inc     xs_pip     xs_pim
#N_inc        N_pip        N_pim     mulp_pip     mulp_pim         Asym        Astat

## Read-In 
for i in range(0, N):
    values=lines[i+1].split()
    values = np.array(values,dtype=float)
    Q2_c12[i]=(values[1])
    x_c12[i]=(values[2])
    W_c12[i]=(values[3])
    z_c12[i]=(values[4])
    pt_c12[i]=(values[5])
    xs_inc_c12[i]=(values[6])
    xs_pip_c12[i]=(values[7])
    xs_pim_c12[i]=(values[8])
    N_inc_c12[i]=(values[9])
    N_pip_c12[i]=(values[10])
    N_pim_c12[i]=(values[11])
    mulp_pip_c12[i]=(values[12])
    mulp_pim_c12[i]=(values[13])
    Asym_c12[i]=(values[14])
    #Asym[i] = (mulp_pip[i]-mulp_pim[i]) #/(mulp_pip[i]+mulp_pim[i])
    Astat_c12[i]=(values[15])
    Ax_c12[i]=(values[16])

    Q2_avg_c12 += Q2_c12[i]
Q2_avg_c12 /= (N-1)


# # HERMES Data

# In[5]:

### Loading the theoretical calculations
filename = 'HERMES_zx_pip.dat'
line1 = open(filename,'r').readlines()
filename = 'HERMES_zx_pim.dat'
line2 = open(filename,'r').readlines()

N = 5
Asym_H = np.zeros(N, dtype=float)   
Asym_Herr = np.zeros(N, dtype=float)   

Q2_p1 = np.zeros(N, dtype=float)
x_p1 = np.zeros(N, dtype=float)
z_p1 = np.zeros(N, dtype=float)
pt_p1 = np.zeros(N, dtype=float)
mulp_p1 = np.zeros(N, dtype=float)   
stat_p1 = np.zeros(N, dtype=float)  
syst_p1 = np.zeros(N, dtype=float) 

Q2_p2 = np.zeros(N, dtype=float)
x_p2 = np.zeros(N, dtype=float)
z_p2 = np.zeros(N, dtype=float)
pt_p2 = np.zeros(N, dtype=float)
mulp_p2 = np.zeros(N, dtype=float)   
stat_p2 = np.zeros(N, dtype=float)  
syst_p2 = np.zeros(N, dtype=float)  

Q2_pip = np.zeros(N, dtype=float)
x_pip = np.zeros(N, dtype=float)
z_pip = np.zeros(N, dtype=float)
pt_pip = np.zeros(N, dtype=float)
Mul_pip = np.zeros(N, dtype=float)   
stat_pip = np.zeros(N, dtype=float)  
syst_pip = np.zeros(N, dtype=float)  
err_pip = np.zeros(N, dtype=float)  

## Read-In 
for i in range(0, N):
    values=line1[i+2].split()
    values = np.array(values,dtype=float)
    mulp_p1[i]=(values[0])
    stat_p1[i]=(values[1])
    stat_p1[i]=(values[2])
    Q2_p1[i]=(values[3])
    x_p1[i]=(values[4])
    z_p1[i]=(values[5])
    pt_p1[i]=(values[6])

    values=line1[i+8].split()
    values = np.array(values,dtype=float)
    mulp_p2[i]=(values[0])
    stat_p2[i]=(values[1])
    stat_p2[i]=(values[2])
    Q2_p2[i]=(values[3])
    x_p2[i]=(values[4])
    z_p2[i]=(values[5])
    pt_p2[i]=(values[6])
    
    Mul_pip[i] = (mulp_p1[i] + mulp_p2[i])/2.0
    syst_pip[i] = (syst_p1[i] + syst_p2[i])/2.0
    stat_pip[i] = (stat_p1[i] + stat_p2[i])/2.0
    err_pip[i] = sqrt(syst_pip[i]**2 + stat_pip[i]**2)
    
    Q2_pip[i] = (Q2_p1[i] + Q2_p2[i])/2.0
    x_pip[i] = (x_p1[i] + x_p2[i])/2.0
    z_pip[i] = (z_p1[i] + z_p2[i])/2.0
    pt_pip[i] = (pt_p1[i] + pt_p2[i])/2.0

Q2_pim = np.zeros(N, dtype=float)
x_pim = np.zeros(N, dtype=float)
z_pim = np.zeros(N, dtype=float)
pt_pim = np.zeros(N, dtype=float)
Mul_pim = np.zeros(N, dtype=float)   
stat_pim = np.zeros(N, dtype=float)  
syst_pim = np.zeros(N, dtype=float)  
err_pim = np.zeros(N, dtype=float)  
err_total = np.zeros(N, dtype=float)  

for i in range(0, N):
    values=line2[i+2].split()
    values = np.array(values,dtype=float)
    mulp_p1[i]=(values[0])
    stat_p1[i]=(values[1])
    stat_p1[i]=(values[2])
    Q2_p1[i]=(values[3])
    x_p1[i]=(values[4])
    z_p1[i]=(values[5])
    pt_p1[i]=(values[6])

    values=line2[i+8].split()
    values = np.array(values,dtype=float)
    mulp_p2[i]=(values[0])
    stat_p2[i]=(values[1])
    stat_p2[i]=(values[2])
    Q2_p2[i]=(values[3])
    x_p2[i]=(values[4])
    z_p2[i]=(values[5])
    pt_p2[i]=(values[6])
    
    Mul_pim[i] = (mulp_p1[i] + mulp_p2[i])/2.0
    syst_pim[i] = (syst_p1[i] + syst_p2[i])/2.0
    stat_pim[i] = (stat_p1[i] + stat_p2[i])/2.0
    err_pim[i] = sqrt(syst_pim[i]**2 + stat_pim[i]**2)

    Q2_pim[i] = (Q2_p1[i] + Q2_p2[i])/2.0
    x_pim[i] = (x_p1[i] + x_p2[i])/2.0
    z_pim[i] = (z_p1[i] + z_p2[i])/2.0
    pt_pim[i] = (pt_p1[i] + pt_p2[i])/2.0
    
    Asym_H[i] = (Mul_pip[i]-Mul_pim[i])/(Mul_pip[i]+Mul_pim[i])
    Asym_Herr[i] = Asym_H[i]*sqrt( (err_pip[i]/Mul_pip[i])**2 + (err_pim[i]/Mul_pim[i])**2)
    err_total[i] = sqrt( err_pip[i]**2 + err_pim[i]**2)
    


# # Renormalize Mulplicity

# In[6]:

sum_pip = mulp_pip.sum()
sum_pim = mulp_pim.sum()
for i in range(0, len(z)):
    mulp_pip[i] /= sum_pip
    mulp_pim[i] /= sum_pim
    mulp_err_pip[i] = mulp_pip[i] * Astat[i]
    mulp_err_pim[i] = mulp_pim[i] * Astat[i]
    Asym[i]=( mulp_pip[i]- mulp_pim[i])/( mulp_pip[i] + mulp_pim[i])

sum_pip_c12 = mulp_pip_c12.sum()
sum_pim_c12 = mulp_pim_c12.sum()
for i in range(0, len(z_c12)):
    mulp_pip_c12[i] /= sum_pip_c12
    mulp_pim_c12[i] /= sum_pim_c12
    mulp_err_pip_c12[i] = mulp_pip_c12[i] * Astat_c12[i]
    mulp_err_pim_c12[i] = mulp_pim_c12[i] * Astat_c12[i]
    Asym_c12[i]=( mulp_pip_c12[i]- mulp_pim_c12[i])/( mulp_pip_c12[i] + mulp_pim_c12[i])
    
sum_pip = Mul_pip.sum()
sum_pim = Mul_pim.sum()
for i in range(0, len(z_pip)):
    Mul_pip[i] /= sum_pip
    Mul_pim[i] /= sum_pim
    err_pip[i] /= sum_pim 
    err_pim[i] /= sum_pim
    
    Asym_H[i] = (Mul_pip[i]-Mul_pim[i])/(Mul_pip[i]+Mul_pim[i])
    Asym_Herr[i] = sqrt( (err_pip[i]/Mul_pip[i])**2 + (err_pim[i]/Mul_pim[i])**2)
    err_total[i] = sqrt( err_pip[i]**2 + err_pim[i]**2)
    
x_Havg = x_pip.sum()/len(x_pip)
Q2_Havg = Q2_pip.sum()/len(Q2_pip)


# # Plotting

# In[7]:

## Plot GEn##########{{{
f1, ax = plt.subplots(3, 1, sharex=False, figsize=(12,12))
f1.subplots_adjust(bottom=0.08, top=0.96, hspace=0.5)
#plt.grid()

#######  Subplot for Pi+ Mulplicity
axis=ax[0]

axis.errorbar(z_pip, Mul_pip, yerr=err_pip,fmt='D', color='red',  label='HERMES')

axis.errorbar(z, mulp_pip, yerr=Astat ,fmt='o', color='blue',  label="$M^{p(e,e' \pi^{+})X}(z)$")
axis.errorbar(z_c12, mulp_pip_c12, yerr=Astat ,fmt='^', color='g',  label="$M^{^{12}C(e,e' \pi^{+})X}(z)$")

axis.set_title(r'$\pi^{+}$ Mulplicity (normalized)')
axis.set_ylabel(r'$M^{\pi^{+}}(z)$')
axis.set_xlabel(r'$z~(GeV/c)$')
axis.set_xlim(0.15,0.95)
axis.set_ylim(-0.1,0.6)
axis.grid()

#######  Subplot for Pi- Mulplicity
axis=ax[1]

axis.errorbar(z_pim, Mul_pim, yerr=err_pim,fmt='D', color='red',  label='ep-HERMES')

axis.errorbar(z, mulp_pim, yerr=Astat ,fmt='o', color='blue',  label='e-p EIC')
axis.errorbar(z_c12, mulp_pim_c12, yerr=Astat_c12 ,fmt='^', color='g',  label=r"e-$^{12}$C~EIC")

axis.set_title(r'$\pi^{-}$ Mulplicity (normalized)')
axis.set_ylabel(r'$M^{\pi^{-}}(z)$')
axis.set_xlabel(r'$z~(GeV/c)$')
axis.set_xlim(0.15,0.95)
axis.set_ylim(-0.1,0.6)
axis.grid()

#######  Subplot for Pi+ Mulplicity
axis=ax[2]

axis.errorbar(z_pip, Mul_pip-Mul_pim, yerr=err_total ,fmt='D', color='red',  label='HERMES')

axis.errorbar(z, mulp_pip-mulp_pim, yerr=Astat ,fmt='o', color='blue',  label='$Different of Mulplicity(\pi^{+}-\pi^{-})$')
axis.errorbar(z_c12, mulp_pip_c12-mulp_pim_c12, yerr=Astat_c12 ,fmt='^', color='g',  label='$Different of Mulplicity(\pi^{+}-\pi^{-})$')

axis.set_ylabel(r'$M^{\pi^{+}}(z)-M^{\pi^{-}}(z)$')
axis.set_title(r'Different of Mulplicity $(\pi^{+}-\pi^{-})$ (normalized)')
axis.set_xlabel(r'$z~(GeV/c)$')
axis.set_xlim(0.15,0.95)
axis.set_ylim(-0.08,0.035)
axis.grid()

ax[0].text(0.4, 0.45, r'$0.08<x_{B}<0.12, <Q^2>=%4.2f ~GeV^2$'%Q2_avg, color = 'b')
ax[0].text(0.4, 0.35, r'$<x_{B}>=%4.3f, <Q^2>=%4.2f ~GeV^2$'%(x_Havg,Q2_Havg), color = 'r')

ax[2].text(0.3, -0.06, r'$P_e = 10 GeV/c, P_{^{12}C}=600GeV/c$', color = 'b')

ax[1].legend(loc='upper right', shadow='True', fontsize='medium', numpoints=1)

#plotname='c12_asym_mulp.pdf'
plotname=filename
#plotname=plotname.replace('.dat','_mulp.png')
plotname=plotname.replace('.dat','_mulp.png')

plt.savefig(plotname,bbox_inches='tight')


# In[10]:

f1, axis = plt.subplots(1, 1, sharex=False, figsize=(8,6))

axis.errorbar(z_pip, Asym_H, yerr=Asym_Herr ,fmt='D', color='red',  label='e-p HERMES')

axis.errorbar(z, Asym, yerr=Astat ,fmt='o', color='blue',  label=r'e-p EIC')
axis.errorbar(z_c12, Asym_c12, yerr=Astat_c12 ,fmt='^', color='g',  label=r'e-$^{12}$C EIC')

axis.set_ylabel(r'$Asym=\frac{ N^{\pi^{+}}(z) - N^{\pi^{-}}(z)  }{ N^{\pi^{+}}(z) + N^{\pi^{-}}(z)  }$')
axis.set_title(r'Charge Asymmetry $(\pi^{+}-\pi^{-})$')
axis.set_xlabel(r'$z~(GeV/c)$')
axis.set_xlim(0.15,0.95)
#axis.set_ylim(-0.01,0.15)
axis.grid()

font_prop = font_manager.FontProperties(size=14)
axis.text(0.20, 0.19, r'EIC: $0.08<x_{B}<0.12, <Q^2>=%4.2f ~GeV^2$'%Q2_avg, color = 'b', fontproperties=font_prop)
axis.text(0.2, 0.22, r'HERMES: $<x_{B}>=%4.3f, <Q^2>=%4.2f ~GeV^2$'%(x_Havg,Q2_Havg), color = 'r', fontproperties=font_prop)
axis.text(0.4, -0.08, r'$P_e = 10 GeV/c, P_{^{12}C}=600GeV/c$', color = 'b', fontproperties=font_prop)

axis.legend(loc='upper left', shadow='True', fontsize='medium', numpoints=1,)
plotname=filename
#plotname=plotname.replace('.dat','_asym.png')
plotname=plotname.replace('.dat','_asym.pdf')
plotname=plotname.replace('pdf','png')

plt.savefig(plotname,bbox_inches='tight')


# In[9]:

f1, axis = plt.subplots(2, 1, sharex=False, figsize=(8,9))

axis[0].errorbar(z, Asym/Ax, yerr=Astat ,fmt='o', color='blue',  label=r'e-p EIC')
axis[0].errorbar(z_c12, Asym_c12/Ax_c12, yerr=Astat_c12 ,fmt='^', color='g',  label=r'e-$^{12}$C EIC')

axis[0].set_ylabel(r'$B(z) = Asym(x,z)/A(x)$')
#axis[0].set_title(r'B(z) = $Asym(x,z)/A(x)$')
axis[0].set_xlabel(r'$z~(GeV/c)$')
axis[0].set_xlim(0.15,0.95)
axis[0].set_ylim(-0.3,0.55)
axis[0].grid()
axis[0].legend(loc='upper left', shadow='True', fontsize='medium', numpoints=1)

axis[1].errorbar(z, Ax, yerr=Astat ,fmt='o', color='blue',  label=r'e-p EIC')
axis[1].errorbar(z_c12, Ax_c12, yerr=Astat_c12 ,fmt='^', color='g',  label=r'e-$^{12}$C EIC')

axis[1].set_ylabel(r'$A(x)$')
axis[1].set_xlabel(r'$z~(GeV/c)$')
axis[1].set_xlim(0.15,0.95)
#axis[1].set_ylim(-0.24,0.55)
axis[1].grid()

axis[0].text(0.34, -0.27, r'$0.08<x_{B}<0.12, <Q^2>=%4.2f ~GeV^2$'%Q2_avg, color = 'b', fontproperties=font_prop)
axis[1].text(0.3, 0.118, r'$P_e = 10 GeV/c, P_{^{12}C}=600GeV/c$', color = 'b', fontproperties=font_prop)


plotname=filename
#plotname=plotname.replace('.dat','_asym.png')
plotname=plotname.replace('.dat','_Bz.pdf')
plotname=plotname.replace('pdf','png')

plt.savefig(plotname,bbox_inches='tight')


# In[ ]:



