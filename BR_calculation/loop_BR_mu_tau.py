import math
from math import log10, exp, cos, acos, sin, asin, sqrt, radians

import pandas as pd
from scipy import interpolate

import csv
#import vegas
import numpy as np
from scipy.integrate import quad
import scipy.special as special

import scipy.interpolate
import pylab
import matplotlib as mpl

from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})
mpl.rcParams['xtick.major.pad'] = 12
mpl.rcParams['ytick.major.pad'] = 12


#########################################################
# Factores de Mixing # 
#########################################################
############ DIRAC ##############
Dfe = pd.read_csv('Dir_e.csv')
Dxe = Dfe['MN'].to_numpy()
Dye = Dfe['N_e'].to_numpy()
DNe = interpolate.InterpolatedUnivariateSpline(Dxe, Dye)

Dfmu = pd.read_csv('Dir_mu.csv')
Dxmu = Dfmu['MN'].to_numpy()
Dymu = Dfmu['N_mu'].to_numpy()
DNmu = interpolate.InterpolatedUnivariateSpline(Dxmu, Dymu)

Dftau = pd.read_csv('Dir_tau.csv')
Dxtau = Dftau['MN'].to_numpy()
Dytau = Dftau['N_tau'].to_numpy()
DNtau = interpolate.InterpolatedUnivariateSpline(Dxtau, Dytau)

############ MAJORANA ##############
Mfe = pd.read_csv('Maj_e.csv')
Mxe = Mfe['MN'].to_numpy()
Mye = Mfe['N_e'].to_numpy()
MNe = interpolate.InterpolatedUnivariateSpline(Mxe, Mye)

Mfmu = pd.read_csv('Maj_mu.csv')
Mxmu = Mfmu['MN'].to_numpy()
Mymu = Mfmu['N_mu'].to_numpy()
MNmu = interpolate.InterpolatedUnivariateSpline(Mxmu, Mymu)

Mftau = pd.read_csv('Maj_tau.csv')
Mxtau = Mftau['MN'].to_numpy()
Mytau = Mftau['N_tau'].to_numpy()
MNtau = interpolate.InterpolatedUnivariateSpline(Mxtau, Mytau)

#print(MNtau(4))

#########################################################
# definimos una funcion que permite cubrir un rango de valores # 
#########################################################

def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L

########################################################
##############  Defino Constantes   #################### 
########################################################
PI = 3.14159;
GF = 1.1664*10**-5.0;

#Here eff stands for the efficiency
eff = 1.0;

########################################################
MBc = 6.275;                          
FBc = 0.322;
Vcb = 0.04;
########################################################

#BeN2 = 10**-8
#BmuN2 = 5*10**-7
#BtauN2 = 5*10**-6



# Acceptance Factor
def AF_Maj(MN):
  return 1.0# - np.exp(-LD*GN_Maj(MN)/albe)
  
def AF_Dir(MN):
  return 1.0# - np.exp(-LD*GN_Dir(MN)/albe)



BR_Maj = np.zeros((81,81,42))
BR_Dir = np.zeros((81,81,42))
ctau_Maj_array = np.zeros((81,81,42))
ctau_Dir_array = np.zeros((81,81,42))

#massPoints = [2.2,3.0,4.0,5.0]

#BR_Maj = np.zeros((81,81,4))
#BR_Dir = np.zeros((81,81,4))
#GN_Maj_array = np.zeros((81,81,4))
#GN_Dir_array = np.zeros((81,81,4))




Me = 0.511*10**-3;
Mmu = 105.6*10**-3;
MT = 1.777;

M1 = Mmu;
M2 = Mmu;
ML = MT;
########################################################
###############  Defino Funciones   #################### 
########################################################
# LD: Largo del detector en GeV^-1, Donde L esta esta mm
L = 1000; #mm
LD = L*5.0*(10**12); #GeV
albe = 2.0; #parametro alfa beta cineticos

#B Mesons Total Decay Width 
GBc = 1.29*10**-12 #GeV
GB = 4.017*10**-13 #GeV

########################################################
########################################################
########################################################
##############  Funcion de Integracion  ################ 
########################################################

def ALAMsqr(x, y, z):
  return ((x**2.0 + y**2.0 + z**2.0 - 2.0*x*y - 2.0*x*z - 2.0*z*y))**(0.5)	

MN_count = [];

#GammaBc_LNV = []
#GammaBc_LNC = []

#GammaBc_Maj = []
#GammaBc_Dir = []

#Br_eff_Majorana = []
#Br_eff_Dirac = []

MN_counter = 0
for MN in frange(1.9, 6.1, 0.1):
#for MN in massPoints:
  MN_counter += 1
  MN_count.append(MN)
  Int_LNV = None
  Int_LNC = None
  #Br_eff_Dir = None
  #Br_eff_Maj = None
  
  ########################################################
  ##################  Funciones a Usar   ################# 
  ########################################################
  
    
  #LNV		
  #ZLNV = ((GF**4*BmuN2**2*Vcb**2*MN*FBc**2)/((2*PI)**4*MBc**3))*ALAMsqr(MBc**2,MN**2,M1**2)
  ZLNV = ((GF**4*Vcb**2*MN*FBc**2)/((2*PI)**4*MBc**3))*ALAMsqr(MBc**2,MN**2,M1**2)
  
  def dGdEL_LNV(x1):
    return 2*ZLNV*(1/(2*MN))*(MBc**2*(MN**2+M1**2)-(MN**2-M1**2)**2)*x1*sqrt(x1**2-ML**2)* \
		((MN**2-2*MN*x1+ML**2-M2**2)**2/(MN**2-2*MN*x1+ML**2))
  Int_LNV = quad(dGdEL_LNV,ML,(MN**2+ML**2-M1**2)/(2*MN))
  
  #GammaBc_LNV.append(Int_LNV[0])
	#print("Br_LNV:", Int_LNV[0], "error:", Int_LNV[1])
  
  #LNC
  #ZLNC = ((GF**4*BmuN2*BtauN2*Vcb**2*MN*FBc**2)/((2*PI)**4*MBc**3))*ALAMsqr(MBc**2,MN**2,M1**2)
  ZLNC = ((GF**4*Vcb**2*MN*FBc**2)/((2*PI)**4*MBc**3))*ALAMsqr(MBc**2,MN**2,M1**2)
  
  #print(MN, ZLNC, ZLNV)
  
  def dGdEL_LNC(x1):
    return (2*ZLNC/(96*MN**2*(ML**2+MN*(MN-2*x1))**3))*( \
			8*sqrt(x1**2-ML**2)*MN*(M2**2-ML**2+(2*x1-MN)*MN)**2 \
			*(-M1**4+MBc**2*MN**2-MN**4+M1**2*(MBc**2+2*MN**2)) \
			*(8*x1**3*MN**2-2*ML**2*MN*(2*M2**2+ML**2+MN**2) \
			-2*x1**2*MN*(M2**2+5*(ML**2+MN**2))+x1*(3*ML**4+10*ML**2*MN**2+3*MN**4+3*M2**2*(ML**2+MN**2))))
  Int_LNC = quad(dGdEL_LNC,ML,(MN**2+ML**2-M1**2)/(2*MN))
  #GammaBc_LNC.append(Int_LNC[0])
  #print("Br_LNC:", Int_LNC[0], "error:", Int_LNC[1])
  
  #if (Int_LNC[0] > 0 and Int_LNV[0] > 0):
  
    
  for logAlphaMu in range(81):
    for logAlphaTau in range(81):
      BeN2 = 0
      BmuN2 = 10**(-logAlphaMu/10)
      BtauN2 = 10**(-logAlphaTau/10)
      
      GN_Maj = (MNe(MN)*BeN2+MNmu(MN)*BmuN2+MNtau(MN)*BtauN2)*(GF**2*MN**5)/(96*PI**3)
      
      GN_Dir = (DNe(MN)*BeN2+DNmu(MN)*BmuN2+DNtau(MN)*BtauN2)*(GF**2*MN**5)/(96*PI**3)
      
      BR_Maj[logAlphaMu][logAlphaTau][MN_counter-1] = (Int_LNC[0]*BmuN2*BtauN2 + Int_LNV[0]*BmuN2**2)/(GN_Maj*GBc)
      BR_Dir[logAlphaMu][logAlphaTau][MN_counter-1] = (Int_LNC[0]*BmuN2*BtauN2)           /(GN_Dir*GBc)
      ctau_Maj_array[logAlphaMu][logAlphaTau][MN_counter-1] = (197*10**(-15))/GN_Maj
      ctau_Dir_array[logAlphaMu][logAlphaTau][MN_counter-1] = (197*10**(-15))/GN_Dir
      #BR_Maj[logAlpha][MN_counter-1] = (Int_LNC[0]+Int_LNV[0])/(GBc)
      #BR_Dir[logAlpha][MN_counter-1] = (Int_LNC[0])           /(GBc)


np.savetxt('BR_Majorana_mu_tau.txt', BR_Maj.reshape(BR_Maj.shape[0], -1), delimiter=',')
np.savetxt('BR_Dirac_mu_tau.txt', BR_Dir.reshape(BR_Dir.shape[0], -1), delimiter=',')
np.savetxt('ctau_Majorana_mu_tau.txt', ctau_Maj_array.reshape(ctau_Maj_array.shape[0], -1), delimiter=',')
np.savetxt('ctau_Dirac_mu_tau.txt', ctau_Dir_array.reshape(ctau_Dir_array.shape[0], -1), delimiter=',')


'''
# Debo agregar GN para HN tipo Dirac, y GBc.
	Br_eff_Dir = eff*AF_Dir(MN)*(Int_LNC[0]/(GN_Dir(MN)*GBc))
	Br_eff_Dirac.append(Br_eff_Dir)
	print("Dirac:", MN, "Br:", Br_eff_Dir)

# Debo agregar GN para HN tipo Majorana, y GBc.
	Br_eff_Maj = eff*AF_Maj(MN)*(Int_LNC[0]+Int_LNV[0])/(GN_Maj(MN)*GBc)
	Br_eff_Majorana.append(Br_eff_Maj)
	print("Majorana:", MN, "Br:", Br_eff_Maj)
####################################
#############  PLOTS  ##############
####################################
# Funciones para graficar
plt.plot(MN_count,Br_eff_Majorana,'-',linewidth=2, color="red", label=r'${\rm Majorana}$')
plt.plot(MN_count,Br_eff_Dirac,'--',linewidth=2, color="Blue", label=r'${\rm Dirac}$')

###################################################################################
###################################################################################
plt.xlim(2.1, 6.1)
#plt.ylim(1.0*10**(-13), 1.0*10**(-7))
	
#	ax.xaxis.set_label_coords(0.5, -0.15)
plt.xlabel(r'${m_N\ (GeV)}}$' ,fontsize=21,fontweight='bold')
plt.ylabel(r'${\rm Br( B^{+})}$' ,fontsize=21,fontweight='bold')
	
plt.legend()
	
#plt.yscale('log')
plt.savefig('Plot_Br_Maryam.pdf')
plt.show()'''