import matplotlib.pyplot as pl
import numpy as np

from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.stats import chisquare
# -----------------------------------------------------------------

# universal constants and units
c = 2.99792458e8        # m/s
h = 6.62606957e-34      # J s
k = 1.3806488e-23       # J/K
pc = 3.08567758e16      # m
Msun = 1.9891e30        # kg

# typical opacities at 350 micron
kappa350_Cortese = 0.192    # m2/kg (Cortese)
kappa350_Zubko = 0.330      # m2/kg (Zubko)

## This static function returns the black body emissivity \f$B_\nu\f$ (W/m2/Hz) for the specified wavelength or
# wavelengths \em wave (micron), and for the given dust temperature \em T (K).
def Bnu(wave, T):
    nu = c / (wave*1e-6)                                    # Hz
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    return Bnu

## This static function returns the grey body flux (Jy) for the specified wavelength or wavelengths \em wave (micron),
# for the given observer distance \em D (pc), power-law exponent \em beta, opacity \em kappa at 350 micron (m2/kg),
# dust temperature \em T (K) and dust mass \em M (Msun).
def greybody(D, beta, kappa350, wave, T, M):
    nu = c / (wave*1e-6)                                    # Hz
    nu350 = c / 350e-6                                      # Hz
    kappa = kappa350 * (nu/nu350)**beta                     # m2/kg
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    flux = M*Msun * kappa * Bnu / (D*pc)**2                 # W/m2/Hz
    return flux * 1e26 * 1e3                                # mJy

def M_d(D, beta, kappa_b3, wave, T, flux):
    flux *= 1e-26*1e-3                                      # mJy -> W/m2/Hz
    nu = c / (wave*1e-6)                                    # Hz
    nu520 = c / 520e-6                                      # Hz
    kappa = kappa_b3 * (nu/nu520)**beta                     # m2/kg
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    M = flux / (kappa * Bnu / (D*pc)**2)                    # kg
    return M/Msun                                           # Msun


#redshift of SPT0348
z = 5.656

# Global fluxes
#SPT
S_3mm = np.nan
e_S_3mm = np.nan
S_150GHz = 4.11 #raw=4.8
e_S_150GHz = 1.42 #raw=1.5
S_220GHz = 17.09#raw = 20.6
e_S_220GHz = 5.81#raw = 4.5
#ALMA?
S_890um = np.nan
e_S_890um = np.nan
#LABOCA
S_870um = 51.5
e_S_870um = S_870um*np.sqrt((4.2/S_870um)**2 + (0.1)**2)
#HERSCHEL-SPIRE
S_500um = 55.0
e_S_500um = S_500um*np.sqrt((6.8/S_500um)**2 +0.05**2)
S_350um = 45.1
e_S_350um = S_350um*np.sqrt((6.4/S_350um)**2 +0.05**2)
S_250um = 29.1
e_S_250um = S_250um*np.sqrt((6.1/S_250um)**2 +0.05**2)
#HERSCHEL-PACS
S_160um = np.nan
e_S_160um = 8.7 #upper limit
S_100um = np.nan
e_S_100um = 1.9 #upper limit

global_wavespace = np.array([c/(150e9)*1e6,c/(220e9)*1e6,870.,500.,350.,250.])/(1+z)
global_fluxes = np.array([S_150GHz,S_220GHz,S_870um,S_500um,S_350um,S_250um])
e_global_fluxes = np.array([e_S_150GHz,e_S_220GHz,e_S_870um,e_S_500um,e_S_350um,e_S_250um])
#USE THIS TO FIND BEST BETA
# test_beta = np.linspace(1.0,3.0,1000)
# chi2_global = []
# for b in test_beta:
#     chi2_global.append(chisquare(global_fluxes,greybody(d,b,kappa,global_wavespace,Td_global[1],Md_global[1]))[0])
# chi2_global = np.array(chi2_global)
# test_beta[np.where(chi2_global == min(chi2_global))[0]][0]





#b3_cyl3 = c/((115+84)/2.0*1e9)*1e6/(1+z)
b3= c/(93.8713e9)*1e6/(1+z)  #c/((86.7+88.6+98.7+100.8)/4*1e9)*1e6/(1+z)    #c/(86.9e9)*1e6 /(1+z) #3.5e3
b6= c/(226.954e9)*1e6/(1+z) #c/((217.9+220+234.0+235.9)/4*1e9)*1e6/(1+z)   #c/(220e9)*1e6 /(1+z) #1.4e3
b7= c/(290.392e9)*1e6/(1+z) #c/((283.3+285.2+295.3+297.4)/4*1e9)*1e6/(1+z) #c/(285.25e9)*1e6 /(1+z) #1.1e3

# opacity at Band-3 from Klaas+2001
beta_E = 2.0
beta_W = 2.0
dbeta = 0.2
beta_global = 2.0
dbeta_global = 0.2
kappa_b3 = 0.865*10*(b3/850.0)**beta_E # m2/kg (kappa_lamba1 = kappa_lambda2 * (lambda1/lambda2)^beta)


b3_E = 0.2
b6_E = 3.04
b7_E = 5.2
e_b3_E = np.sqrt((b3_E*0.1)**2 + 0.03**2)
e_b6_E = np.sqrt((b6_E*0.1)**2 + (0.1*np.sqrt(265.0/82.055))**2)
e_b7_E = np.sqrt((b7_E*0.1)**2 + 0.5**2)
# test_beta = np.linspace(1.0,3.0,1000)
# chi2_E = []
# for b in test_beta:
#     chi2_E.append(chisquare([b7_E,b6_E,b3_E],greybody(d,b,kappa,np.array([b7,b6,b3]),Td_E[1],Md_E[1]))[0])
# chi2_E = np.array(chi2_E)
# print test_beta[np.where(chi2_E == min(chi2_E))[0]][0]

b3_W = 0.54
b6_W = 9.6
b7_W = 18.9
e_b3_W = np.sqrt((b3_W*0.1)**2 + 0.04**2)
e_b6_W = np.sqrt((b6_W*0.1)**2 + (0.1*np.sqrt(408.0/82.055))**2)
e_b7_W = np.sqrt((b7_W*0.1)**2 + (0.6)**2)
beta_W=2.0
# test_beta = np.linspace(1.0,3.0,1000)
# chi2_W = []
# for b in test_beta:
#     chi2_W.append(chisquare([b7_W,b6_W,b3_W],greybody(d,b,kappa,np.array([b7,b6,b3]),Td_W[1],Md_W[1]))[0])
# chi2_W = np.array(chi2_W)
# print test_beta[np.where(chi2_W == min(chi2_W))[0]][0]


d = 55057.3e6
DL = d*3.086e16 #pc to meters
Mdust = 3.52139402499e40 / (2e30)
kappa = (kappa350_Cortese + kappa350_Zubko)/2
wavespace = np.logspace(1,4,500)
def greybody_W(wave,T,M):
    #beta = beta_W
    return greybody(d, beta, kappa, wave, T, M)


Td_E = []
e_Td_E = []
Md_E = []
e_Md_E = []
for i in range(3):
    if i == 0:
        beta=beta_E-dbeta
    if i == 1:
        beta=beta_E
    if i == 2:
        beta=beta_E+dbeta
    def greybody_E(wave,T,M):
        return greybody(d, beta, kappa, wave, T, M)

    params_E,cov_E = curve_fit(greybody_E,[b3,b6,b7],[b3_E,b6_E,b7_E],sigma=[e_b3_E,e_b6_E,e_b7_E],p0=[45,1e10]) #
    Td_E_tmp,Md_E_tmp = params_E[0],params_E[1]
    e_Td_E_tmp,e_Md_E_tmp = np.sqrt(np.diag(cov_E))
    Td_E.append(Td_E_tmp)
    e_Td_E.append(e_Td_E_tmp)
    Md_E.append(Md_E_tmp)
    e_Md_E.append(e_Md_E_tmp)
Td_E = np.array(Td_E)
e_Td_E = np.array(e_Td_E)
Md_E = np.array(Md_E)
e_Md_E = np.array(e_Md_E)
greybody_spt0348_E = greybody(d,beta_E,kappa,wavespace,Td_E[1],Md_E[1])
greybody_spt0348_E_upper = greybody(d,beta_E-dbeta,kappa,wavespace,Td_E[0],Md_E[0])
greybody_spt0348_E_lower = greybody(d,beta_E+dbeta,kappa,wavespace,Td_E[2],Md_E[2])

Td_E_upper = (Td_E - Td_E[1])[0]
Td_E_lower = (Td_E[1] - Td_E)[-1]


print 'SPT0348-E:'
print 'T_d = ' + str(int(round(Td_E[1],0))) + ' +' + str(int(round(Td_E_upper,0))) + '-' + str(int(round(Td_E_lower,0))) + ' K'
print 'peak flux = ' + str(round(np.max(greybody_spt0348_E),1)) + 'mJy at wavelength = ' + str(round(wavespace[[np.where(greybody_spt0348_E == np.max(greybody_spt0348_E))[0][0]]][0],1)) + 'um \n'

Td_W = []
e_Td_W = []
Md_W = []
e_Md_W = []
for i in range(3):
    if i == 0:
        beta=beta_W-dbeta
    if i == 1:
        beta=beta_W
    if i == 2:
        beta=beta_W+dbeta
    def greybody_W(wave,T,M):
        return greybody(d, beta, kappa, wave, T, M)

    params_W,cov_W = curve_fit(greybody_W,[b3,b6,b7],[b3_W,b6_W,b7_W],sigma=[e_b3_W,e_b6_W,e_b7_W],p0=[45,1e10]) #
    Td_W_tmp,Md_W_tmp = params_W[0],params_W[1]
    e_Td_W_tmp,e_Md_W_tmp = np.sqrt(np.diag(cov_W))
    Td_W.append(Td_W_tmp)
    e_Td_W.append(e_Td_W_tmp)
    Md_W.append(Md_W_tmp)
    e_Md_W.append(e_Md_W_tmp)
Td_W = np.array(Td_W)
e_Td_W = np.array(e_Td_W)
Md_W = np.array(Md_W)
e_Md_W = np.array(e_Md_W)
greybody_spt0348_W = greybody(d,beta_W,kappa,wavespace,Td_W[1],Md_W[1])
greybody_spt0348_W_upper = greybody(d,beta_W-dbeta,kappa,wavespace,Td_W[0],Md_W[0])
greybody_spt0348_W_lower = greybody(d,beta_W+dbeta,kappa,wavespace,Td_W[2],Md_W[2])
Td_W_upper = (Td_W - Td_W[1])[0]
Td_W_lower = (Td_W[1] - Td_W)[-1]

print 'SPT0348-W:'
print 'T_d = ' + str(int(round(Td_W[1],0))) + ' +' + str(int(round(Td_W_upper,0))) + '-' + str(int(round(Td_W_lower,0))) + ' K'
print 'peak flux = ' + str(round(np.max(greybody_spt0348_W),1)) + 'mJy at wavelength = ' + str(round(wavespace[[np.where(greybody_spt0348_W == np.max(greybody_spt0348_W))[0][0]]][0],1)) + 'um \n'


greybody_tot = greybody(d,beta_W,kappa,wavespace,Td_W[1],Md_W[1]) + greybody(d,beta_E,kappa,wavespace,Td_E[1],Md_E[1])
Td_tot = (Td_E[1]+Td_W[1])/2.0
Td_tot_upper = np.sqrt(Td_E_upper**2 + Td_W_upper**2)
Td_tot_lower = np.sqrt(Td_E_lower**2 + Td_W_lower**2)
print 'SPT0348-E + SPT0348-W:'
print 'T_d = ' + str(int(round(Td_tot,0))) + ' +' + str(int(round(Td_tot_upper,0))) + '-' + str(int(round(Td_tot_lower,0))) + ' K'
print 'peak flux = ' + str(round(np.max(greybody_tot),1)) + 'mJy at wavelength = ' + str(round(wavespace[[np.where(greybody_tot == np.max(greybody_tot))[0][0]]][0],1)) + 'um \n'

Td_global = []
e_Td_global = []
Md_global = []
e_Md_global = []
for i in range(3):
    if i == 0:
        beta=beta_global-dbeta_global
    if i == 1:
        beta=beta_global
    if i == 2:
        beta=beta_global+dbeta_global
    def greybody_global(wave,T,M):
        return greybody(d, beta, kappa, wave, T, M)

    params_global,cov_global = curve_fit(greybody_global,global_wavespace[:5],global_fluxes[:5],sigma=e_global_fluxes[:5],p0=[38,3e10]) #
    Td_global_tmp,Md_global_tmp = params_global[0],params_global[1]
    e_Td_global_tmp,e_Md_global_tmp = np.sqrt(np.diag(cov_global))
    Td_global.append(Td_global_tmp)
    e_Td_global.append(e_Td_global_tmp)
    Md_global.append(Md_global_tmp)
    e_Md_global.append(e_Md_global_tmp)
Td_global = np.array(Td_global)
e_Td_global = np.array(e_Td_global)
Md_global = np.array(Md_global)
e_Md_global = np.array(e_Md_global)

greybody_spt0348_global = greybody(d,beta_global,kappa,wavespace,Td_global[1],Md_global[1])
greybody_spt0348_global_upper = greybody(d,beta_global-dbeta_global,kappa,wavespace,Td_global[0],Md_global[0])
greybody_spt0348_global_lower = greybody(d,beta_global+dbeta_global,kappa,wavespace,Td_global[2],Md_global[2])
Td_global_upper = (Td_global - Td_global[1])[0]
Td_global_lower = (Td_global[1] - Td_global)[-1]


print 'Global SPT0348:'
print 'T_d = ' + str(int(round(Td_global[1],0))) + ' +' + str(int(round(Td_global_upper,0))) + '-' + str(int(round(Td_global_lower,0))) + ' K'
print 'peak flux = ' + str(round(np.max(greybody_spt0348_global),1)) + 'mJy at wavelength = ' + str(round(wavespace[[np.where(greybody_spt0348_global == np.max(greybody_spt0348_global))[0][0]]][0],1)) + 'um \n'




fig = pl.figure()
pl.subplots_adjust(hspace=0,wspace=0)
ax1 = fig.add_subplot(1,1,1)
ax1.set_xscale('log')
ax1.set_yscale('log')

label_W = r'${:}_{{-{:}}}^{{+{:}}}$ K'.format(int(round(Td_W[1],0)),int(round(Td_W_lower,0)),int(round(Td_W_upper,0)))
ax1.plot([b7,b6,b3], [b7_W,b6_W,b3_W],marker='.',c='r',zorder=10, linestyle='None',label='SPT0348-W')
ax1.errorbar([b7,b6,b3], [b7_W,b6_W,b3_W], yerr=[e_b7_W,e_b6_W,e_b3_W],linestyle='None',c='r',zorder=10)
ax1.loglog(wavespace, greybody_spt0348_W,alpha=0.5, c='r')
ax1.fill_between(wavespace,greybody_spt0348_W_upper,greybody_spt0348_W_lower,color='r',alpha=0.2, label=r'$\beta=$'+str(beta_W)+r'$\pm$'+str(dbeta) + r', $T_d = $' + label_W)
#ax1.scatter(c/(343e9)*1e6/(1+z),22,marker='o',c='m',zorder=10,linestyle='None',label='Cycle3 Band3 W')

label_E = r'${:}_{{-{:}}}^{{+{:}}}$ K'.format(int(round(Td_E[1],0)),int(round(Td_E_lower,0)),int(round(Td_E_upper,0)))
ax1.plot([b7,b6,b3],[b7_E,b6_E,b3_E],marker='.',c='b',linestyle='None',zorder=10, label='SPT0348-E')
ax1.errorbar([b7,b6,b3],[b7_E,b6_E,b3_E], yerr=[e_b7_E,e_b6_E,e_b3_E],linestyle='None',c='b',zorder=10)
ax1.loglog(wavespace, greybody_spt0348_E,alpha=0.5, c='b')
ax1.fill_between(wavespace,greybody_spt0348_E_upper,greybody_spt0348_E_lower,color='b',alpha=0.2, label=r'$\beta=$'+str(beta_E)+r'$\pm$'+str(dbeta)+ r', $T_d = $' + label_E)
#ax1.scatter(c/(343e9)*1e6/(1+z),5.7,marker='o',c='m',zorder=10,linestyle='None',label='Cycle3 Band3 E')

label_tot = r'${:}_{{-{:}}}^{{+{:}}}$ K'.format(int(round(Td_tot,0)),int(round(Td_tot_lower,0)),int(round(Td_tot_upper,0)))
ax1.loglog(wavespace, greybody_tot,alpha=1, c='c')
ax1.errorbar([b7,b6,b3],[b7_W+b7_E,b6_W+b6_E,b3_W+b3_E], yerr=[np.sqrt(e_b7_E**2+e_b7_W**2),np.sqrt(e_b6_E**2+e_b6_W**2),np.sqrt(e_b3_E**2+e_b3_W**2)],linestyle='None',c='c',zorder=10)
ax1.scatter([b7,b6,b3], [b7_W+b7_E,b6_W+b6_E,b3_W+b3_E],marker='s',c='c',linestyle='None',zorder=10,label=r'W+E')#: $T_d = $' + label_tot)

label_global = r'${:}_{{-{:}}}^{{+{:}}}$ K'.format(int(round(Td_global[1],0)),int(round(Td_global_lower,0)),int(round(Td_global_upper,0)))
#ax1.scatter(global_wavespace,global_fluxes,marker='*',c='purple', zorder=10,linestyle='None')#, label=r'global: $\beta=$'+str(beta_global)+r'$\pm$'+str(dbeta))
#ax1.errorbar(global_wavespace,global_fluxes,yerr=e_global_fluxes,linestyle='None',c='purple')
ax1.scatter(global_wavespace,global_fluxes,marker='*',c='k', zorder=10,linestyle='None', label=r'global: $\beta=$'+str(beta_global)+r'$\pm$'+str(dbeta))
ax1.errorbar(global_wavespace,global_fluxes,yerr=e_global_fluxes,linestyle='None',c='k',zorder=10)
# ax1.scatter(global_wavespace[:5],global_fluxes[:5],marker='*',c='k', zorder=10,linestyle='None', label=r'global: $\beta=$'+str(beta_global)+r'$\pm$'+str(dbeta))
# ax1.errorbar(global_wavespace[:5],global_fluxes[:5],yerr=e_global_fluxes[:5],linestyle='None',c='k')
ax1.loglog(wavespace, greybody_spt0348_global,alpha=1, c='k',zorder=10)
ax1.fill_between(wavespace,greybody_spt0348_global_upper,greybody_spt0348_global_lower,color='k',alpha=0.2, label=r'$T_d = $' + label_global)


'''plot global flux points'''

ax1.set_xlabel(r'rest wavelength ($\mu$m)')
ax1.set_ylabel(r'flux (mJy)')
ax1.legend(ncol=1,fontsize=8,loc=0)
low_xlim,high_xlim = 1e1,2000
ax1.set_ylim(1e-2,200)
ax1.set_xlim(low_xlim,high_xlim)

# ax1.axvline(x=8,ls='--',c='k')
# ax1.axvline(x=1000,ls='--',c='k')

ax2 = ax1.twiny()
ax2.set_xscale('log')
ax2.set_xlim(low_xlim*(1+z),high_xlim*(1+z))
ax2.set_xlabel(r'observed wavelength ($\mu$m)')

pl.savefig('../Figures/SED/SPT0348-EW_greybody_Td20.pdf')
pl.close()



'''integrate from 8-1000um to get L_TIR'''

def greybody_int(nu,D, beta, kappa350, T, M):
    #nu = c / (wave*1e-6)                                    # Hz
    nu350 = c / 350e-6                                      # Hz
    kappa = kappa350 * (nu/nu350)**beta                     # m2/kg
    Bnu = 2*h*nu**3/ c**2 / (np.exp((h*nu)/(k*T)) - 1)      # W/m2/Hz
    flux = M*Msun * kappa * Bnu / (D*pc)**2                 # W/m2/Hz
    return flux                                             # W/m2/Hz

# nuspace = np.linspace(nu_high,nu_low,1000)
# pl.loglog(nuspace, greybody_int(nuspace,d, 1.8, kappa, Td_W[1], Md_W[1]),alpha=1, c='orange')
# pl.savefig('../Figures/test_nu.pdf')
# pl.close()


Lsun = 3.828e26 #W
#Bounds given by Kennicutt & Evans (2012) - eq. 12 & Table 1
nu_low = c/(3e-6)#/(1+z)
nu_high = c/(1100e-6)#/(1+z)
L_E = quad(greybody_int,nu_high,nu_low,args=(d,beta_E,kappa,Td_E[1],Md_E[1]))[0] *DL**2*4*np.pi /(1+z)#W
L_E_low = quad(greybody_int,nu_high,nu_low,args=(d,beta_E+dbeta,kappa,Td_E[2],Md_E[2]))[0] *DL**2*4*np.pi /(1+z)#W
L_E_high = quad(greybody_int,nu_high,nu_low,args=(d,beta_E-dbeta,kappa,Td_E[0],Md_E[0]))[0] *DL**2*4*np.pi /(1+z)#W

L_W = quad(greybody_int,nu_high,nu_low,args=(d,beta_W,kappa,Td_W[1],Md_W[1]))[0] *DL**2*4*np.pi /(1+z)#W
L_W_low = quad(greybody_int,nu_high,nu_low,args=(d,beta_W+dbeta,kappa,Td_W[2],Md_W[2]))[0] *DL**2*4*np.pi /(1+z)#W
L_W_high = quad(greybody_int,nu_high,nu_low,args=(d,beta_W-dbeta,kappa,Td_W[0],Md_W[0]))[0] *DL**2*4*np.pi /(1+z)#W

L_global = quad(greybody_int,nu_high,nu_low,args=(d,beta_global,kappa,Td_global[1],Md_global[1]))[0] *DL**2*4*np.pi /(1+z)#W
L_global_low = quad(greybody_int,nu_high,nu_low,args=(d,beta_global+dbeta_global,kappa,Td_global[2],Md_global[2]))[0] *DL**2*4*np.pi /(1+z)#W
L_global_high = quad(greybody_int,nu_high,nu_low,args=(d,beta_global-dbeta_global,kappa,Td_global[0],Md_global[0]))[0] *DL**2*4*np.pi /(1+z)#W


print 'L_E = (' + str(round(L_E/Lsun*1e-13,2)) + ' +' + str(round((L_E_high-L_E)/Lsun*1e-13,1)) + '-' + str(round((L_E-L_E_low)/Lsun*1e-13,1)) + ') x10^13 Lsun'
print 'L_W = (' + str(round(L_W/Lsun*1e-13,2)) + ' +' + str(round((L_W_high-L_W)/Lsun*1e-13,1)) + '-' + str(round((L_W-L_W_low)/Lsun*1e-13,1)) + ') x10^13 Lsun'
print 'L_E+L_W = (' + str(round((L_E+L_W)/Lsun*1e-13,2)) + ' +' + str(round(np.sqrt((L_E_high-L_E)**2 + (L_W_high-L_W)**2)/Lsun*1e-13,1)) + '-' +str(round(np.sqrt((L_E_low-L_E)**2 + (L_W_low-L_W)**2)/Lsun*1e-13,1)) + ') x10^13 Lsun'
print 'L_global = (' + str(round(L_global/Lsun*1e-13,2)) + ' +' + str(round((L_global_high-L_global)/Lsun*1e-13,1)) + '-' + str(round((L_global-L_global_low)/Lsun*1e-13,1)) + ') x10^13 Lsun'

print 'L_global / L_W+E = ' + str(round(L_global/(L_W + L_E),2)) +'\n'#+ ' +' + str(round(np.sqrt((7.1/3.24)**2 + (4.0/4.1)**2) * L_global/(L_W + L_E),1)) + '-' + str(round(np.sqrt((1.8/3.24)**2 + (1.9/4.1)**2) * L_global/(L_W + L_E),1)) + '\n'

# Kennicutt & Evans (2012) - eq. 12 & Table 1
def SFR(L_TIR):
    L_TIR *= 1e7
    return 10**(np.log10(L_TIR) - 43.41) #Msun/yr

SFR_E = SFR(L_E)
SFR_E_low = SFR(L_E)-SFR(L_E_low)
SFR_E_high = SFR(L_E_high)-SFR(L_E)
SFR_W = SFR(L_W)
SFR_W_low = SFR(L_W)-SFR(L_W_low)
SFR_W_high = SFR(L_W_high)-SFR(L_W)
SFR_global = SFR(L_global)
SFR_global_low = SFR(L_global)-SFR(L_global_low)
SFR_global_high = SFR(L_global_high)-SFR(L_global)

print 'SFR(E) = ' + str(int(round(SFR_E,0))) + ' +' + str(int(round(SFR_E_high,0))) + '-' + str(int(round(SFR_E_low,0))) + ' Msun/yr'
print 'SFR(W) = ' + str(int(round(SFR_W,0))) + ' +' + str(int(round(SFR_W_high,0))) + '-' + str(int(round(SFR_W_low,0))) + ' Msun/yr'
print 'SFR(E+W) = ' + str(int(round(SFR_E+SFR_W,0))) + ' +' + str(int(round(np.sqrt(SFR_E_high**2 + SFR_W_high**2),0))) + '-' + str(int(round(np.sqrt(SFR_E_low**2 + SFR_W_low**2),0))) + ' Msun/yr'
print 'SFR(global) = ' + str(int(round(SFR_global,0))) + ' +' + str(int(round(SFR_global_high,0))) + '-' + str(int(round(SFR_global_low,0))) + ' Msun/yr \n'

Md_spt0348_E = M_d(d,beta_E,kappa_b3,b3,Td_E[1],b3_E)
e_Md_spt0348_E = Md_spt0348_E - M_d(d,beta_E,kappa_b3,b3,Td_E[1],b3_E-e_b3_E)
Md_spt0348_W = M_d(d,beta_W,kappa_b3,b3,Td_W[1],b3_W)
e_Md_spt0348_W = Md_spt0348_W - M_d(d,beta_W,kappa_b3,b3,Td_W[1],b3_W-e_b3_W)
print 'M_d(E) = (' + str(round(Md_spt0348_E*1e-8,1)) + ' +/-' + str(round(e_Md_spt0348_E*1e-8,1)) +') x10^8 Msun'
print 'M_d(W) = (' + str(round(Md_spt0348_W*1e-8,1)) + ' +/-' + str(round(e_Md_spt0348_W*1e-8,1)) +') x10^8 Msun \n'

H2_W = 6.1e10 *(1/0.39)
H2_E = 1.4e10 *(1/0.39)
e_H2_W = np.array([9.32e+10,4.18e+11])
e_H2_E = np.array([1.66e+10,7.70e+10])

print 'H2/Md(E) = ' + str(round(H2_E / Md_spt0348_E,1)) + ' +' + str(int(round(H2_E /Md_spt0348_E * np.sqrt((e_H2_E[1]/H2_E)**2 + (e_Md_spt0348_E/Md_spt0348_E)**2),0))) + '-' + str(int(round(H2_E /Md_spt0348_E * np.sqrt((e_H2_E[0]/H2_E)**2 + (e_Md_spt0348_E/Md_spt0348_E)**2),0)))
print 'H2/Md(W) = ' + str(round(H2_W / Md_spt0348_W,1)) + ' +' + str(int(round(H2_W /Md_spt0348_W * np.sqrt((e_H2_W[1]/H2_W)**2 + (e_Md_spt0348_W/Md_spt0348_W)**2),0))) + '-' + str(int(round(H2_W /Md_spt0348_W * np.sqrt((e_H2_W[0]/H2_W)**2 + (e_Md_spt0348_W/Md_spt0348_W)**2),0)))
