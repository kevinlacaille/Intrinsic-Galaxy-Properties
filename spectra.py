import numpy as np
import aplpy as ap
#import montage_wrapper as montage
import pyfits as py
import pylab as pl
import img_scale
from astropy.visualization import make_lupton_rgb
import coords as co





# moving average function
def movingaverage(interval, window_size):
	window = np.ones(int(window_size))/float(window_size)
	return np.convolve(interval, window, 'same')


########################
# where to look for data
########################

CO54_dir = '../../../../../Research/SPT0348/Data/ALMA/Cycle_4/Lines/spt0348_CO54_dirty_contsub_briggs_robust05.fits'
NII_dir = '../../../../../Research/SPT0348/Data/ALMA/Cycle_4/Lines/spt0348_N+_bin5_dirty_contsub_briggs_robust05.fits'
CII_dir = '../../../../../Research/SPT0348/Data/ALMA/Cycle_4/Lines/spt0348_C+_dirty_contsub_briggs_robust05.fits'

################
# get data cubes
################

CO54 = np.nan_to_num(py.getdata(CO54_dir))[0]
CO54_vel = np.linspace(-4203.73,2667.95,128)
CO54_freq = np.linspace(87.7771,85.8082,128)

NII = np.nan_to_num(py.getdata(NII_dir))[0]
NII_vel = np.linspace(636.209,-1924.42,25)
NII_freq = np.linspace(219.055,220.93,25)

CII = np.nan_to_num(py.getdata(CII_dir))[0]
CII_vel = np.linspace(-823.058,1260.55,128)
CII_freq = np.linspace(286.321,284.337,128)


###############
# get spectra W
###############
CO54_spectra_W = []
nchans_CO54 = np.shape(CO54)[0]
for i in range(nchans_CO54):
	CO54_spectra_W.append(CO54[i][447][448])

NII_spectra_W = []
nchans_NII = np.shape(NII)[0]
for i in range(nchans_NII):
	NII_spectra_W.append(NII[i][361][362])

CII_spectra_W = []
nchans_CII = np.shape(CII)[0]
for i in range(nchans_CII):
	CII_spectra_W.append(CII[i][360][359])

CO54_spectra_W = np.array(CO54_spectra_W)
NII_spectra_W = np.array(NII_spectra_W)
CII_spectra_W = np.array(CII_spectra_W)

CO54_spectra_smooth_W = movingaverage(CO54_spectra_W,4)
NII_spectra_smooth_W = movingaverage(NII_spectra_W,4)
CII_spectra_smooth_W = movingaverage(CII_spectra_W,4)



###############
# get spectra E
###############
CO54_spectra_E = []
nchans_CO54 = np.shape(CO54)[0]
for i in range(nchans_CO54):
	CO54_spectra_E.append(CO54[i][454][425])

NII_spectra_E = []
nchans_NII = np.shape(NII)[0]
for i in range(nchans_NII):
	NII_spectra_E.append(NII[i][377][315])

CII_spectra_E = []
nchans_CII = np.shape(CII)[0]
for i in range(nchans_CII):
	CII_spectra_E.append(CII[i][376][318])

CO54_spectra_E = np.array(CO54_spectra_E)
NII_spectra_E = np.array(NII_spectra_E)
CII_spectra_E = np.array(CII_spectra_E)

CO54_spectra_smooth_E = movingaverage(CO54_spectra_E,4)
NII_spectra_smooth_E = movingaverage(NII_spectra_E,3)
CII_spectra_smooth_E = movingaverage(CII_spectra_E,4)




################
# plot spectra W
################
ax1 = pl.subplot(111)
fig = pl.figure(figsize=(8,5))
pl.rc('font',size=18)
pl.rc('mathtext', default='regular')

pl.xlim(-2000,2000)
pl.ylim(-0.4,1.2)

# pl.plot(CO54_vel, CO54_spectra_smooth_W/max(CO54_spectra_smooth_W), 'r-', label='CO(5-4)')
# pl.plot(NII_vel, NII_spectra_smooth_W/max(NII_spectra_smooth_W), 'g-', label='[NII]')
# pl.plot(CII_vel, CII_spectra_smooth_W/max(CII_spectra_smooth_W), 'b-', label='[CII]')
pl.plot(CO54_vel, CO54_spectra_smooth_W/max(CO54_spectra_smooth_W), c='r',ls='steps', label='CO(5-4)')
pl.plot(NII_vel, NII_spectra_smooth_W/max(NII_spectra_smooth_W), c='g',ls='steps', label='[NII]')
pl.plot(CII_vel, CII_spectra_smooth_W/max(CII_spectra_smooth_W), c='b',ls='steps', label='[CII]')


pl.text(-1750, 1, 'SPT0348-W', size=18)

pl.ylabel('Fraction of peak')
pl.xlabel(r'Radio velocity (km s$^{-1}$)')

pl.legend(fontsize=14,loc=0,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

pl.savefig('spectra_smooth_W.pdf', bbox_inches='tight')
pl.close()



################
# plot spectra E
################

ax2 = pl.subplot(111)
fig = pl.figure(figsize=(8,5))
pl.rc('font',size=18)
pl.rc('mathtext', default='regular')

pl.xlim(-2000,2000)
pl.ylim(-0.4,1.2)

pl.plot(CO54_vel, CO54_spectra_smooth_E/max(CO54_spectra_smooth_E), c='r',ls='steps', label='CO(5-4)')
pl.plot(NII_vel, NII_spectra_smooth_E/max(NII_spectra_smooth_E), c='g',ls='steps', label='[NII]')
pl.plot(CII_vel, CII_spectra_smooth_E/max(CII_spectra_smooth_E), c='b',ls='steps', label='[CII]')

pl.text(-1750, 1, 'SPT0348-E', size=18)

pl.ylabel('Fraction of peak')
pl.xlabel(r'Radio velocity (km s$^{-1}$)')

pl.legend(fontsize=14,loc=0,ncol=1,frameon=True,numpoints=1,scatterpoints=1)

pl.savefig('spectra_smooth_E.pdf', bbox_inches='tight')
pl.close()






