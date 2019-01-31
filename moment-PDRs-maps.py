import aplpy as ap
import numpy as np
import pylab as pl
import coords as co


mom0_map = '../moment_maps/spt0348_C+_dirty_contsub_briggs_robust05.image.mom0.fits'
mom1_map = '../moment_maps/spt0348_C+_dirty_contsub_briggs_robust05.image.mom1.fits'
mom2_map = '../moment_maps/spt0348_C+_dirty_contsub_briggs_robust05.image.mom2.fits'



############
#MOMENT MAPS
############

#set up figure
fig = pl.figure()
pl.rc('font',size=30)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('03:48:42.271')
DEC_c = co.convDMS('-62:20:50.85')


#moment-0 map
mom0 = ap.FITSFigure(mom0_map, figure = fig, subplot = [0,0,0.75,1])
mom0.tick_labels.set_xformat('hh:mm:ss')
mom0.tick_labels.set_yformat('dd:mm:ss')
mom0.hide_xaxis_label()
mom0.hide_xtick_labels()
mom0.hide_yaxis_label()
mom0.hide_ytick_labels()
mom0.recenter(RA_c, DEC_c,width = 5./3600, height= 5./3600)
mom0.show_colorscale(cmap='jet',vmax=12,vmin=0.0,interpolation='nearest')
mom0.add_colorbar(axis_label_text=r'S$_{[CII]}$ (Jy beam$^{-1}$ km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

mom0.add_beam()
mom0.beam.set_edgecolor('black')
mom0.beam.set_facecolor('none')
mom0.beam.set_hatch('//')
mom0.add_scalebar(1.0/3600,linewidth=3)
mom0.scalebar.set_label('6kpc')
mom0.scalebar.set_color('black')


#moment-1 map
mom1 = ap.FITSFigure(mom1_map, figure = fig, subplot = [1,0,0.75,1])
mom1.tick_labels.set_xformat('hh:mm:ss')
mom1.tick_labels.set_yformat('dd:mm:ss')
mom1.hide_xaxis_label()
mom1.hide_xtick_labels()
mom1.hide_yaxis_label()
mom1.hide_ytick_labels()
mom1.recenter(RA_c, DEC_c,width = 5./3600, height= 5./3600)
mom1.show_colorscale(cmap='jet',vmax=100,vmin=-500,interpolation='nearest')
mom1.add_colorbar(axis_label_text=r'S$_{[CII]}$ (km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

mom1.add_beam()
mom1.beam.set_edgecolor('black')
mom1.beam.set_facecolor('none')
mom1.beam.set_hatch('//')
mom1.add_scalebar(1.0/3600,linewidth=3)
mom1.scalebar.set_label('6kpc')
mom1.scalebar.set_color('black')


#moment-2 map
mom2 = ap.FITSFigure(mom2_map, figure = fig, subplot = [2,0,0.75,1])
mom2.tick_labels.set_xformat('hh:mm:ss.s')
mom2.tick_labels.set_yformat('dd:mm:ss')
mom2.hide_xaxis_label()
mom2.hide_xtick_labels()
mom2.hide_yaxis_label()
mom2.hide_ytick_labels()
mom2.recenter(RA_c, DEC_c,width = 5./3600, height= 5./3600)
mom2.show_colorscale(cmap='jet',vmax=300,vmin=0,interpolation='nearest')
mom2.add_colorbar(axis_label_text=r'S$_{[CII]}$ (km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

mom2.add_beam()
mom2.beam.set_edgecolor('black')
mom2.beam.set_facecolor('none')
mom2.beam.set_hatch('//')
mom2.add_scalebar(1.0/3600,linewidth=3)
mom2.scalebar.set_label('6kpc')
mom2.scalebar.set_color('black')


mom0.save('spt0348_moms.pdf',dpi=250)
#pl.show()
pl.close()
mom0.close()
mom1.close()
mom2.close()











band3_map = '../moment_maps/mom0s/spt0348_band3_spw0_clean1000_contsub_spw0_2sig.image.mom0_2sig.fits'
band6_map = '../moment_maps/mom0s/spt0348_band6_spw1_clean1000_contsub_bin5.image.mom0_2sig.fits'
band7_map = '../moment_maps/mom0s/spt0348_band7_spw1_clean1000_contsub_2sig.image.mom0.fits'

#set up figure
fig = pl.figure()
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('03:48:42.271')
DEC_c = co.convDMS('-62:20:50.85')


#band3
band3 = ap.FITSFigure(band3_map, figure = fig, subplot = [0,0,0.75,1])
band3.tick_labels.set_xformat('hh:mm:ss')
band3.tick_labels.set_yformat('dd:mm:ss')
band3.hide_xaxis_label()
band3.hide_xtick_labels()
band3.hide_yaxis_label()
band3.hide_ytick_labels()
band3.recenter(RA_c, DEC_c,width = 5./3600, height= 5./3600)
band3.show_colorscale(cmap='jet',vmax=1.38251,vmin=0.0,interpolation='nearest')
band3.add_colorbar(axis_label_text=r'S$_{CO(5-4)}$ (Jy beam$^{-1}$ km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

band3.add_beam()
band3.beam.set_edgecolor('black')
band3.beam.set_facecolor('none')
band3.beam.set_hatch('//')
band3.add_scalebar(1.0/3600,linewidth=3)
band3.scalebar.set_label('6kpc')
band3.scalebar.set_color('black')


#band6
band6 = ap.FITSFigure(band6_map, figure = fig, subplot = [0,-1,0.75,1])
band6.tick_labels.set_xformat('hh:mm:ss')
band6.tick_labels.set_yformat('dd:mm:ss')
band6.hide_xaxis_label()
band6.hide_xtick_labels()
band6.hide_yaxis_label()
band6.hide_ytick_labels()
band6.recenter(RA_c, DEC_c,width = 5./3600, height= 5./3600)
band6.show_colorscale(cmap='jet',vmax=0.408967,vmin=0.05,interpolation='nearest')
band6.add_colorbar(axis_label_text=r'S$_{N[II]}$ (Jy beam$^{-1}$ km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

band6.add_beam()
band6.beam.set_edgecolor('black')
band6.beam.set_facecolor('none')
band6.beam.set_hatch('//')
band6.add_scalebar(1.0/3600,linewidth=3)
band6.scalebar.set_label('6kpc')
band6.scalebar.set_color('black')


#band7
band7 = ap.FITSFigure(mom0_map, figure = fig, subplot = [0,-2,0.75,1])
band7.tick_labels.set_xformat('hh:mm:ss')
band7.tick_labels.set_yformat('dd:mm:ss')
band7.hide_xaxis_label()
band7.hide_xtick_labels()
band7.hide_yaxis_label()
band7.hide_ytick_labels()
band7.recenter(RA_c, DEC_c,width = 5./3600, height= 5./3600)
band7.show_colorscale(cmap='jet',vmax=13,vmin=-1,interpolation='nearest')
band7.add_colorbar(axis_label_text=r'S$_{C[II]}$ (Jy beam$^{-1}$ km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

band7.add_beam()
band7.beam.set_edgecolor('black')
band7.beam.set_facecolor('none')
band7.beam.set_hatch('//')
band7.add_scalebar(1.0/3600,linewidth=3)
band7.scalebar.set_label('6kpc')
band7.scalebar.set_color('black')


band3.save('spt0348_bands.pdf',dpi=250)
#pl.show()
pl.close()
band3.close()
band6.close()
band7.close()



#PDR

cp_fir = '../../../../Research/SPT0348/Data/PDR/Models/cp_fir.fits'
cp_co54 = '../../../../Research/SPT0348/Data/PDR/Models/cp_co54.fits'

#set up figure
fig = pl.figure()
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')


#C+ to L_IR
PDR = ap.FITSFigure(cp_fir, figure = fig)
# PDR.tick_labels.set_xformat('hh:mm:ss')
# PDR.tick_labels.set_yformat('dd:mm:ss')
# PDR.hide_xaxis_label()
# PDR.hide_xtick_labels()
# PDR.hide_yaxis_label()
# PDR.hide_ytick_labels()
pl.xlabel('log n (cm$^{-3}$)')
pl.ylabel('log G$_0$')
PDR.show_colorscale(cmap='jet',vmax=0.025,vmin=5e-6,interpolation='nearest', stretch='log')
PDR.add_colorbar(axis_label_text=r'C[II]/L$_{IR}$',axis_label_rotation=-90,axis_label_pad=50)
PDR.colorbar.set_ticks(np.array(['1e-1','1e-2','1e-3','1e-4','1e-5','1e-6']).astype(np.float))

PDR.show_contour(cp_co54,levels = [1,1e1,1e2,1e3,1e4,1e5,5e5], colors = 'w',linewidths = 2, linestyles='--')
PDR.show_contour(cp_co54,levels = [23,40], colors = 'w',linewidths = 2,linestyles='-')
PDR.show_contour(cp_fir,levels = [0.7e-4,2.7e-4], colors = 'm',linewidths = 2,linestyles='-')


PDR.save('cp_fir.pdf', dpi=250)
#pl.show()
pl.close()
PDR.close()




laboca = '../../../../../Research/SPT0348/Data/LABOCA/SPT0348-62-laboca-flux.fits'

#set up figure
fig = pl.figure()
pl.rc('font',size=22)
pl.rc('mathtext', default='regular')

RA_c = co.convHMS('03:48:42.271') #3:48:42.002
DEC_c = co.convDMS('-62:20:50.85') #-62:20:51.00

#C+ to L_IR
laboca_map = ap.FITSFigure(laboca, figure = fig)
laboca_map.tick_labels.set_xformat('hh:mm:ss')
laboca_map.tick_labels.set_yformat('dd:mm:ss')
laboca_map.hide_xaxis_label()
laboca_map.hide_xtick_labels()
laboca_map.hide_yaxis_label()
laboca_map.hide_ytick_labels()
laboca_map.show_colorscale(cmap='jet',vmax=0.05,vmin=-0.03,interpolation='nearest')
laboca_map.add_colorbar(axis_label_text=r'S$_{870\mu m}$ (Jy beam$^{-1}$)',axis_label_rotation=-90,axis_label_pad=30)

laboca_map.recenter(RA_c, DEC_c+2./3600,width = 1.2/60, height= 1.2/60)

band3_map = '../moment_maps/mom0s/spt0348_band3_spw0_clean1000_contsub_spw0_2sig.image.mom0_2sig.fits'
band6_map = '../moment_maps/mom0s/spt0348_band6_spw1_clean1000_contsub_bin5.image.mom0_2sig.fits'
band7_map = '../moment_maps/mom0s/spt0348_band7_spw1_clean1000_contsub_2sig.image.mom0.fits'
spire500 = '../../../../../Research/SPT0348/Data/SPIRE/SPT0348-62_spire_500umImG2.fits'
spire350 = '../../../../../Research/SPT0348/Data/SPIRE/SPT0348-62_spire_350umImG2.fits'
spire250 = '../../../../../Research/SPT0348/Data/SPIRE/SPT0348-62_spire_250umImG2.fits'

# laboca_map.show_contour(band3_map,levels = [0.2], colors = 'w',linewidths = 2, linestyles='-')
#laboca_map.show_contour(band6_map,levels = [0.1], colors = 'k',linewidths = 2, linestyles='-')
laboca_map.show_contour(band7_map,levels = [1], colors = 'w',linewidths = 2, linestyles='-')
# laboca_map.show_contour(spire500,levels = [0.02,0.028], colors = 'w',linewidths = 2, linestyles='-')
# laboca_map.show_contour(spire350,levels = [0.015,0.02], colors = 'k',linewidths = 2, linestyles='-')
# laboca_map.show_contour(spire250,levels = [0.01,0.015], colors = 'b',linewidths = 2, linestyles='-')

# laboca_map.add_beam()
# laboca_map.beam.set_edgecolor('black')
# laboca_map.beam.set_facecolor('none')
# laboca_map.beam.set_hatch('//')
laboca_map.add_scalebar(15.0/3600,linewidth=3)
laboca_map.scalebar.set_label('15-arcsec')
laboca_map.scalebar.set_color('black')


laboca_map.save('laboca.pdf', dpi=250)
#laboca_map.save('laboca_alma.pdf', dpi=250)
#laboca_map.save('laboca_spire.pdf', dpi=250)
pl.close()
laboca_map.close()



######################
# GalPak velocity maps
######################

#set up figure
fig = pl.figure()
pl.rc('font',size=30)
pl.rc('mathtext', default='regular')

# SPT0348-W
SPT_W = ap.FITSFigure('galpak_SPT0348_W_run6_gauss_obs_vel_map.fits', figure = fig, subplot = [0,0,0.75,1])
SPT_W.hide_xaxis_label()
SPT_W.hide_xtick_labels()
SPT_W.hide_yaxis_label()
SPT_W.hide_ytick_labels()
SPT_W.show_colorscale(cmap='jet',interpolation='nearest')
SPT_W.add_colorbar(axis_label_text=r'[CII] velocity (km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)

# SPT0348-W
SPT_E = ap.FITSFigure('galpak_SPT0348_E_run4_gauss_obs_vel_map.fits', figure = fig, subplot = [0,-1,0.75,1])
SPT_E.hide_xaxis_label()
SPT_E.hide_xtick_labels()
SPT_E.hide_yaxis_label()
SPT_E.hide_ytick_labels()
SPT_E.show_colorscale(cmap='jet',interpolation='nearest')
SPT_E.add_colorbar(axis_label_text=r'[CII] velocity (km s$^{-1}$)',axis_label_rotation=-90,axis_label_pad=50)


SPT_W.save('galpak_SPT0348_obs_vel.pdf', dpi=250)
pl.close()
SPT_W.close()
SPT_E.close()
