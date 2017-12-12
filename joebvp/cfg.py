import os
import numpy as np
from astropy.constants import c as C

wave=0
flux=0
sigup=0
normwave=0
normflux=0
normsig=0
fitidx=0
fitcoeff=0
fitcovmtx=0
fiterrors=0
origpars=0
origparinfo=0
filename=''
field=''
homedir=os.path.expanduser('~')
spectrum = []
bad_pixels = []  # this will store bad pixel indices
min_sn = 0.2  # float, minimum signal to noise to include pixels in fit; leave as 0 and nothing will happen
# spectral_gaps = [[0,1174], [1300,1307], [1312.5,1330], [1560,1575], [1750,2000]]  # these gaps can be used to define bad pixels



Todd = False
if Todd:
    lsf='COS_LP1'
    instr=['COS','COS','COS','COS','STIS']
    gratings=['G130M','G160M','G185M','G225M','E230M']
    slits=['NA','NA','NA','NA','0.2x0.2']
    lsfranges=np.array([[1100,1400],[1400,1800],[1800,2100],[2100,2278],[2278,4100]])
    lps=['1','1','1','1','1']
    cen_wave=['1327','1600','1953','2250','2707']

else:  # this is for the casual user
    lsf='COS_LP3'
    instr=['COS','COS']
    gratings=['G130M','G160M']
    slits=['NA','NA']
    lsfranges=np.array([[1160,1450],[1450,1800]])
    lps=['3','3']
    cen_wave=['1318','1600']
    # spectral gaps for defyning bad pixels
    # spectral_gaps = [[0,1385], [1387,1387.6], [1558.4,1575.4], [1577.77,1578.3], [1748,2000]]  # BAL
    spectral_gaps = [[0,1174], [1300,1307], [1312.5,1330], [1560,1575], [1750,2000]]  # J0321



# fundamental constants
echarge = 4.803204505713468e-10
m_e = 9.10938291e-28
c = C.to('cm/s').value
c2 = 29979245800.0

# store LSFs and FGs
lsfs=[]
fgs=[]


outputdir = './'

VPparoutfile = outputdir + field + '_VP.dat'
VPmodeloutfile = outputdir + field + 'VPmodel.fits'
contoutfile = outputdir + 'continua.dat'
largeVPparfile = outputdir + '_VP_log.dat'
defaultcol = 13.1
defaultb = 20.0
defaultvlim = 1000.
lowblim = 4.
upperblim = 85.
upperblim_HI = 210.

# Plotting
ylim = (-0.2, 1.4)
xtick_fontsize = 'small'
ytick_fontsize = 'small'
xy_fontsize = 'small'
x_labelpad = 0
y_labelpad = -3
label_ypos = 0.2 * (ylim[0]+ylim[1])
label_fontsize = 8


# Handy definitions
lyseries=np.array([913.8260,914.0390,914.2860,914.5760,914.919,915.329,915.824,916.429,917.1806,918.1294,919.3514,920.9631,923.1504,926.2257,930.7483,937.8035,949.7431,972.5368,1025.7223,1215.6701])
