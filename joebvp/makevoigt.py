#######################################################################
###
###  Main function for generating Voigt profiles.  
###  convenience function: cosvoigt(waves,params)
###  params format (each entry a list): [restwave,coldens,bval,z,vel]
###
#######################################################################


import joebgoodies as jbg
import numpy as np  
from scipy.signal import convolve
from scipy.special import wofz
import scipy.stats
from astropy import constants as const
import sys,os
import cfg
from sklearn.cluster import MeanShift, estimate_bandwidth
import joebvpfit
import astropy.units as u
from astropy.io import ascii
from astropy.table import Table
from linetools.spectra.lsf import LSF

ln2=np.log(2)
c=const.c.value/1e3

modpath=os.path.dirname(__file__)


# if cfg.lsf=='COS_LP1':
# 	G130file=modpath+'/LSF/COS_LP1/G130M.REBIN.LSF.dat'
# 	relpix,w1150,w1200,w1250,w1300,w1350,w1400,w1450=np.loadtxt(G130file,unpack=True,skiprows=1)
# 	G160file=modpath+'/LSF/COS_LP1/G160M.REBIN.LSF.dat'
# 	relpix,w1450,w1500,w1550,w1600,w1650,w1700,w1750=np.loadtxt(G160file,unpack=True,skiprows=1)
# else:
# 	raise('Error: Only COS LP1 is supported at this time')


def Hfunc(x,a):
	z=x+1j*a
	I=wofz(z).real
	return I


def cosvoigt(vwave,vpars):
	pars,info=joebvpfit.initlinepars(vpars[3],vpars[0],vpars,initinfo=[0,0,0])
	cfg.fitidx=joebvpfit.fitpix(vwave,pars)
	cfg.wave=vwave
	vflux=np.zeros(len(vwave))+1.
	factor=voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
	convfactor=convolvecos(vwave,factor,vpars[0],vpars[3])
	vflux*=convfactor
	return vflux

def voigt(waves,line,coldens,bval,z,vels):
	tautot=np.zeros(len(waves))
	if len(cfg.lams)==0:
		lam,fosc,gam=setatomicdata(restwaves)
		cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam
	for i in range(len(coldens)):
		#thatfactor=(1.+z[i])*(1.+vels[i]/c)
		thatfactor=(1.+z[i])

		lam0=cfg.lams[i]
		gam=cfg.gam[i]
		fosc=cfg.fosc[i]

		#thatfactor=1
		#lam.append(arange(lam0-5.,lam0+5.,step=.001))
		lam = waves/thatfactor
		dlam=bval[i]*lam0/c  #Doppler param in wavelength
		x=(lam-lam0-lam0*vels[i]/c)/dlam
		a=gam/(4.*np.pi*(c*1e13/lam0**2*dlam))
		vp=Hfunc(x,a)
		tauval=np.sqrt(np.pi)*cfg.echarge**2/cfg.m_e/cfg.c**2*(lam0*1e-8)**2/dlam*1e8*(10**coldens[i])*fosc*vp
		#tau.append(tauval)
		tautot+=tauval
	return np.exp(-tautot)

def get_lsfs():

	lsfobjs=[]
	for i,inst in enumerate(cfg.instr):
		lsfobjs.append(LSF(dict(name=inst,grating=cfg.gratings[i],life_position=cfg.lps[i], cen_wave=cfg.cen_wave[i])))
	cfg.lsfs=[]
	for fg in cfg.fgs:

		if isinstance(fg,int):
			lamobs=cfg.wave[fg]
			lsfmatch = jbg.wherebetween(lamobs, cfg.lsfranges[:, 0], cfg.lsfranges[:, 1])
			lsf = lsfobjs[lsfmatch].interpolate_to_wv_array(cfg.wave[cfg.fgs] * u.AA)
			cfg.lsfs.append(lsf['kernel'])
			break
		else:
			lamobs=np.median(cfg.wave[fg])
			lsfmatch = jbg.wherebetween(lamobs, cfg.lsfranges[:, 0], cfg.lsfranges[:, 1])
			lsf = lsfobjs[lsfmatch].interpolate_to_wv_array(cfg.wave[fg] * u.AA)
			cfg.lsfs.append(lsf['kernel'])
	'''
	for ll in cfg.uqwgidxs:

		lamobs=cfg.wavegroups[ll,0]
		if (lamobs<=1175): lsf=w1150
		elif (lamobs>1175) & (lamobs<=1225): lsf=w1200
		elif (lamobs>1225) & (lamobs<=1275): lsf=w1250
		elif (lamobs>1275) & (lamobs<=1325): lsf=w1300
		elif (lamobs>1325) & (lamobs<=1375): lsf=w1350
		elif (lamobs>1375) & (lamobs<=1425): lsf=w1400
		elif (lamobs>1425) & (lamobs<=1475): lsf=w1450
		elif (lamobs>1475) & (lamobs<=1525): lsf=w1500
		elif (lamobs>1525) & (lamobs<=1575): lsf=w1550
		elif (lamobs>1575) & (lamobs<=1625): lsf=w1600
		elif (lamobs>1625) & (lamobs<=1675): lsf=w1650
		elif (lamobs>1675) & (lamobs<=1725): lsf=w1700
		elif (lamobs>1725): lsf=w1750
		cfg.lsfs.append(lsf)
	'''

def convolvecos(wave,profile,lines,zs):
	if len(wave)>len(cfg.fitidx):
		fitwaves=wave[cfg.fitidx]
	else:
		fitwaves=wave
	if cfg.wavegroups==[]:
		X = np.array(zip(fitwaves,np.zeros(len(fitwaves))), dtype=float)
		ms = MeanShift(bandwidth=25.)
		ms.fit(X)
		cfg.wgidxs = ms.labels_
		cfg.wavegroups = ms.cluster_centers_
		cfg.uqwgidxs=np.unique(cfg.wgidxs)
		# Identify groups of fitidxs
		buf=4
		df=cfg.fitidx[1:]-cfg.fitidx[:-1]
		dividers = np.where(df > buf)[0] #These are the last indices of each group
		if len(dividers)==0:
			cfg.fgs=cfg.fitidx
		else:
			cfg.fgs=[np.arange(cfg.fitidx[0],cfg.fitidx[dividers[0]])] #1st group
			for i, idx in enumerate(dividers[:-1]):
				cfg.fgs.append(np.arange(cfg.fitidx[idx+1],cfg.fitidx[dividers[i+1]])) # 'i+1' b/c 1st group handled separately
			cfg.fgs.append(np.arange(cfg.fitidx[dividers[-1]+1],cfg.fitidx[-1]+1))  #last group
			newfgs=[]
			### Deal with wavegroups that may have huge jumps in wavelength between pixels
			for i,fg in enumerate(cfg.fgs):

				diffs=cfg.wave[fg[1:]]-cfg.wave[fg[:-1]]
				outlier=np.where(diffs>(30.*np.median(diffs)))[0]

				if len(outlier)>0:
					if len(outlier)==1:
						cfg.fgs[i]=np.arange(fg[0],fg[outlier+1])
						newfgs.append(np.arange(fg[outlier]+1,fg[-1]))
						if len(newfgs[-1])<15:
							nummorepix=15-len(newfgs)
							newfgs[-1]=np.concatenate([newfgs[-1],np.arange(fg[-1],fg[-1]+nummorepix+1)])  #add 1 to 2nd arg to make sure this val is included
					else:
						cfg.fgs[i]=np.arange(fg[0],fg[outlier[0]])
						for j,out in enumerate(outlier):
							if out==outlier[-1]:
								newfgs.append(np.arange(fg[out]+1,fg[-1]+1)) #add 1 to 2nd arg to make sure this val is included
							else:
								newfgs.append(np.arange(fg[out]+1,fg[outlier[j+1]]+1))   #add 1 to 2nd arg to make sure this val is included
							if len(newfgs[-1])<15:
								nummorepix=15-len(newfgs)
								newfgs[-1]=np.concatenate([newfgs[-1],np.arange(fg[-1],fg[-1]+nummorepix+1)]) #add 1 to 2nd arg to make sure this val is included

			for fg in newfgs:
				cfg.fgs.append(fg)
		get_lsfs()
	convprof=profile
	'''
	for i,ll in enumerate(cfg.uqwgidxs):
		matches=cfg.fitidx[np.where(cfg.wgidxs==ll)[0]]
		paddedprof=np.insert(profile[matches],0,[1.]*16)
		paddedprof=np.append(paddedprof,[1.]*16)
		#convprof[matches]=convolve(profile[matches],lsf,mode='same')[16:-16]
		convprof[matches] = convolve(paddedprof, cfg.lsfs[i], mode='same')[16:-16]
		#convprof[matches[:16]]=[1.]*16 ; convprof[matches[-16:]]=[1.]*16
	'''
	for i,ll in enumerate(cfg.fgs):
		if isinstance(ll,int):
			lsfwidth=len(cfg.fgs)/2+1
			paddedprof = np.insert(profile[cfg.fgs], 0, [1.] * lsfwidth)
			paddedprof = np.append(paddedprof, [1.] * lsfwidth)
			convprof[cfg.fgs] = convolve(paddedprof, cfg.lsfs[i], mode='same')[lsfwidth:-lsfwidth]
			break
		else:

			lsfwidth=len(ll)/2+1
			paddedprof = np.insert(profile[ll], 0, [1.] * lsfwidth)
			paddedprof=np.append(paddedprof,[1.]*lsfwidth)
			'''
			tab2write=Table([cfg.wave[ll]],names=['wave'])
			ascii.write(tab2write, output='wave_'+str(i)+'.dat')
			tab2write=Table([cfg.lsfs[i]],names=['lsf'])
			ascii.write(tab2write, output='lsf'+str(i)+'.dat')
			tab2write=Table([paddedprof],names=['paddedprof'])
			ascii.write(tab2write, output='padprof_'+str(i)+'.dat')
			tab2write=Table([cfg.normflux[ll]],names=['flux'])
			ascii.write(tab2write, output='flux_'+str(i)+'.dat')
			'''
			convprof[ll] = convolve(paddedprof, cfg.lsfs[i], mode='same')[lsfwidth:-lsfwidth]


	return convprof

'''
def convolvecos(profile,line,z):
    lamobs=restwave(line)*(1.+z)
    if (lamobs<=1175): lsf=w1150
    elif (lamobs>1175) 	& (lamobs<=1225): lsf=w1200
    elif (lamobs>1225) & (lamobs<=1275): lsf=w1250
    elif (lamobs>1275) & (lamobs<=1325): lsf=w1300
    elif (lamobs>1325) & (lamobs<=1375): lsf=w1350
    elif (lamobs>1375) & (lamobs<=1425): lsf=w1400
    elif (lamobs>1425) & (lamobs<=1475): lsf=w1450
    elif (lamobs>1475) & (lamobs<=1525): lsf=w1500
    elif (lamobs>1525) & (lamobs<=1575): lsf=w1550
    elif (lamobs>1575) & (lamobs<=1625): lsf=w1600
    elif (lamobs>1625) & (lamobs<=1675): lsf=w1650
    elif (lamobs>1675) & (lamobs<=1725): lsf=w1700
    elif (lamobs>1725): lsf=w1750
    return convolve(profile,lsf,mode='same')
    '''