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
import joebvp.atomicdata as atomicdata
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

def cosvoigt_cont(vwave,cont,vpars):
	pars,info=joebvpfit.initlinepars(vpars[3],vpars[0],vpars,initinfo=[0,0,0])
	cfg.fitidx=joebvpfit.fitpix(vwave,pars)
	cfg.wave=vwave
	vflux=np.zeros(len(vwave))+1.
	factor=voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
	convfactor=convolvecos(vwave,factor*cont,vpars[0],vpars[3])
	vflux*=convfactor
	return vflux

def voigt(waves,line,coldens,bval,z,vels):
	tautot=np.zeros(len(waves))
	if len(cfg.lams)==0:
		lam,fosc,gam=atomicdata.setatomicdata(line)
		cfg.lams=lam ; cfg.fosc=fosc ; cfg.gam=gam
	for i in range(len(coldens)):
		#thatfactor=(1.+z[i])*(1.+vels[i]/c)
		thatfactor=(1.+z[i])

		lam0=cfg.lams[i]
		gam=cfg.gam[i]
		fosc=cfg.fosc[i]
		lam = waves/thatfactor
		dlam=bval[i]*lam0/c  #Doppler param in wavelength
		x=(lam-lam0-lam0*vels[i]/c)/dlam
		a=gam/(4.*np.pi*(c*1e13/lam0**2*dlam))
		vp=Hfunc(x,a)
		tauval=np.sqrt(np.pi)*cfg.echarge**2/cfg.m_e/cfg.c**2*(lam0*1e-8)**2/dlam*1e8*(10**coldens[i])*fosc*vp
		tautot+=tauval
	return np.exp(-tautot)

def get_lsfs():

	lsfobjs=[]
	for i,inst in enumerate(cfg.instr):
		lsfobjs.append(LSF(dict(name=inst,grating=cfg.gratings[i],
								life_position=cfg.lps[i], cen_wave=cfg.cen_wave[i],
						   		slit=cfg.slits[i])))
	cfg.lsfs=[]
	for fg in cfg.fgs:
		if isinstance(fg,int):
			lamobs=cfg.wave[fg]
			lsfmatch = jbg.wherebetween(lamobs, cfg.lsfranges[:, 0], cfg.lsfranges[:, 1])
			lsf = lsfobjs[lsfmatch[0]].interpolate_to_wv_array(cfg.wave[cfg.fgs] * u.AA,kind='cubic')
			cfg.lsfs.append(lsf['kernel'])
			break
		else:
			lamobs=np.median(cfg.wave[fg])
			lsfmatch = jbg.wherebetween(lamobs, cfg.lsfranges[:, 0], cfg.lsfranges[:, 1])
			lsf = lsfobjs[lsfmatch[0]].interpolate_to_wv_array(cfg.wave[fg] * u.AA,kind='cubic')
			cfg.lsfs.append(lsf['kernel'])

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
			convprof[ll] = convolve(paddedprof, cfg.lsfs[i], mode='same')[lsfwidth:-lsfwidth]


	return convprof