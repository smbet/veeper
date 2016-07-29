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
from astropy.table import Table
import pdb
from linetools.spectra.lsf import LSF

ln2=np.log(2)
c=const.c.value/1e3

modpath=os.path.dirname(__file__)

if cfg.lsf=='COS_LP1':
	G130file=modpath+'/LSF/COS_LP1/G130M.REBIN.LSF.dat'
	relpix,w1150,w1200,w1250,w1300,w1350,w1400,w1450=np.loadtxt(G130file,unpack=True,skiprows=1)
	G160file=modpath+'/LSF/COS_LP1/G160M.REBIN.LSF.dat'
	relpix,w1450,w1500,w1550,w1600,w1650,w1700,w1750=np.loadtxt(G160file,unpack=True,skiprows=1)
else:
	raise('Error: Only COS LP1 is supported at this time')


def Hfunc(x,a):
	z=x+1j*a
	I=wofz(z).real
	return I


def cosvoigt(vwave,vpars):
	if isinstance(cfg.fitidx,int):
		pars,info=joebvpfit.initlinepars(vpars[3],vpars[0],vpars,initinfo=[0,0,0])
		cfg.fitidx=joebvpfit.fitpix(vwave,pars)

	vflux=np.zeros(len(vwave))+1.
	factor=voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
	convfactor=convolvecos(vwave,factor,vpars[0],vpars[3])
	vflux*=convfactor
	return vflux

def voigt(waves,line,coldens,bval,z,vels):
	tau=[]
	lam=[]
	#lams,fs,gs=getdata(line)
	for i in range(len(coldens)):
		#thatfactor=(1.+z[i])*(1.+vels[i]/c)
		thatfactor=(1.+z[i])
		#lam0=restwave(line[i])
		#gam=gamma(line[i])
		#fosc=osc(line[i])
		try:
			lam0=cfg.lams[i]
			gam=cfg.gam[i]
			fosc=cfg.fosc[i]
		except:
			lam0,gam,fosc=joebvpfit.setatomicdata(line)
		#thatfactor=1
		#lam.append(arange(lam0-5.,lam0+5.,step=.001))
		lam.append(waves/thatfactor)
		dlam=bval[i]*lam0/c  #Doppler param in wavelength
		x=(lam[i]-lam0-lam0*vels[i]/c)/dlam
		a=gam/(4.*np.pi*(c*1e13/lam0**2*dlam))
		vp=Hfunc(x,a)
		tauval=np.sqrt(np.pi)*cfg.echarge**2/cfg.m_e/cfg.c**2*(lam0*1e-8)**2/dlam*1e8*(10**coldens[i])*fosc*vp
		tau.append(tauval)
	taus=tau
	tautot=np.zeros(len(waves))
	for i in range(len(taus)):
		tautot+=taus[i]
	return np.exp(-tautot)

#TODO: retrieve lsf from linetools, compare with tabulated versions
def get_lsfs():
	cfg.lsfs=[]
	for ll in cfg.uqwgidxs:
		#if cfg.
		#lsfobj=LSF()
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
		dividers = np.where(df > buf)[0]
		if len(dividers)==0:
			fgs=cfg.fitidx
		else:
			fgs=[np.arange(0,dividers[0])]
			for i, idx in enumerate(dividers[:-1]):
				fgs.append(np.arange(idx,dividers[i+1]))
			fgs.append(np.arange(dividers[-1],len(cfg.fitidx)))
			# Check joebvpfit.fitpix line groups to see if any fall in separate wavegroups
			for i, gg in enumerate(fgs):
				if len(np.unique(cfg.wgidxs[gg]))!=1:
					domuq=scipy.stats.mode(cfg.wgidxs[gg])[0][0]
					tochange=np.where(cfg.wgidxs[gg]!=domuq)[0]
					cfg.wgidxs[gg[tochange]] = domuq
		get_lsfs()
	convprof=profile
	for i,ll in enumerate(cfg.uqwgidxs):
		matches=cfg.fitidx[np.where(cfg.wgidxs==ll)[0]]
		paddedprof=np.insert(profile[matches],0,[1.]*16)
		paddedprof=np.append(paddedprof,[1.]*16)
		#convprof[matches]=convolve(profile[matches],lsf,mode='same')[16:-16]
		convprof[matches] = convolve(paddedprof, cfg.lsfs[i], mode='same')[16:-16]
		#convprof[matches[:16]]=[1.]*16 ; convprof[matches[-16:]]=[1.]*16
	#pdb.set_trace()
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