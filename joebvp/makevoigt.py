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
from astropy import constants as const
import sys,os
import cfg
from sklearn.cluster import MeanShift, estimate_bandwidth


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

atomdata=np.genfromtxt(modpath+'/atomicdata/LineListUVgam.dat',dtype=None,delimiter='|')
atomlam=jbg.arrfromcol(atomdata,0)
atomtrans=jbg.arrfromcol(atomdata,1)
atomosc=jbg.arrfromcol(atomdata,2)
atomgam=jbg.arrfromcol(atomdata,3)

def osc(line):
	match=jbg.closest(atomlam,line)
	return atomosc[match]

def gamma(line):
	match=jbg.closest(atomlam,line)
	return atomgam[match]

def restwave(line):
	match=jbg.closest(atomlam,line)
	return atomlam[match]

def Hfunc(x,a):
	z=x+1j*a
	I=wofz(z).real
	return I


def cosvoigt(vwave,vpars):
	vflux=np.zeros(len(vwave))+1.
	factor=voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
	convfactor=convolvecos(vwave,factor,vpars[0],vpars[3])
	vflux*=convfactor
	return vflux

def voigt(waves,line,coldens,bval,z,vels):
	tau=[]
	lam=[]
	for i in range(len(coldens)):
		#thatfactor=(1.+z[i])*(1.+vels[i]/c)
		thatfactor=(1.+z[i])
		lam0=restwave(line[i])
		gam=gamma(line[i])
		fosc=osc(line[i])
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


def convolvecos(wave,profile,lines,zs):
	fitwaves=wave
	if cfg.wavegroups==[]:
		X = np.array(zip(fitwaves,np.zeros(len(fitwaves))), dtype=float)
		ms = MeanShift(bandwidth=25.)
		ms.fit(X)
		cfg.wgidxs = ms.labels_
		cfg.wavegroups = ms.cluster_centers_
		cfg.uqwgidxs=np.unique(cfg.wgidxs)

	convprof=profile
	for ll in cfg.uqwgidxs[:-3]:
		matches=np.where(cfg.wgidxs==ll)[0]
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
		convprof[matches]=convolve(profile[matches],lsf,mode='same')
		convprof[matches[:16]]=[1.]*16 ; convprof[matches[-16:]]=[1.]*16

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