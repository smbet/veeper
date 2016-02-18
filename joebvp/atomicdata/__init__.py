from .. import joebgoodies as jbg
import numpy as np
import os

modpath=os.path.dirname(__file__)

vernerlist=np.genfromtxt(modpath+'/verner6.txt',dtype=None,delimiter=[10,8,3,4,3,2,9,6])
vernlam=jbg.arrfromcol(vernerlist,0)
vernion=jbg.arrfromcol(vernerlist,1)
vernzatom=jbg.arrfromcol(vernerlist,2)
vernnume=jbg.arrfromcol(vernerlist,3)
verngl=jbg.arrfromcol(vernerlist,4)
verngu=jbg.arrfromcol(vernerlist,5)
vernosc=jbg.arrfromcol(vernerlist,6)
vernp=jbg.arrfromcol(vernerlist,7)

for i in range(len(vernion)):
    vernion[i]=vernion[i].strip()

def lam2ion(restwave):
    if (isinstance(restwave,int))|(isinstance(restwave,float)):
        return vernion[jbg.closest(vernlam,restwave)].strip()
    else:
        ions=[]
        for rw in restwave: ions.append(vernion[jbg.closest(vernlam,rw)].strip())
        return ions

def lam2osc(restwave):
    return round(vernosc[jbg.closest(vernlam,restwave)],3)

def lam2vernp(restwave):
    if (isinstance(restwave,int)|isinstance(restwave,float)):
        return round(vernp[jbg.closest(vernlam,restwave)],2)
    else:
        return vernp[jbg.closest(vernlam,restwave)].round(2)

def ion2lam(ion):
    return vernlam[np.where(vernion==ion)[0]]

def ion2laminrange(ion,wave1,wave2,z=0.,frame='obs',pthresh=9.5):
    '''
    Usage: ion2laminrange(ion,wave1,wave2,z=0.,frame='obs',pthresh=9.5)
    Given ion (string) and wavelength range observed, return lambdas to search.
    Use frame='rest' to return restframe wavelengths even if supplying a redshift.
    '''
    try: lamidx=np.where(vernion==ion)[0]
    except:
        splitname=ion[0:2]+' '+ion[2:]
        try: lamidx=np.where(vernion==ion)[0]
        except:
            print 'Ion name doesn\'t match database.'
        else:
            restlams=vernlam[lamidx]
            obswaves=restlams*(1.+z)
            linesinrange=lamidx[np.where(jbg.between(obswaves,wave1,wave2)&(lam2vernp(restlams)>pthresh))]
            if frame=='obs':
               return vernlam[linesinrange]*(1.+z)
            if frame=='rest':
               return vernlam[linesinrange]
    else:
        restlams=vernlam[lamidx]
        obswaves=restlams*(1.+z)
        linesinrange=lamidx[np.where(jbg.between(obswaves,wave1,wave2)&(lam2vernp(restlams)>pthresh))]
        if frame=='obs':
           return vernlam[linesinrange]*(1.+z)
        if frame=='rest':
           return vernlam[linesinrange]
