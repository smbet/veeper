import numpy as np
#from angles import *
#from velcorr import *
import sys
import subprocess

### Return index of array element closest to value
def closest(arr,value):
    if (isinstance(value,int)|isinstance(value,float)):
        idx = (np.abs(arr-value)).argmin()
    else:
        idx=[]
        for val in value: idx.append((np.abs(arr-val)).argmin())
    return idx

### Transform wavelength into velocity space centered on some line
def veltrans(redshift,waves,line):
    c=299792.458
    if (isinstance(line,int)|isinstance(line,float)):
        transline =c*(waves-(1+redshift)*line)/line/(1+redshift)
    else:
        transline=[]
        for ll in line: transline.append(c*(waves-(1+redshift)*ll)/ll/(1+redshift))
    return transline 

### Return row of pipe-delimited data
def pipedelimrow(data):
    numcols=len(data)
    rowstring=str(data[0])
    for i in range(1,numcols):
        rowstring=rowstring+'|'+str(data[i])
    return rowstring +' \n'

### Quickly load pipe-delimited data file
def loadpipe(filename,names=False):
    if names: tablearr=np.genfromtxt(filename,names=True,dtype=None,delimiter='|')
    else: tablearr=np.genfromtxt(filename,dtype=None,delimiter='|')
    return tablearr

### Return row of comma-delimited data
def commadelimrow(data):
    numcols=len(data)
    rowstring=str(data[0])
    for i in range(1,numcols):
        rowstring=rowstring+','+str(data[i])
    return rowstring +' \n'

### Quickly load comma-delimited data file
def loadcomma(filename,names=False):
    if names: tablearr=np.genfromtxt(filename,names=True,dtype=None,delimiter=',')
    else: tablearr=np.genfromtxt(filename,dtype=None,delimiter=',')
    return tablearr

### Return row of tab-delimited data
def tabdelimrow(data):
    numcols=len(data)
    rowstring=str(data[0])
    for i in range(1,numcols):
        rowstring=rowstring+'\t '+str(data[i])
    return rowstring +' \n'

### Make array from table column
def arrfromcol(datatable,col):
    arr=np.array([row[col] for row in datatable])
    return arr

### Return fraction of L_star for abs. mag.
def llstar(mag,lstar=-21.12):
    return 10.**((lstar-mag)/2.5)



### Convert hhmmss to degrees
def hhmmssdeg(coord):
    if coord[2]=='h':
        hours=float(coord.partition('h')[0])
        mins=float(coord.partition('h')[2].partition('m')[0])
        secs=float(coord.partition('h')[2].partition('m')[2].partition('s')[0])
    elif coord[2]==' ':
        parted=coord.partition(' ')
        hours=float(parted[0])
        mins=float(parted[2].partition(' ')[0])
        secs=float(parted[2].partition(' ')[2])
    elif coord[2]==':':
        parted=coord.partition(':')
        hours=float(parted[0])
        mins=float(parted[2].partition(':')[0])
        secs=float(parted[2].partition(':')[2])
    if coord[0]=='-':
        degs=-1.*(hours*15.+15.*mins/60.+15*secs/3600.)
	print 'RED FLAG! You shouldn\'t have negative hours!'	    
    else: degs=hours*15.+15.*mins/60.+15.*secs/3600.
    return degs

### Convert ddmmss to decimal degrees
def ddmmssdeg(coord):
    if ((coord[2]=='d') | (coord[3]=='d')):
	degrees=np.abs(float(coord.partition('d')[0]))
        mins=float(coord.partition('d')[2].partition('m')[0])
        secs=float(coord.partition('d')[2].partition('m')[2].partition('s')[0])
    elif ((coord[2]==' ') | (coord[3]==' ')):
        parted=coord.partition(' ')
	degrees=abs(float(parted[0]))
        mins=float(parted[2].partition(' ')[0])
        secs=float(parted[2].partition(' ')[2])
    elif ((coord[2]==':') | (coord[3]==':')):
        parted=coord.partition(':')
        degrees=np.abs(float(parted[0]))
        mins=float(parted[2].partition(':')[0])
        secs=float(parted[2].partition(':')[2])
    degs=degrees+mins/60.+secs/3600.
    if (coord[0]=='-'): degs=-1.*degs
    return degs

### Round to n significant digits
def roundto(x,n,retstring=1):
    if np.isnan(x)==True: return '-'
    if (x==0): number=0
    else: number=round(x, -int(np.floor(np.log10(abs(x))))+(n-1))
    numstr=str(number)
    if retstring==1: 
	if ((len(numstr)-1)!=n) & (numstr[0]!='0'):
		leftdig=len(numstr.partition('.')[0])
		rightdig=len(numstr.partition('.')[2])
		zerostoadd=n-leftdig-rightdig
		numstr=numstr+'0'*zerostoadd
	return numstr
    else: return number;print 'yo'

### Round number to n significant figures
def decimalplaces(x,n):
    #numstr=str(round(x,n))
    numstr="{:f}".format(float(round(x,n)))
    lhs=numstr.split('.')[0]
    rhs=numstr.split('.')[1]
    if len(str(rhs))<n: rhs=rhs+'0'*(n-len(str(rhs)))
    return lhs+'.'+rhs
    
### Volume of a sphere
def volsphere(r):
    return 4./3.*np.pi*r**3

### Relativistic velocity separation between z1,z2
def velsep(z1,z2):
    c=2.99792458e5
    #return c*(((1.+z1)**2-1.)/((1.+z1)**2+1.)-(((z2+1)**2-1)/((z2+1)**2+1)))
    return c*((1.+z1)**2-(1.+z2)**2)/((1.+z1)**2+(1.+z2)**2)

### Combined uncertainty of velocity separation between z1,z2
def sigvelsep(z1,z2,sigz1,sigz2):
    c=2.99792458e5
    difz1=2.*c*(1.+z1)*((1.+z1)**2+(1.+z2)**2)**-1 + -2.*c*(1.+z1)*((1.+z1)**2+(1.+z2)**2)**-2
    difz2=-2.*c*(1.+z2)*((1.+z1)**2+(1.+z2)**2)**-1 + -2.*c*(1.+z2)*((1.+z1)**2+(1.+z2)**2)**-2
    return np.sqrt((difz1*sigz1)**2+(difz2*sigz2)**2)


# '''
# ### Angular separation on the sky between two objects
# def angsep(ra1,dec1,ra2,dec2,units='degrees'):
#     '''
#     Takes two coordinates and returns the angular distance in degrees.
#     '''
#     if units=='degrees':
#         r1=Angle(ra1)
#         d1=Angle(dec1)
#         r2=Angle(ra2)
#         d2=Angle(dec2)
#         pos1=AngularPosition(alpha=r1.h,delta=d1.d)
#         pos2=AngularPosition(alpha=r2.h,delta=d2.d)
#         return r2d(pos1.sep(pos2))
#
# ### Return apparent mag for L_star fraction at given z
# def appmaglstar(ra,dec,llstar,z,lstar=-21.12):
#    absmag=lstar-2.5*np.log10(llstar)
#    distmod=mould_distance(ra,dec,z)['distmod']
#    return absmag+distmod
#
# ### Return L_star fraction at given z for given apparent mag
# def lstarappmag(appmag,z,lstar=-21.12):
#     dm=distmod(z)
#     absmag=appmag-dm
#     return 10.**((lstar-absmag)/2.5)
#
#
# ### Take an SDSS coord-derived name and return coords
# def SDSStodeg(sdssname):
#     if sdssname[0]=='J':
#         trimj=sdssname[1:]
#         if len(trimj.partition('+'))==3:
#             rapart=trimj.partition('+')[0]
#             raleftdec=rapart.partition('.')[0]
#             rarightdec=rapart.partition('.')[2]
#             if len(raleftdec)==6:
#                 ra=hhmmssdeg(rapart[:2]+'h'+rapart[2:4]+'m'+rapart[4:]+'s')
#             else:
#                 ra=hhmmssdeg('0'+rapart[:1]+'h'+rapart[1:3]+'m'+rapart[3:]+'s')
#             decpart=trimj.partition('+')[2]
#             dec=float(decpart[:2])+float(decpart[2:4])/60.+float(decpart[4:])/3600.
#         else:
#             rapart=trimj.partition('-')[0]
#             raleftdec=rapart.partition('.')[0]
#             rarightdec=rapart.partition('.')[2]
#             if len(raleftdec)==6:
#                 ra=hhmmssdeg(rapart[:2]+'h'+rapart[2:4]+'m'+rapart[4:]+'s')
#             else: ra=hhmmssdeg('0'+rapart[:1]+'h'+rapart[1:3]+'m'+rapart[3:]+'s')
#             decpart=trimj.partition('-')[2]
#             dec=-1*(float(decpart[:2])+float(decpart[2:4])/60.+float(decpart[4:])/3600.)
#         return ra,dec
#     else:
#         print 'Input must be string with format "JHHMMSS.ss+DDMMSS.ss" or "JHMMSS.ss+DDMMSS.ss"'
# '''

def SDSSnav(ra,dec):
    pre='http://skyserver.sdss3.org/dr9/en/tools/chart/navi.asp?ra='
    rastr=str(ra)
    decstr=str(dec)
    totstr=pre+rastr+'&dec='+decstr
    if 'linux' in sys.platform:  proc=subprocess.Popen(['xdg-open',totstr])
    else: proc=subprocess.Popen(['open',totstr])

def changefocus(window):
    '''
    Change focus between windows.  Accepts: 'terminal' and 'figure'
    '''
    if window=='terminal':
        if 'linux' in sys.platform:
            subprocess.call(['wmctrl','-a','cardinal'])
        else: 
            subprocess.call(['open','-a','terminal'])    
    if window=='figure':
        if 'linux' in sys.platform:
            subprocess.call(['wmctrl','-a','Figure'])
        else: 
            subprocess.call(['open','-a','XQuartz'])

### Test whether a value is two others
def between(testval,lim1,lim2):
	return ((testval>lim1)&(testval<lim2))

'''
### Calculate impact parameter (in kpc) corrected for peculiar velocities
def impact(qra,qdec,galra,galdec,galz):
    ang=angsep(qra,qdec,galra,galdec)*3600.
    scale=kpcarcsec(galz)
    return ang*scale    

### Calculate impact parameter (in kpc) corrected for peculiar velocities
def impactvelcorr(qra,qdec,galra,galdec,galz):
    ang=angsep(qra,qdec,galra,galdec)*3600.
    scale=mould_distance(galra,galdec,galz)['kpcarcsec']
    return ang*scale
'''


def wherebetween(testval,lim1,lim2):
	'''
	Given multiple ranges of values (lim1 and lim2 are same-sized),
	return indices where testval falls between lim1[idx] and lim2[idx].
	'''
	lim1=np.array(lim1);lim2=np.array(lim2)
	result=np.where((testval>lim1)&(testval<lim2))[0]
	return result

# '''
# ### Return sorted indices of matched coordinates from a list
# def coordmatch(ralist,declist,ra,dec,tol=0.0001):
#     radiff=np.abs(ralist-ra)
#     decdiff=np.abs(declist-dec)
#     match=where((radiff<=tol)&(decdiff<=tol))[0]
#     if len(match)>1:
#         seps=np.zeros(len(match))
#         for i in range(len(match)):
#             seps[i]=angsep(ralist[match[i]],declist[match[i]],ra,dec)
#         order=np.argsort(seps)
#         return match[order]
#     else:
#         return np.where((radiff<=tol)&(decdiff<=tol))[0]
#
# ### Estimate total H I mass given flat line profile
# def himass(dist,S,dv):
# 	'''
# 	Units: [dist] = Mpc, [S] = Jy, [dv] = km/s
# 	'''
# 	return 2.356e5*dist**2*S*dv
# '''

### Generate appropriate calls for a grid of subplots
def subplotgrid(num):
  if num==1: calls=[[1,1,1]]
  elif num==2: calls=[[2,1,1],[2,1,2]]
  elif num==3: calls=[[3,1,1],[3,1,2],[3,1,3]]
  elif num==4: calls=[[4,1,1],[4,1,2],[4,1,3],[4,1,4]]
  elif num==5: calls=[[3,2,1],[3,2,3],[3,2,5],[3,2,2],[3,2,4]]
  elif num==6: calls=[[3,2,1],[3,2,3],[3,2,5],[3,2,2],[3,2,4],[3,2,6]]
  elif num==7: calls=[[4,2,1],[4,2,3],[4,2,5],[4,2,2],[4,2,4],[4,2,6],[4,2,7]]
  elif num==8: calls=[[4,2,1],[4,2,3],[4,2,5],[4,2,7],[4,2,2],[4,2,4],[4,2,6],[4,2,8]]
  elif num==9: calls=[[5,2,1],[5,2,3],[5,2,5],[5,2,7],[5,2,2],[5,2,4],[5,2,6],[5,2,8],[5,2,9]]
  elif num==10: calls=[[5,2,1],[5,2,3],[5,2,5],[5,2,7],[5,2,2],[5,2,4],[5,2,6],[5,2,8],[5,2,9],[5,2,10]]
  elif num==11: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3]]
  elif num==12: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6]]
  elif num==13: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6],[5,3,9]]
  elif num==14: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6],[5,3,9],[5,3,12]]
  elif num==15: calls=[[5,3,1],[5,3,4],[5,3,7],[5,3,10],[5,3,13],[5,3,2],[5,3,5],[5,3,8],[5,3,11],[5,3,14],[5,3,3],[5,3,6],[5,3,9],[5,3,12],[5,3,15]]
  return calls

def subplotgridspec(num):
    if num==1: calls=[[0]]
    elif num==2: calls=[[0],[1]]
    elif num==3: calls=[[0],[1],[2]]
    elif num==4: calls=[[0],[1],[2],[3]]
    elif num==5: calls=[[0,0],[0,1],[0,2],[1,0],[1,1]]
    elif num==6: calls=[[0,0],[1,0],[2,0],[0,1],[1,1],[2,1]]
    elif num==7: calls=[[0,0],[1,0],[2,0],[3,0],[0,1],[1,1],[2,1]]
    elif num==8: calls=[[0,0],[1,0],[2,0],[3,0],[0,1],[1,1],[2,1],[3,1]]
    return calls

# '''
# ### Calculate Wilson confidence intervals
# def confinterval(hits,total,sig=1):
# 	'''
# 	Calculate Wilson confidence intervals
# 	Usage: confintervals(hits,total,sig=1)
# 	'sig' value corresponds to confidence level in multiples of sigma
# 	'''
# 	from astropy.stats.funcs import binom_conf_interval as conf
# 	if sig==1:
# 		uplim,lowlim=conf(hits,total)
# 	else:
# 		if sig==2: uplim,lowlim=conf(hits,total,conf=0.9545)
# 		elif sig==3: uplim,lowlim=conf(hits,total,cong=0.9973)
# 		else:
# 			print 'Invalid significance level'
# 			return
# 	return float(hits)/total,uplim,lowlim
#
# ### Make apparent column density profile
# def ACDprofile(wave,flux,err,cont,z,restlam,velrange):
#     from atomicdata import *
#     osc=lam2osc(restlam)
#     vel=veltrans(z,wave,restlam)
#     sigsqN=0.
#     sigsqEW=0.
#     velidx=np.where((vel>=-velrange)&(vel<=velrange))[0]
#     EW=0.
#     N=0.
#     sigsqN=0.
#     sigsqEW=0.
#     ACD=np.zeros(len(velidx))
#     for ii in range(len(velidx)):
#         i=velidx[ii]
#         EW=EW+(wave[i+1]-wave[i])*((cont[i]-flux[i])/cont[i])/(1.+z)
#         sigsqEW=sigsqEW+((wave[i+1]-wave[i])/cont[i]/(1.+z)*err[i])**2
#         tauv=np.log(cont[i]/(flux[i]))
#         tauverr=np.log(cont[i]/(cont[i]-err[i]))
#         ACD[ii]=N+1./2.654e-15/restlam/osc*(vel[i+1]-vel[i])*tauv
#         sigsqN=sigsqN+(1./2.654e-15/restlam/osc*(vel[i+1]-vel[i])*tauverr)**2
#     return ACD
# '''