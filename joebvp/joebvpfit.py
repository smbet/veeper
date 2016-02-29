### Fit a voigt profile to a normalized spectrum
### Parameters: [lam_rest,column_density,doppler_param,z,velocity_offset]

import numpy as np
import joebgoodies as jbg
from stsci.tools import nmpfit
import makevoigt
import cfg

c=299792.458

def foldpars(pars,numpars=5):
	rows=len(pars)/numpars
	fpars=[]
	for i in range(numpars):
		fparrow=[]
		for j in range(rows):
			fparrow.append(pars[i+numpars*j])
		fpars.append(fparrow)
	return fpars

def unfoldpars(pars,numpars=5):
	rows=len(pars[0])
	ufpars=[]
	for j in range(rows):
		for i in range(numpars):
			ufpars.append(pars[i][j])
	return ufpars

def voigtfunc(vwave,vpars):
	vflux=np.zeros(len(vwave))+1.
	factor=makevoigt.voigt(vwave,vpars[0],vpars[1],vpars[2],vpars[3],vpars[4])
	convfactor=makevoigt.convolvecos(vwave,factor,vpars[0],vpars[3])
	vflux*=convfactor
	return vflux

def voigterrfunc(p,x,y,err,fjac=None):
	fp=foldpars(p)
	model=voigtfunc(x,fp)
	status=0
	return([status,(y-model)/err])

def fitpix(wave,pararr):
	ll=pararr[0]
	lz=pararr[3]
	lv1=pararr[5]
	lv2=pararr[6]
	relpix=[]
	for i in range(len(ll)):
		w1=ll[i]*(1.+lz[i]+lv1[i]/c)
		w2=ll[i]*(1.+lz[i]+lv2[i]/c)
		p1=jbg.closest(wave,w1)
		p2=jbg.closest(wave,w2)
		relpix.extend(range(p1-10,p2+10))
		rp=np.unique(np.array(relpix))
	return rp

def prepparinfo(linepars,parflags):
	parinfo=[]
	parflags=np.array(parflags)
	numpars=5
	for i in range(len(parflags[0])):
		parinfo.extend([{'fixed':1},{'fixed':0},{'fixed':0},{'fixed':1},{'fixed':0}])
		for j in range(1,len(parflags)):
			if parflags[j][i]==1: parinfo[i*numpars+j]['fixed']=1
			elif parflags[j][i]<=0: parinfo[i*numpars+j]['fixed']=0
			else:
				matches=np.where(np.array(parflags[j])==parflags[j][i])[0]
				if matches[0]!=i:
					tiedpar=int(matches[0]*numpars+j)
					parinfo[i*numpars+j]['tied']='p['+str(tiedpar)+']'
		col=round(linepars[1][i],2)
		vel=round(linepars[4][i],2)
		bpar=round(linepars[2][i],2)
		parinfo[numpars*i+1]['limited']=[1,1]
		parinfo[numpars*i+1]['limits']=[round(col-5.,2),round(col+5.,2)]
		parinfo[numpars*i+2]['limited']=[1,1]
		### adjust b-value limits based on species
		lydiff=abs(linepars[0][i]-cfg.lyseries)
		lymatch = np.where(abs(lydiff)<=0.05)[0]
		if lymatch:
			parinfo[numpars*i+2]['limits']=[max([6.,bpar-10.]),min([bpar+10],150.)]
		else:
			parinfo[numpars*i+2]['limits']=[max([6.,bpar-10.]),min([bpar+10],85.)]
		#else: parinfo[numpars*i+2]['limits']=[1.,150.]
	#parinfo[numpars*i+2]['maxstep']=5.
		parinfo[numpars*i+2]['step']=0.5
		parinfo[numpars*i+2]['mpside']=2
		#parinfo[numpars*i+2]['relstep']=0.0001
		parinfo[numpars*i+4]['limited']=[1,1]
		### Allow velocity to flop around
		if parflags[4][i]<0:
			flopamt=abs(parflags[4][i])
			parinfo[numpars*i+4]['limits']=[round(vel-flopamt,2),round(vel+flopamt,2)]
		elif len(linepars)>5:
			v1=round(linepars[5][i],2) ; v2=round(linepars[6][i],2)
			parinfo[numpars*i+4]['limits']=[v1,v2]
		else:
			parinfo[numpars*i+4]['limits']=[round(vel-50.,2),round(vel+50.,2)]
		parinfo[numpars*i+4]['step']=1.
		#parinfo[numpars*i+4]['relstep']=0.01
	return parinfo

def joebvpfit(wave,flux,sig,linepars,flags):

	xtol=1e-11
	gtol=1e-11
	# Only feed to the fitter the parameters that go into the model
	partofit=linepars[:5]
	parinfo=prepparinfo(partofit,flags)
	# Save the velocity windows to add back to the parameter array
	vlim1=linepars[5] ; vlim2=linepars[6]
	partofit=unfoldpars(partofit)
	modelvars={'x':wave,'y':flux,'err':sig}
	m=nmpfit.mpfit(voigterrfunc,partofit,functkw=modelvars,parinfo=parinfo,nprint=1,quiet=0,fastnorm=1,ftol=1e-10,xtol=xtol,gtol=gtol)
	if m.status <= 0: print 'Fitting error:',m.errmsg
	fitpars=foldpars(m.params)
	# Add velocity windows back to parameter array
	fitpars.append(vlim1) ; fitpars.append(vlim2)
	fiterrors=foldpars(m.perror)

	print '\nFit results: \n'
	for i in range(len(fitpars[0])):
		print jbg.tabdelimrow([round(fitpars[0][i],2),jbg.decimalplaces(fitpars[3][i],5),jbg.roundto(fitpars[1][i],5),jbg.roundto(fitpars[2][i],5),jbg.roundto(fitpars[4][i],5)])[:-2]
		print jbg.tabdelimrow([' ',' ',' ',round(fiterrors[1][i],3),round(fiterrors[2][i],3),round(fiterrors[4][i],3)])

	return fitpars,fiterrors


def initlinepars(zs,restwaves,initvals=[],initinfo=[]):

	### Look for mult\iplet membership of each line
	seriesassoc=np.zeros(len(restwaves))-99
	for i in range(len(restwaves)):
		for j in range(len(cfg.multiplets)):
			currmult=cfg.multiplets[j]
			if (abs(restwaves[i]-currmult[jbg.closest(currmult,restwaves[i])]) < 0.01):
				seriesassoc[i]=j

	initpars=[[],[],[],[],[]]
	parinfo=[[],[],[],[],[]]
	matches=[]
	defaultcol=cfg.defaultcol
	defaultb=cfg.defaultb
	if initvals==[]:
		for i in range(len(restwaves)):
			initpars[0].extend([restwaves[i]])
			initpars[1].extend([defaultcol])
			initpars[2].extend([defaultb])
			initpars[3].extend([zs[i]])
			initpars[4].extend([0.])
			initpars[5].extend([-cfg.defaultvlim])
			initpars[6].extend([cfg.defaultvlim])
	else:
		if len(initvals)==5:
			for i in range(len(restwaves)):
				initpars[0].extend([initvals[0][i]])
				initpars[1].extend([initvals[1][i]])
				initpars[2].extend([initvals[2][i]])
				initpars[3].extend([initvals[3][i]])
				initpars[4].extend([initvals[4][i]])
				initpars[5].extend([-cfg.defaultvlim])
				initpars[6].extend([cfg.defaultvlim])
		else:
			initpars=[[],[],[],[],[],[],[]]
			for i in range(len(restwaves)):
				initpars[0].extend([initvals[0][i]])
				initpars[1].extend([initvals[1][i]])
				initpars[2].extend([initvals[2][i]])
				initpars[3].extend([initvals[3][i]])
				initpars[4].extend([initvals[4][i]])
				initpars[5].extend([initvals[5][i]])
				initpars[6].extend([initvals[6][i]])

	parinfo=np.zeros([5,len(restwaves)],dtype=int)
	parinfo[0]=parinfo[0]+1
	parinfo[3]=parinfo[3]+1
	if ((initinfo==[])&(initvals==[])):
		uqions=np.unique(seriesassoc).tolist()

		if -99 in uqions: uqions.remove(-99)
		flagctr=2
		for uqion in uqions:
			ionmatch=np.where(seriesassoc == uqion)[0]
			uqzs=np.unique(zs[ionmatch])
			for uz in uqzs:
				matchcrit=(zs==uz)&(seriesassoc==uqion)
				rellines=np.where(matchcrit)[0]
				uqlams=np.unique(restwaves[rellines])
				if len(uqlams)>1:
					complist=[]
					numcomps=[]
					for ul in uqlams:
						matchcrit2=matchcrit&(restwaves==ul)
						matches=np.where(matchcrit2)[0]
						complist.append(matches)
						numcomps.append(len(matches))
					numcomps=np.array(numcomps)
					complist=np.array(complist)
					compidxsort=sorted(range(len(numcomps)),key = lambda x: numcomps[x],reverse=True)
					numcompsort=numcomps[compidxsort]
					complistsort=complist[compidxsort]
					maxcomps=numcompsort[0]
					for compidx in range(maxcomps):
						for li in range(len(complistsort)):
							if compidx<numcompsort[li]:
								parinfo[1][complistsort[li][compidx]]=flagctr
								parinfo[2][complistsort[li][compidx]]=flagctr
								parinfo[4][complistsort[li][compidx]]=flagctr
							else: continue
						flagctr+=1

	elif initinfo!=[]:
		parinfo[1]=initinfo[0]
		parinfo[2]=initinfo[1]
		parinfo[4]=initinfo[2]
	elif ((initinfo==[])&(initvals!=[])):
		### Fix measurements that are imported
		for i in range(len(restwaves)):
			if ((initpars[1][i]!=defaultcol)&(initpars[2][i]!=defaultb)):
				parinfo[1][i]=1 ; parinfo[2][i]=1 ; parinfo[4][i]=1
		uqions=np.unique(seriesassoc).tolist()
		if -99 in uqions: uqions.remove(-99)
		flagctr=2
		for uqion in uqions:
			ionmatch=np.where(seriesassoc == uqion)[0]
			uqzs=np.unique(zs[ionmatch])
			for uz in uqzs:
				matchcrit=(zs==uz)&(seriesassoc==uqion)
				rellines=np.where(matchcrit)[0]
				uqlams=np.unique(restwaves[rellines])
				if len(uqlams)>1:
					complist=[]
					numcomps=[]
					for ul in uqlams:
						matchcrit2=matchcrit&(restwaves==ul)
						matches=np.where(matchcrit2)[0]
						complist.append(matches.tolist())
						numcomps.append(len(matches))
					numcomps=np.array(numcomps)
					complist=np.array(complist)
					compidxsort=sorted(range(len(numcomps)),key = lambda x: numcomps[x])
					complistsort=complist[compidxsort]
					complistsort=complistsort.tolist()
					for i in range(len(complistsort)-1):
						for idx in complistsort[i]:
							parinfo[1][idx]=flagctr
							parinfo[2][idx]=flagctr
							parinfo[4][idx]=flagctr
							for j in range(i+1,len(complistsort)):
								idxidx=jbg.closest(initvals[4][complistsort[j]],initvals[4][idx])
								clidx=complistsort[j][idxidx]
								parinfo[1][clidx]=flagctr
								parinfo[2][clidx]=flagctr
								parinfo[4][clidx]=flagctr
								complistsort[j].remove(clidx)
							flagctr+=1
	return initpars,parinfo