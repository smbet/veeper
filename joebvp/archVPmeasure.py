import numpy as np
import numpy.polynomial.legendre as L
from numpy.linalg import svd
import matplotlib
import sys,os
if 'linux' in sys.platform: matplotlib.use('wxagg')  # and also try macosx and wxagg and tkagg
else: matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import pyfits as pf
#import datetime, time
#from astropy import constants as const
#from time import sleep
from joebvp import joebgoodies as jbg
#from joebvpfit import *
from joebvp import joebvpfit
from joebvp import cfg
from joebvp import atomicdata
from astropy.io import ascii
import pdb

plt.rcParams['font.family']='serif'

c=299792.458

def closewave(wavearr,wave):
    closeidx=jbg.closest(wavearr,wave)
    return wavearr[closeidx]


def onclick(event):
    global fittog,EWtog,veltog
    global fitidx
    global wave
    global flux
    global EWwave1,newvel,newrestwave,newz,newvel
    global zEW,EWrestwave
    if fittog==1:
        if event.button==1:
            fitidx=np.append(fitidx,jbg.closest(wave,event.xdata))
            fitidx.sort()
            updateplot()
        if event.button==2:
            closeidx=jbg.closest(wave[fitidx],event.xdata)
            fitidx=np.delete(fitidx,closeidx)
            updateplot()
    #if EWtog==1:
    #    if event.button==1:
    #        EWwave1=jbg.closest(wave,event.xdata)
    #        print "Right-click line right bound."
    #    if event.button==2:
    #        EWwave2=jbg.closest(wave,event.xdata)
    #        calcEW(EWwave1,EWwave2,zEW,EWrestwave)
    if veltog==1:
        if isinstance(newrestwave,float): nrw=newrestwave
        else: nrw=newrestwave[0]
        newloc=wave[jbg.closest(wave,event.xdata)]
        newvel=(newloc-nrw*(1.+newz))/(nrw*(1.+newz))*c
        print 'Velocity:',round(newvel,3)
        veltog=0
        finishaddnewline(newvel)

def keypress(event):
    global cont,fittog,EWtog,wave1,wave2,wrange,ewid,cid,filename,fitcoeff,fitcovmtx,zs,restwaves,normtog,VPtog,veltog,fitpars,fiterrors,parinfo
    global wave,flux,sigup,fitidx
    global normflux,normsig
    if event.key=='t':
       if fittog==0:
           print 'Selecting points for fit'
           fittog=1
       elif fittog==1:
           print 'Not selecting points for fit'
           fittog=0
    #if event.key=='f':
    #   print 'Fitting'
    #   fitcont(wave,flux,sigup,fitidx)
    #   updateplot()
    if event.key=='h':
        normtog=0
        updateplot()
    if event.key=='c':
        kernsize=3
        flux=np.convolve(flux,np.ones(kernsize)/kernsize)
        updateplot()
    if event.key=='b':
        #jbg.changefocus('terminal')
        wave1=float(raw_input('Enter beginning wavelength: \n'))
        wave2=float(raw_input('Enter ending wavelength: \n'))
        wrange=range(jbg.closest(wave,wave1),jbg.closest(wave,wave2))
        updateplot()
    #if event.key=='d':
       #dumpcontinuum(filename,wave1,wave2)
    #if event.key=='w':
    #   if EWtog==0:
    #       enterEWmode()
    #   elif EWtog==1:
    #       print 'Leaving equivalent width mode'
    #       EWtog=0
    #if event.key=='u':
    #   print 'Measure upper limit'
    #   changefocus('terminal')
    #   velbound=raw_input('Enter +/- velocity range. \n')
    #   calcEWupper(float(velbound))
    if event.key=='n':
       normspec()
       jbg.changefocus('terminal')
       answer=raw_input('Write out normalized spectrum? \n')
       if ((answer=='y')|(answer=='yes')): dumpnormspec()
    if event.key=='v':
        enterVPmode()
    if event.key=='p':
        normspec()
        VPtog=1
        print 'Fitting line profile(s)...'
        #fitpars,fiterrors=joebvpfit.joebvpfit(wave[wrange],normflux[wrange],normsig[wrange],fitpars,parinfo)
        fitpars,fiterrors=joebvpfit.joebvpfit(wave[cfg.fitidx],normflux[cfg.fitidx],normsig[cfg.fitidx],fitpars,parinfo)
        updateplot()
        #from PyQt4.QtCore import pyqtRemoveInputHook

        #pyqtRemoveInputHook()
        #pdb.set_trace()
    if event.key=='r':
        jbg.changefocus('terminal')
        answer=raw_input('Enter name of line parameter file. \n')
        try:
            parfiletable=np.genfromtxt('lineparfiles/'+answer,dtype=None,delimiter='|')
            if len(parfiletable[0])==5:
                restwaves,zs,cols,bs,vels=np.loadtxt('lineparfiles/'+answer,dtype=None,delimiter='|',unpack=True)
                sorted=np.argsort(restwaves)
                zs=zs[sorted]
                restwaves=restwaves[sorted]
                cols=cols[sorted]
                bs=bs[sorted]
                vels=vels[sorted]
                fitpars,parinfo=joebvpfit.initlinepars(zs,restwaves,initvals=[restwaves,cols,bs,zs,vels])
            else:
                restwaves,zs,cols,bs,vels,cflag,bflag,vflag=np.loadtxt('lineparfiles/'+answer,dtype=None,delimiter='|',unpack=True)
                cflag=cflag.astype(int) ; bflag=bflag.astype(int) ; vflag=vflag.astype(int)
                fitpars,parinfo=joebvpfit.initlinepars(zs,restwaves,initvals=[restwaves,cols,bs,zs,vels],initinfo=[cflag,bflag,vflag])
        except IOError: print 'Invalid filename!'
        jbg.changefocus('figure')
    updateplot()

def enterVPmode():
    global normflux,normsig,normtog,VPtog
    #normspec()
    obswaves=np.array(restwaves)*(1.+np.array(zs))
    wave1idx=jbg.closest(wave,wave1)
    wave2idx=jbg.closest(wave,wave2)
    prange=np.arange(wave1idx,wave2idx)
    updateplot(prange)
    VPtog=1
    VPoptions()

 
def VPoptions():
    global newrestwave,newz,veltog,newcol,newb,fitpars,parinfo,fiterrors
    jbg.changefocus('terminal')
    doubs=np.array(['C IV','Si IV','O VI','Si II','N V'])
    doublams=np.array([[1548.195,1550.77],[1393.755,1402.77],[1031.9261,1037.6167],[1190.4158,1193.2897],[1238.821,1242.804]])
    if len(fitpars)!=0:
            ans=raw_input('Add, change, delete, or write out lines and fit parameters? (a, c, d, or w) \n')
    else: ans='a'
    if ((ans=='a') | (ans=='A')):
        ans=raw_input('Approximate rest frame wavelength of the line to add or species to add doublet? \n')
        if ans in doubs:
           newrestwave=doublams[np.where(doubs==ans)[0]][0]
           print newrestwave
           lamidx1=jbg.closest(atomicdata.vernlam,newrestwave[0])
           lamidx2=jbg.closest(atomicdata.vernlam,newrestwave[1])
           ans=raw_input('Add the following lines:  '+atomicdata.vernion[lamidx1]+ ' ' +round(atomicdata.vernlam[lamidx1],2)+' '+round(atomicdata.vernlam[lamidx2],2) +' \n?')
        else:
           lamidx=jbg.closest(atomicdata.vernlam,float(ans))
           newrestwave=atomicdata.vernlam[lamidx]
           ans=raw_input('Add the following line:  '+atomicdata.vernion[lamidx]+ ' ' +round(atomicdata.vernlam[lamidx],2)+' \n?')
        if ((ans=='y') | (ans=='Y')):
          ans=raw_input('Redshift of new line? \n') ; newz=float(ans)
          ans=raw_input('Column density? \n') ; newcol=float(ans)
          ans=raw_input('b parameter? \n') ; newb=float(ans)
          ans=raw_input('Click on location of line.')
          veltog=1
          jbg.changefocus('figure')
          ### This operation continues in new function to intercept click

    elif ((ans=='c') | (ans=='C')):
       listlines()
       ans1=raw_input('Which line(s) to modify, or \'a\' to modify all? \n')
       if ans1=='a':
           ans2=raw_input('Vary all parameters, vary velocities, velocity flop, fix velocities, or restore parameters? (va, vv, vf, fv, r)')
           if ans2=='va':
               for i in range(len(parinfo[4])):
                  parinfo[1][i]=0; parinfo[2][i]=0; parinfo[4][i]=0
           elif ans2=='vv':
               for i in range(len(parinfo[4])):
                  parinfo[4][i]=0
           elif ans2=='vf':
               for i in range(len(parinfo[4])):
                  parinfo[4][i]=-5
           elif ans2=='fv':
               for i in range(len(parinfo[4])):
                  parinfo[4][i]=1
           elif ans2=='r':
               parinfo=cfg.origparinfo
       else:
           idx=ans1.split()
           idx=map(int,idx)
           ans=raw_input('Change the column density, b parameter, velocity or parameter associations? (c, b, v, or a) \n')
           if ((ans=='c') | (ans=='C')):
                  try:
                        newcol=float(raw_input('New column density? \n'))
                        for i in idx: fitpars[1][i]=float(newcol)
                  except ValueError: print 'Invalid value!'
           elif ((ans=='b') | (ans=='B')):
                        try:
                                newb=float(raw_input('New b parameter? \n'))
                                for i in idx: fitpars[2][i]=float(newb)
                        except ValueError: print 'Invalid value!'
           elif ((ans=='v') | (ans=='V')):
                        try:
                                newvel=float(raw_input('New component velocity? \n'))
                                for i in idx: fitpars[4][i]=float(newvel)
                        except ValueError: print 'Invalid value!'
           elif ((ans=='a') | (ans=='A')):
                        try:
                                newb=raw_input('Enter column density, b, and velocity flags separated by spaces. \n')
                                flags=newb.split()
                                parinfo[1][idx]=int(flags[0])
                                parinfo[2][idx]=int(flags[1])
                                parinfo[4][idx]=int(flags[2])
                        except ValueError: print 'Invalid value!'
           listlines()
    elif ((ans=='d') | (ans=='D')):
            listlines()
            dellines=raw_input('Which line(s) to delete? (separated by spaces) \n')
            dellist=dellines.split()
            delintlist=[]
            for idx in dellist:
                delintlist.append(int(idx))
            if len(fitpars[0])==1:
                print 'Hey now!'
                fitpars=[[],[],[],[],[]]; fiterrors=[[],[],[],[],[]]; parinfo=[[],[],[],[],[]]
            else:
                fitpars=np.delete(fitpars,delintlist,axis=1)
                parinfo=np.delete(parinfo,delintlist,axis=1)
            listlines()
    elif ((ans=='w') | (ans=='W')):
        writelinepars(filename,fitpars,fiterrors)
        updateplot()

def finishaddnewline(newvel):
    global newrestwave,newz,veltog,newcol,newb,parinfo,fitpars,fiterrors
    ### Show list of lines before concatenating to prevent errors
    #listlines()
    ### Initialize fitpars and fiterrors if necessary
    if len(fitpars)==0:
        fitpars=[[],[],[],[],[]] ; fiterrors=[[],[],[],[],[]] ; parinfo=[[],[],[],[],[]]

    ### Now concatenate and carry on
    if isinstance(newrestwave,float):
        fitpars[0].extend([newrestwave])
        fitpars[1].extend([newcol])
        fitpars[2].extend([newb])
        fitpars[3].extend([newz])
        fitpars[4].extend([newvel])
        fiterrors[0].extend([0.]) ;  fiterrors[1].extend([0.]) ;  fiterrors[2].extend([0.])
        fiterrors[3].extend([0.]) ;  fiterrors[4].extend([0.])
    else:
        for i in range(len(newrestwave)):
            fitpars[0].extend([newrestwave[i]])
            fitpars[1].extend([newcol])
            fitpars[2].extend([newb])
            fitpars[3].extend([newz])
            fitpars[4].extend([newvel])
            fiterrors[0].extend([0.]) ;  fiterrors[1].extend([0.]) ;  fiterrors[2].extend([0.])
            fiterrors[3].extend([0.]) ;  fiterrors[4].extend([0.])
    jbg.changefocus('terminal')
    ans=raw_input('Tie this line to any other? \n')
    if isinstance(ans,int):
        if int(ans) in range(len(fitpars[0])):
            idx=int(ans)
            if parinfo[1][idx] > 1:
                parinfo=np.concatenate([parinfo,[[1],[idx],[idx],[1],[idx]]],axis=1)
            else:
                nexttie=max[parinfo[1]]+1
                parinfo=np.concatenate([parinfo,[[1],[nexttie],[nexttie],[1],[nexttie]]],axis=1)
                parinfo[1][idx]=nexttie ; parinfo[2][idx]=nexttie ; parinfo[4][idx]=nexttie
        else:  parinfo=np.concatenate([parinfo,[[1],[0],[0],[1],[0]]],axis=1)
    elif isinstance(newrestwave,np.ndarray):
        for i in range(len(newrestwave)):
            parinfo=np.concatenate([parinfo,[[1],[0],[0],[1],[0]]],axis=1)
   
    else:
        parinfo=np.concatenate([parinfo,[[1],[0],[0],[1],[0]]],axis=1)
    listlines()
    jbg.changefocus('figure')

def writelinepars(filename,vfitpars,vfiterrors,writeoutfile='default'):
    global fitpars,fiterrors
    fitpars=vfitpars ; fiterrors=vfiterrors
    bigfiletowrite=cfg.largeVPparfile
    if writeoutfile=='default':
        filetowrite=cfg.VPparoutfile
    else:
        filetowrite=writeoutfile
    if os.path.isfile(filetowrite):
            VPparfile=open(filetowrite,'ab')
            bigparfile=open(bigfiletowrite,'ab')
    else:
            VPparfile=open(filetowrite,'wb')
            bigparfile=open(bigfiletowrite,'wb')
    for i in range(len(fitpars[0])):
        towrite=jbg.pipedelimrow([filename,fitpars[0][i],fitpars[3][i],fitpars[1][i],fiterrors[1][i],fitpars[2][i],fiterrors[2][i],fitpars[4][i],fiterrors[4][i]])
        VPparfile.write(towrite)
        bigparfile.write(towrite)
    VPparfile.close()
    bigparfile.close()
    print 'Line parameters written to:'
    print filetowrite

def listlines():
    print 'Current lines'
    for i in range(len(fitpars[0])):
        print i,fitpars[0][i],jbg.decimalplaces(fitpars[3][i],5),jbg.decimalplaces(fitpars[1][i],3),jbg.decimalplaces(fitpars[2][i],3),jbg.decimalplaces(fitpars[4][i],3),'\t',parinfo[1][i],parinfo[2][i],parinfo[4][i]

def normspec():
    global wrange,cont,normflux,normsig,normtog
    normtog=1
    #normflux=flux[wrange]/cont
    #normsig=sigup[wrange]/cont
    updateplot()

def dumpnormspec():
    global wrange,cont,normflux,normsig,normtog
    global wave1,wave2
    normflux=flux[wrange]/cont
    normsig=sigup[wrange]/cont
    newfilename='normspec-'+cfg.field+'-'+str(round(wave1,1))+'-'+str(round(wave2,1))+'.dat'
    np.savetxt(newfilename,np.transpose((wave[wrange],normflux,normsig)))
    print 'Normalized spectrum written to:',newfilename
    jbg.changefocus('figure')

def fitlegendre(wavepts,fluxpts,sigpts,order):
    vander=L.legvander(wavepts,order)
    design=np.zeros(vander.shape)
    for i in range(len(wavepts)):
        design[i]=vander[i]/sigpts[i]
    U,s,v=svd(design,full_matrices=False,compute_uv=True)
    V=v.transpose()
    solvec=np.zeros(order+1)
    for i in range(order+1):
        solvec+=np.dot(U[:,i],fluxpts/sigpts)/s[i]*V[:,i]
    ### Build covariance matrix
    covmtx=np.zeros([order+1,order+1])
    for j in range(order+1):
        for k in range(order+1):
            covmtx[j,k]=sum(V[:,j]*V[:,k]/s)
    return solvec,covmtx

def redchisq(fluxpts1,fitpts1,sigpts,order):
    Xsq=0.
    diff=sum(fluxpts1-fitpts1)
    df=len(fluxpts1)-order-1
    Xsq=1./df*sum(diff**2/sigpts**2)
    return Xsq,df

def updateplot(plotrange='initial',numchunks=8):
    global wave,flux,fitidx,wrange,fitcoeff,normtog,VPtog,ax,fitpars
    global fig,spls
    if plotrange=='initial': prange=wrange
    else: prange=plotrange
    #cont=L.legval(wave[prange],fitcoeff)
    wlen=(waveidx2-waveidx1)/numchunks
    if normtog==1:
        for i,sp in enumerate(spls):
            sp.clear()
            prange=range(waveidx1+i*wlen,waveidx1+(i+1)*wlen)

            if ((VPtog==1)&(len(fitpars[0])>0)):
                #sp.plot(wave[prange],normflux[prange],linestyle='steps')
                #model=joebvpfit.voigtfunc(wave[prange],fitpars)
                #res=normflux[prange]-model
                #sp.plot(wave[prange],model,'r')
                #sp.plot(wave[prange],res,'.',color='black')
                #sp.plot(wave[prange],[0]*len(prange),color='gray')

                sp.plot(wave,normflux,linestyle='steps')
                model=joebvpfit.voigtfunc(wave,fitpars)
                res=normflux-model
                sp.plot(wave,model,'r')
                sp.plot(wave,res,'.',color='black')
                sp.plot(wave,[0]*len(wave),color='gray')
                sp.set_ylim(-0.5,1.3)
                sp.set_xlim(wave[prange[0]],wave[prange[-1]])
                ### label lines we are trying to fit

                for j in range(len(fitpars[0])):
                    labelloc=fitpars[0][j]*(1.+fitpars[3][j])+fitpars[4][j]/c*fitpars[0][j]*(1.+fitpars[3][j])
                    label='- '+str(round(fitpars[0][j],1))+'  z = '+str(round(fitpars[3][j],5))
                    sp.text(labelloc,1.0,label,rotation='vertical',withdash=True,horizontalalignment='center',va='bottom',clip_on=True)

            #sp.plot(wave[prange],normsig[prange],linestyle='steps',color='red')
            #sp.plot(wave[prange],-normsig[prange],linestyle='steps',color='red')

            sp.plot(wave,normsig,linestyle='steps',color='red')
            sp.plot(wave,-normsig,linestyle='steps',color='red')


        #ax.clear()
        #ax.plot(wave[prange],normflux[prange],linestyle='steps')
        #ax.plot(wave[prange],normsig[prange],linestyle='steps',color='red')
        #if ((VPtog==1)&(len(fitpars[0])>0)):
        #    model=voigtfunc(wave[prange],fitpars)
        #    ax.plot(wave[prange],model,'r')
        #    ### label lines we are trying to fit
        #    logobswaves=restwaves*(1.+zs)
        #    for i in range(len(fitpars[0])):
        #        labelloc=fitpars[0][i]*(1.+fitpars[3][i])+fitpars[4][i]/c*fitpars[0][i]*(1.+fitpars[3][i])
        #        label='- '+str(round(fitpars[0][i],1))+'  z = '+str(round(fitpars[3][i],5))
        #        text(labelloc,1.0,label,rotation='vertical',withdash=True,horizontalalignment='center',va='bottom')

    else:
        ax.clear()
        ax.plot(wave[prange],normflux[prange],linestyle='steps')
        if plotrange=='initial': ax.plot(wave[fitidx],flux[fitidx],'gd')
        ax.plot(wave[prange],sigup[prange],linestyle='steps',color='red')
        plt.xlabel('wavelength')
    plt.ylabel('flux')
    plt.title(cfg.field)
    #ticklocs=[wave[jbg.closest(vel,-velrange)],wave[jbg.closest(vel,-velrange/2.)],wave[jbg.closest(vel,0)],wave[jbg.closest(vel,velrange/2.)],wave[jbg.closest(vel,velrange)]]
    #plt.xticks(ticklocs,[-velrange,-velrange/2.,0,velrange,velrange/2.])
    plt.draw()
    plt.tight_layout()
    #jbg.changefocus('figure')

def initplot(fig, wave1, wave2,numchunks=8):
    global wave,flux,normflux,spls,wrange,waveidx1,waveidx2
    #waveidx1=jbg.closest(wave,wave1)
    #waveidx2=jbg.closest(wave,wave2)
    #wrange=range(waveidx1,waveidx2)
    wlen=(waveidx2-waveidx1)/numchunks
    spls=[]
    sg=jbg.subplotgrid(numchunks)
    for i in range(numchunks):
        spls.append(fig.add_subplot(sg[i][0],sg[i][1],sg[i][2]))
        pixs=range(waveidx1+i*wlen,waveidx1+(i+1)*wlen)
        spls[i].plot(wave[pixs],normflux[pixs],linestyle='steps')
    plt.tight_layout()

###################################################################

fittog=0
EWtog=0
normtog=0
VPtog=0
veltog=0

###################################################################

#########################

# Begin main program

#########################

### This block was used when we had continuum fitting was part of the code
# def activefit(zs_arg,restwaves_arg,waveone,wavetwo):
#     global zs,restwaves,wave,flux,sigup,fitidx,wrange,ax,filename
#     global normflux,normsig
#     global fitpars,fiterrors,parinfo
#     global wave1,wave2
#     global fig
#
#     zs=zs_arg ; restwaves=restwaves_arg
#     wave=cfg.wave ; flux=cfg.flux ; sigup=cfg.sigup ; filename=cfg.filename
#     wave1=waveone ; wave2=wavetwo
#     linepars,parinfo=joebvpfit.initlinepars(zs,restwaves)
#     waveidx1=jbg.closest(wave,waveone)
#     waveidx2=jbg.closest(wave,wavetwo)
#     wrange=range(waveidx1,waveidx2)
#     fitcoeff,fitcovmtx=initcont(wave,flux,sigup,wave1,wave2)
#     contwave=cfg.wave[wrange]
#     contflux=cfg.flux[wrange]
#     conterr=cfg.sigup[wrange]
#     cont=L.legval(contwave,fitcoeff)
#     normflux=contflux/cont
#     normsig=conterr/cont
#     fitpars,fiterrors=joebvpfit.joebvpfit(contwave,normflux,normsig,linepars,parinfo)
#     fig=plt.figure()
#     cid=fig.canvas.mpl_connect('button_press_event', onclick)
#     kid=fig.canvas.mpl_connect('key_press_event', keypress)
#     ax=fig.add_subplot(111)
#     ax.plot(wave[wrange],flux[wrange],linestyle='steps')
#     ax.plot(wave[fitidx],flux[fitidx],'gd')
#     plt.xlabel('wavelength')
#     plt.ylabel('flux')
#     plt.title(cfg.field)
#     ax.plot(wave[wrange],cont)
#     plt.show()
#     #dumpcontinuum(filename,wave1,wave2)
#     writelinepars(filename,fitpars,fiterrors)


if __name__=='__main__':

    global fitidx,fig
    velrange=1800.
    print sys.argv
    ### Handle command line arguments
    if len(sys.argv) < 4:
        print "Fit a continuum for a spectrum region and Voigt profile fit absorption lines."
        print "Parameters: filename restwave z"
        print "Or: filename 'waverange' startwave endwave"
        print "Or: filename 'inputlines' startwave endwave lineinputfile"
        sys.exit()
    elif sys.argv[2]=='inputlines':
        print 'hello'
        filename=str(sys.argv[1])
        parfilename=str(sys.argv[5])
        if filename[-5:]=='.fits':
            spectrum=pf.open(filename)
            if len(spectrum)==4:
                wave=spectrum[2].data
                flux=spectrum[0].data
                sigup=spectrum[1].data
                cont=spectrum[3].data
                normflux=flux/cont
                normsig=sigup/cont
                normtog=1
            else:
                wave=spectrum[1].data[0][0]
                flux=spectrum[1].data[0][1]
                sigup=spectrum[1].data[0][3]

        else:
            sys.exit()
            #wave,flux,sigup=loadtxt(filename,dtype=None,unpack=True)
            #filename=filename[:-4]+'.fits'
        wave1=float(sys.argv[3])
        wave2=float(sys.argv[4])

        ### Read in lines to fit
        inputparfile=ascii.read(parfilename)

        linefile=inputparfile['specfile'].data
        linerestwave=inputparfile['restwave'].data
        linez=inputparfile['zsys'].data
        linecol=inputparfile['col'].data
        lineb=inputparfile['bval'].data
        linevel=inputparfile['vel'].data
        linevlim1=inputparfile['vlim1'].data
        linevlim2=inputparfile['vlim2'].data
        colflag=inputparfile['nflag'].data
        bflag=inputparfile['bflag'].data
        velflag=inputparfile['vflag'].data
        pix1=inputparfile['pix1'].data
        pix2=inputparfile['pix2'].data

        ### Find lines in wavelength range and initialize
        lineobswave=linerestwave*(1.+linez)
        restwave=0.
        z=0.
        lineshere=np.where((lineobswave>wave1)&(lineobswave<wave2))[0]
        zs=linez[lineshere]
        restwaves=linerestwave[lineshere]
        linecol=linecol[lineshere]
        lineb=lineb[lineshere]
        linevel=linevel[lineshere]
        linevlim1=linevlim1[lineshere] ; linevlim2=linevlim2[lineshere]
        colflag=colflag[lineshere];bflag=bflag[lineshere];velflag=velflag[lineshere]
        allpars=np.core.records.fromarrays([restwaves,linecol,lineb,zs,linevel,atomicdata.lam2ion(restwaves),linevlim1,linevlim2,colflag,bflag,velflag],names='lamrest,col,b,z,vel,ion,vlim1,vlim2,colflag,bflag,velflag',formats='f8,f8,f8,f8,f8,a4,f8,f8,i4,i4,i4')
        allpars.sort(order=['ion','z','vel','lamrest'])
        linerestwave=allpars['lamrest']
        zs=allpars['z']
        linecol=allpars['col']
        lineb=allpars['b']
        linevel=allpars['vel']
        linevlim1=allpars['vlim1'] ; linevlim2=allpars['vlim2']
        colflag=allpars['colflag'] ; bflag=allpars['bflag'] ; velflag=allpars['velflag']
        #pix1=allpars['pix1'] ; pix2=allpars['pix2']
        restwaves=linerestwave
        initinfo=[colflag,bflag,velflag]
        initpars=[restwaves,linecol,lineb,zs,linevel,linevlim1,linevlim2]
        linepars,parinfo=joebvpfit.initlinepars(zs,restwaves,initpars,initinfo=initinfo)
        origpars=linepars ; origparinfo=parinfo
        vel=[]
        linelistidx=-99.
        cfg.fitidx=joebvpfit.fitpix(wave,linepars)

        waveidx1=jbg.closest(wave,wave1)
        waveidx2=jbg.closest(wave,wave2)
        wrange=range(waveidx1,waveidx2)

        fitpars=linepars

        fig=plt.figure(figsize=(13.5,12))
        cid=fig.canvas.mpl_connect('button_press_event', onclick)
        kid=fig.canvas.mpl_connect('key_press_event', keypress)

        #ax=fig.add_subplot(111)
        initplot(fig,wave1,wave2)

        #ax.plot(wave[wrange],flux[wrange],linestyle='steps')
        #ax.plot(wave[wrange],normflux[wrange],linestyle='steps')
        #xlabel('wavelength')
        #ylabel('flux')

        plt.show()
