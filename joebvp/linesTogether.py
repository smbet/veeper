import numpy as np
from astropy.io import ascii
from goodies import closest
import pdb

buf = 20  #number of pixels between pix2 and pix1 for lines to be unblended

def picksyslines(z,comps):
    zval=comps['z'][closest(comps['z'],z)]
    zmatches=np.where(np.abs(comps['z']-zval)<0.0001)[0]
    return [comps['idxs'][idx] for idx in zmatches]

def picktranslines(trans,comps):
    transmatches=np.where(comps['trans']==trans)[0]
    return [comps['idxs'][idx] for idx in transmatches]


def getblended(lineidx,blends):
    for i in range(len(blends)):
        blmatch=None
        try:
            idx=blends[i].index(lineidx)
            blmatch=blends[i]
            break
        except:
            pass
    return blmatch

def findblends(linepars):
    blends = []
    for i in range(len(linepars)-1):
        if i == 0: 
            thisblend = [i]
        if (linepars['pix1'][i+1] - linepars['pix2'][i]) < buf: 
            thisblend.append(i+1)
        else:
            blends.append(thisblend)
            thisblend = [i+1]
    return blends

def compilecomps(linepars):
    uqzs=np.unique(linepars['zsys'])
    comps=[]
    compzs=[]
    comptrans=[]
    for idx,zz in enumerate(uqzs):
        thisz=np.where(linepars['zsys']==zz)[0]
        uqtrans=np.unique(linepars['trans'][thisz])
        for j,tt in enumerate(uqtrans):
                sametrans=np.where((linepars['zsys']==zz)&(linepars['trans']==tt))[0]
                comps.append(sametrans.tolist())
                compzs.append(zz)
                comptrans.append(tt)
    compdict = {'z': np.array(compzs), 'trans':np.array(comptrans), 'idxs':comps}
    return compdict
            
def compilelist_z(z,comps,blends):
    idxs=picksyslines(z,comps)
    idxs=np.concatenate(idxs)
    compilation=[]
    for i,idx in enumerate(idxs):
        #print idx
        bl=getblended(idx,blends)
        if bl!=None: compilation.append(np.array(bl))
    conc=np.concatenate(compilation)
    return np.unique(conc).tolist()

def compilelist_trans(trans,comps,blends):
    idxs=picktranslines(trans,comps)
    idxs=np.concatenate(idxs)
    compilation=[]
    for i,idx in enumerate(idxs):
        #print idx
        bl=getblended(idx,blends)
        if bl!=None: compilation.append(np.array(bl))
    conc=np.concatenate(compilation)
    return np.unique(conc).tolist()

def writeparfile(linepars,idx,filename):
    #outlinefile=open(filename,'wb')
    #outlinefile.write('specfile|restwave|zsys|col|bval|vel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|pix1|pix2|trans \n')
    ascii.write(linepars[idx],output=filename)

    
    

linepars=ascii.read('pg1407testsyslines.dat')
linepars.sort('pix1')

blends=findblends(linepars)
comps=compilecomps(linepars)





pdb.set_trace()

