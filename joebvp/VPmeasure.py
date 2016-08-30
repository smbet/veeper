# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 22:50:45 2016

@author: burchett
"""

from PyQt4.uic import loadUiType
from PyQt4 import QtGui
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import joebvp.atomicdata as atomicdata
import joebvp.joebgoodies as jbg
from joebvp import cfg
from joebvp import joebvpfit
import os
from linetools.spectra.io import readspec
from astropy.io import ascii
import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

modpath=os.path.abspath(os.path.dirname(__file__))
print os.path.abspath(os.path.dirname(__file__))
Ui_MainWindow, QMainWindow = loadUiType(modpath+'/mainvpwindow.ui')

c=cfg.c/1e5

class LineParTableModel(QAbstractTableModel):
    def __init__(self,fitpars,fiterrors,parinfo,parent=None):
            QAbstractTableModel.__init__(self,parent)
            self.fitpars=fitpars
            self.fiterrors=fiterrors
            self.parinfo=parinfo
            self.headers=['Wavelength','Species','z','N','sig(N)','b','sig(b)','v','sig(v)']

    def rowCount(self,parent):
        return len(self.fitpars[0])

    def columnCount(self, parent):
        return 9
        #return len(self.fitpars)

    def writelinepars(self, outfilename,specfilename):
        fitpars = self.fitpars
        fiterrors = self.fiterrors
        parinfo = self.parinfo
        bigfiletowrite = cfg.largeVPparfile
        filetowrite = outfilename

        if os.path.isfile(filetowrite):
            VPparfile = open(filetowrite, 'wb')
            bigparfile = open(bigfiletowrite, 'ab')
        else:
            VPparfile = open(filetowrite, 'wb')
            bigparfile = open(bigfiletowrite, 'wb')
        header='specfile|restwave|zsys|col|sigcol|bval|sigbval|vel|sigvel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|pix1|pix2|trans \n'
        VPparfile.write(header)
        bigparfile.write(header)
        for i in range(len(fitpars[0])):
            zline=fitpars[3][i]
            vlim1=fitpars[5][i] ; vlim2=fitpars[6][i]
            restwave=fitpars[0][i]
            wobs1=restwave*(1+zline+vlim1/299792.458)
            wobs2=restwave*(1+zline+vlim2/299792.458)
            pix1=jbg.closest(cfg.wave,wobs1)
            pix2=jbg.closest(cfg.wave,wobs2)
            trans=atomicdata.lam2ion(fitpars[0][i])
            towrite=jbg.pipedelimrow([cfg.filename,restwave,round(zline,5),round(fitpars[1][i],3),round(fiterrors[1][i],3),round(fitpars[2][i],3),round(fiterrors[2][i],3),round(fitpars[4][i],3), round(fiterrors[4][i],3),parinfo[1][i],parinfo[2][i],parinfo[4][i],vlim1,vlim2,wobs1,wobs2,pix1,pix2,trans])
            VPparfile.write(towrite)
            bigparfile.write(towrite)
        VPparfile.close()
        bigparfile.close()
        print 'Line parameters written to:'
        print filetowrite

    def writeVPmodel(self, outfile, wave,fitpars,normflux,normsig):
        from astropy.table import Table
        model = joebvpfit.voigtfunc(wave, fitpars)
        modeltab=Table([wave,model,normflux,normsig],names=['wavelength','model','normflux','normsig'])
        modeltab.write(outfile,format='fits')
        print 'Voigt profile model written to:'
        print outfile


    def data(self,index,role):

        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            if column == 1: toret=atomicdata.lam2ion(self.fitpars[0][row])
            elif column == 2: toret=round(self.fitpars[3][row],5)
            elif column == 3: toret=round(self.fitpars[1][row],3)
            elif column == 4: toret=round(self.fiterrors[1][row],3)
            elif column == 5: toret=round(self.fitpars[2][row],3)
            elif column == 6: toret=round(self.fiterrors[2][row],3)
            elif column == 7: toret=round(self.fitpars[4][row],3)
            elif column == 8: toret=round(self.fiterrors[4][row],3)
            else: toret=QVariant(round(self.fitpars[column][row],3))
            return toret

        if index.isValid() and role == Qt.DisplayRole:
            row = index.row()
            column = index.column()
            if column == 1: toret=atomicdata.lam2ion(self.fitpars[0][row])
            elif column == 2: toret=round(self.fitpars[3][row],5)
            elif column == 3: toret=round(self.fitpars[1][row],3)
            elif column == 4: toret=round(self.fiterrors[1][row],3)
            elif column == 5: toret=round(self.fitpars[2][row],3)
            elif column == 6: toret=round(self.fiterrors[2][row],3)
            elif column == 7: toret=round(self.fitpars[4][row],3)
            elif column == 8: toret=round(self.fiterrors[4][row],3)
            else: toret=QVariant(round(self.fitpars[column][row],3))
            return toret
        else: return None

    def flags(self,index):
        return Qt.ItemIsEditable | Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def setData(self,index,value,role=Qt.EditRole):
        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            if value.toFloat()[1]==True:
                val = round(float(value.toFloat()[0]),2)
                if column == 1: pass
                if column == 2: self.fitpars[3][row] = val
                elif column == 3: self.fitpars[1][row] = val
                elif column == 4: self.fiterrors[1][row] = val
                elif column == 5: self.fitpars[2][row] = val
                elif column == 6: self.fiterrors[2][row] = val
                elif column == 7: self.fitpars[4][row] = val
                elif column == 8: self.fiterrors[4][row] = val
                else: self.fitpars[column][row] = val
                self.dataChanged.emit(index,index)
                return True

        return False


    def addLine(self,wave,newrestwave,newz,newcol,newb,newvel,newvel1,newvel2):
        ### Setup new parameters and append them to main arrays of the data model
        newpars=[[newrestwave],[newcol],[newb],[newz],[newvel],[newvel1],[newvel2]]
        fitpars=np.hstack((self.fitpars,newpars))
        newerrors=[[-99],[-99],[-99],[-99],[-99]]
        fiterrors=np.hstack((self.fiterrors,newerrors))
        newindex=np.max(self.parinfo[1])+1
        newinfo=[[1],[newindex],[newindex],[1],[newindex]]
        parinfo=np.hstack((self.parinfo,newinfo))
        ### Call initlinepars to set atomic data in cfg.fosc, etc.
        junk,junk=joebvpfit.initlinepars(fitpars[3],fitpars[0],initvals=fitpars,initinfo=parinfo)
        ### Do the update
        midx=QModelIndex()  # Needed for 'beginInsertRows'
        self.beginInsertRows(midx,len(self.fitpars[0]),len(self.fitpars[0]))
        self.updatedata(fitpars,fiterrors,parinfo)
        self.endInsertRows()
        ### Reset pixels for fit and wavegroups for convolution
        cfg.fitidx=joebvpfit.fitpix(wave,fitpars) #Reset pixels for fit
        cfg.wavegroups=[]
    '''
    def insertRows(self,position,rows,parent=QModelIndex()):
        self.beginInsertRows(parent,position,position+rows-1)
        if isinstance(newrestwave,float):
            self.fitpars[0].extend([newrestwave])
            self.fitpars[1].extend([newcol])
            self.fitpars[2].extend([newb])
            self.fitpars[3].extend([newz])
            self.fitpars[4].extend([newvel])
            self.fitpars[5].extend([-cfg.defaultvlim])
            self.fitpars[6].extend([cfg.defaultvlim])
            self.fiterrors[0].extend([0.]) ;  self.fiterrors[1].extend([0.]) ;  self.fiterrors[2].extend([0.])
            self.fiterrors[3].extend([0.]) ;  self.fiterrors[4].extend([0.])
        for i in range(rows):
            defaultFitpars=[[0]*len(self.fitpars)]
            defaultFiterrors = [[0] * len(self.fiterrors)]
            defaultParinfo = [[0] * len(self.parinfo)]
            self.fitpars.insert(position,defaultFitpars)
        self.endInsertRows()
    '''

    def updatedata(self,fitpars,fiterrors,parinfo):
        self.fitpars=fitpars
        self.fiterrors=fiterrors
        self.parinfo=parinfo
        self.dataChanged.emit(self.index(0,0),self.index(self.rowCount(self)-1,self.columnCount(self)-1))


    def headerData(self,section,orientation,role):
        if role == Qt.DisplayRole:
            if orientation==Qt.Horizontal:
                return self.headers[section]

class newLineDialog(QDialog):
    def __init__(self, parent = None):
        super(newLineDialog,self).__init__(parent)

        layout = QGridLayout(self)

        self.lamLabel=QLabel(self)
        self.lamLabel.setText('Rest Wavelength:')
        layout.addWidget(self.lamLabel,0,0)
        self.lamBox = QLineEdit(self)
        #self.lamBox.setMaximumWidth = 30
        self.lamBox.editingFinished.connect(self.validateWavelength)
        layout.addWidget(self.lamBox,1,0)
        self.zLabel=QLabel(self)
        self.zLabel.setText('z:')
        layout.addWidget(self.zLabel,0,1)
        self.zBox = QLineEdit(self)
        #self.lamBox.setMaximumWidth = 30
        layout.addWidget(self.zBox,1,1)
        self.ionLabel = QLabel(self)
        layout.addWidget(self.ionLabel, 1, 2)
        
        self.colLabel=QLabel(self)
        self.colLabel.setText('N:')
        layout.addWidget(self.colLabel,2,0)
        self.colBox = QLineEdit(self)
        self.colBox.setMaximumWidth = 30
        layout.addWidget(self.colBox,3,0)
        self.bLabel=QLabel(self)
        self.bLabel.setText('b:')
        layout.addWidget(self.bLabel,2,1)
        self.bBox = QLineEdit(self)
        layout.addWidget(self.bBox,3,1)
        self.velLabel=QLabel(self)
        self.velLabel.setText('vel:')
        layout.addWidget(self.velLabel,2,2)
        self.velBox = QLineEdit(self)
        layout.addWidget(self.velBox,3,2)
        
        

        buttons = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel,
            Qt.Horizontal, self)
        layout.addWidget(buttons)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        #import pdb
        #pdb.set_trace()

    def validateWavelength(self):
        ### Handle approximate wavelength entries by looking up precise values on focus exit
        if self.lamBox.text()!='':
            restlam=self.lamBox.text()
            restlam=float(restlam)
            try:
                lam,fosc,gam=joebvpfit.setatomicdata([restlam],precise=False)
                self.lamBox.setText(str(lam[0]))
                self.ionLabel.setText(atomicdata.lam2ion(lam[0]))
            except:
                pass

    def lineParams(self):
        vel=float(self.velBox.text())
        vel1=vel-cfg.defaultvlim ; vel2=vel+cfg.defaultvlim
        return self.lamBox.text(),self.zBox.text(),self.colBox.text(),self.bBox.text(),self.velBox.text(),vel1,vel2

    @staticmethod
    def get_newline(parent = None):
        dialog = newLineDialog(parent)
        result = dialog.exec_()
        if result == 1:
            newlam,newz,newcol,newb,newvel,newvel1,newvel2 = dialog.lineParams()
            return newlam,newz,newcol,newb,newvel,newvel1,newvel2
        else:
            return 0

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self,specfilename,parfilename=None,wave1=None,wave2=None,numchunks=8,parent=None):
        QtGui.QMainWindow.__init__(self, parent)
        #super(Main,self).__init__()
        self.setupUi(self)

        ### Initialize stuff
        self.line_dict = {}
        self.fitpars = None
        self.parinfo = None
        self.wave1 = wave1
        self.wave2 = wave2
        self.numchunks = numchunks
        self.spls = []
        self.labeltog = 1
        self.pixtog = 0
        self.restog = 1
        self.lastclick=1334.

        ### Read in spectrum and list of lines to fit
        self.specfilename=specfilename
        self.spectrum = readspec(specfilename)
        self.wave=self.spectrum.wavelength.value
        self.normflux=self.spectrum.flux/self.spectrum.co
        self.normsig=self.spectrum.sig/self.spectrum.co
        cfg.wave=self.wave
        cfg.normflux=self.normflux
        cfg.filename=self.specfilename

        if not parfilename==None:
            self.linelist = self.initialpars(parfilename)


        ### Connect signals to slots
        self.fitButton.clicked.connect(self.fitlines)
        self.boxLineLabel.clicked.connect(self.toglabels)
        self.boxFitpix.clicked.connect(self.togfitpix)
        self.boxResiduals.clicked.connect(self.togresiduals)
        self.loadParsButton.clicked.connect(self.openParFileDialog)
        self.addLineButton.clicked.connect(self.addLineDialog)
        self.writeParsButton.clicked.connect(self.writeParFileDialog)
        self.writeModelButton.clicked.connect(self.writeModelFileDialog)

        ### Initialize spectral plots
        fig=Figure()
        self.fig=fig
        self.initplot(fig)

        ### Initializea side plot
        sidefig=Figure(figsize=(5.25,2))
        self.sidefig = sidefig
        self.addsidempl(self.sidefig)
        self.sideplot(self.lastclick)  #Dummy initial cenwave setting



    '''
    def poptable(self):
        ### Populate the line parameter table
        for i in range(len(self.fitpars[0])):
            self.lineList.setItem(i,0,QTableWidgetItem(str(self.fitpars[0][i])))
            self.lineList.setItem(i,1,QTableWidgetItem(jbg.decimalplaces(self.fitpars[3][i],3)))
            self.lineList.setItem(i,2,QTableWidgetItem(atomicdata.lam2ion(self.fitpars[0][i])))
            self.lineList.setItem(i,3,QTableWidgetItem(jbg.decimalplaces(self.fitpars[1][i],3)))
            self.lineList.setItem(i,4,QTableWidgetItem(jbg.decimalplaces(self.fiterrors[1][i],3)))
            self.lineList.setItem(i,5,QTableWidgetItem(jbg.decimalplaces(self.fitpars[2][i],3)))
            self.lineList.setItem(i,6,QTableWidgetItem(jbg.decimalplaces(self.fiterrors[2][i],3)))
            self.lineList.setItem(i,7,QTableWidgetItem(jbg.decimalplaces(self.fitpars[4][i],3)))
            self.lineList.setItem(i,8,QTableWidgetItem(jbg.decimalplaces(self.fiterrors[4][i],3)))
    '''

    def initplot(self,fig,numchunks=8):
        wlen=len(self.spectrum.wavelength)/numchunks
        self.spls=[]
        if self.wave1==None:  waveidx1=0  # Default to plotting entire spectrum for now
        else: waveidx1=jbg.closest(wave,wave1)
        if self.fitpars!=None:
                model=joebvpfit.voigtfunc(self.wave,self.datamodel.fitpars)
        sg=jbg.subplotgrid(numchunks)
        for i in range(numchunks):
            self.spls.append(fig.add_subplot(sg[i][0],sg[i][1],sg[i][2]))
            pixs=range(waveidx1+i*wlen,waveidx1+(i+1)*wlen)
            self.spls[i].plot(self.wave[pixs],self.normflux[pixs],linestyle='steps-mid')
            if self.fitpars!=None:
                self.spls[i].plot(self.wave,model,'r')
            self.spls[i].set_xlim(self.wave[pixs[0]],self.wave[pixs[-1]])
            self.spls[i].set_ylim(cfg.ylim)
            self.spls[i].set_xlabel('wavelength',labelpad=0)
            self.spls[i].set_ylabel('relative flux',labelpad=-5)
        fig.subplots_adjust(top=0.98,bottom=0.05,left=0.08,right=0.97,wspace=0.15,hspace=0.25)
        self.addmpl(fig)

    def initialpars(self,parfilename):
        ### Deal with initial parameters from line input file
        self.linelist = ascii.read(parfilename)
        linerestwave=self.linelist['restwave'].data
        linez=self.linelist['zsys'].data
        lineobswave=linerestwave*(1.+linez)
        if self.wave1==None:
            lineshere=np.where((lineobswave>self.wave[0])&(lineobswave<self.wave[-1]))[0]
        else:
            lineshere=np.where((lineobswave>wave1)&(lineobswave<wave2))[0]
        zs=self.linelist['zsys'][lineshere]
        restwaves=linerestwave[lineshere]
        linecol=self.linelist['col'][lineshere]
        lineb=self.linelist['bval'][lineshere]
        linevel=self.linelist['vel'][lineshere]
        linevlim1=self.linelist['vlim1'][lineshere] ; linevlim2=self.linelist['vlim2'][lineshere]
        colflag=self.linelist['nflag'][lineshere];bflag=self.linelist['bflag'][lineshere];velflag=self.linelist['vflag'][lineshere]
        allpars=np.core.records.fromarrays([restwaves,linecol,lineb,zs,linevel,atomicdata.lam2ion(restwaves),linevlim1,linevlim2,colflag,bflag,velflag],names='lamrest,col,b,z,vel,ion,vlim1,vlim2,colflag,bflag,velflag',formats='f8,f8,f8,f8,f8,a4,f8,f8,i4,i4,i4')
        allpars.sort(order=['ion','z','vel','lamrest'])
        linerestwave=allpars['lamrest']
        zs=allpars['z']
        linecol=allpars['col']
        lineb=allpars['b']
        linevel=allpars['vel']
        linevlim1=allpars['vlim1'] ; linevlim2=allpars['vlim2']
        colflag=allpars['colflag'] ; bflag=allpars['bflag'] ; velflag=allpars['velflag']
        restwaves=linerestwave
        initinfo=[colflag,bflag,velflag]
        initpars=[restwaves,linecol,lineb,zs,linevel,linevlim1,linevlim2]
        self.fitpars,self.parinfo=joebvpfit.initlinepars(zs,restwaves,initpars,initinfo=initinfo)
        self.fiterrors=np.ones([5,len(self.fitpars[0])])*-99 #Initialize errors to junk
        cfg.fitidx=joebvpfit.fitpix(self.wave,self.fitpars) #Set pixels for fit
        cfg.wavegroups=[]
        self.datamodel = LineParTableModel(self.fitpars,self.fiterrors,self.parinfo)
        self.tableView.setModel(self.datamodel)

    def sideplot(self,cenwave,wavebuf=3):
        if len(self.sidefig.axes)==0:
            self.sideax=self.sidefig.add_subplot(111)
        self.sideax.clear()
        self.sideax.plot(self.wave, self.normflux, linestyle='steps-mid')
        if self.pixtog == 1:
            self.sideax.plot(self.wave[cfg.fitidx], self.normflux[cfg.fitidx], 'gs', markersize=4, mec='green')
        model = joebvpfit.voigtfunc(self.wave, self.fitpars)
        res = self.normflux - model
        self.sideax.plot(self.wave, model, 'r')
        if self.restog == 1:
            self.sideax.plot(self.wave, -res, '.', color='black', ms=2)
        self.sideax.plot(self.wave, [0] * len(self.wave), color='gray')

        ### label lines we are trying to fit
        if self.labeltog == 1:
            for j in range(len(self.fitpars[0])):
                labelloc = self.fitpars[0][j] * (1. + self.fitpars[3][j]) + self.fitpars[4][j] / c * \
                                                                            self.fitpars[0][j] * (
                                                                            1. + self.fitpars[3][j])
                label = ' {:.1f}_\nz{:.4f}'.format(self.fitpars[0][j], self.fitpars[3][j])
                self.sideax.text(labelloc, cfg.label_ypos, label, rotation=90, withdash=True, ha='center', va='bottom',
                        clip_on=True, fontsize=cfg.label_fontsize)

        self.sideax.plot(self.wave, self.normsig, linestyle='steps-mid', color='red', lw=0.5)
        self.sideax.plot(self.wave, -self.normsig, linestyle='steps-mid', color='red', lw=0.5)
        try:
            self.sideax.set_xlim(cenwave-wavebuf,cenwave+wavebuf)
            self.sideax.set_ylim(cfg.ylim)
            self.changesidefig(self.sidefig)
        except TypeError:
            pass
    def fitlines(self):
        print 'Fitting line profile(s)...'
        #self.fitpars,self.fiterrors=joebvpfit.joebvpfit(self.wave[cfg.fitidx],self.normflux[cfg.fitidx],self.normsig[cfg.fitidx],self.datamodel.fitpars,self.datamodel.parinfo)
        print len(self.fitpars[0]),len(self.fiterrors[0])
        self.fitpars, self.fiterrors = joebvpfit.joebvpfit(self.wave, self.normflux,self.normsig, self.datamodel.fitpars,self.datamodel.parinfo)
        self.datamodel.updatedata(self.fitpars,self.fiterrors,self.parinfo)
        self.tableView.resizeColumnsToContents()
        self.updateplot()
        self.sideplot(self.lastclick)


    def toglabels(self):
        if self.labeltog==1: self.labeltog=0
        else: self.labeltog=1
        self.updateplot()

    def togfitpix(self):
        if self.pixtog == 1:
            self.pixtog = 0
        else:
            self.pixtog = 1
        self.updateplot()

    def togresiduals(self):
        if self.restog == 1:
            self.restog = 0
        else:
            self.restog = 1
        self.updateplot()

    def openParFileDialog(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open line parameter file','.')
        fname = str(fname)
        if fname != '':
            self.initialpars(fname)

        self.updateplot()

    def writeParFileDialog(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save line parameter file', cfg.VPparoutfile)
        fname = str(fname)
        if fname != '':
            self.datamodel.writelinepars(fname,self.specfilename)

    def writeModelFileDialog(self):
        fname = QtGui.QFileDialog.getSaveFileName(self, 'Save model to file', cfg.VPmodeloutfile)
        fname = str(fname)
        if fname != '':
            self.datamodel.writeVPmodel(fname, self.wave, self.fitpars, self.normflux, self.normsig)

    def addLineDialog(self):
        dlgOutput=newLineDialog.get_newline()
        if (dlgOutput != 0):
            if '' not in dlgOutput:
                newlam,newz,newcol,newb,newvel,newvel1,newvel2 = dlgOutput
                self.datamodel.addLine(self.wave,float(newlam), float(newz), float(newcol), float(newb), float(newvel), float(newvel1), float(newvel2))        #dialog=newLineDialog(parent=None)
                self.fitpars = self.datamodel.fitpars
                self.fiterrors = self.datamodel.fiterrors
                self.parinfo = self.datamodel.parinfo
                self.tableView.setModel(self.datamodel)


    def updateplot(self):
        if self.wave1==None:  waveidx1=0  # Default to plotting entire spectrum for now
        else: waveidx1=jbg.closest(self.wave,self.wave1)
        wlen=len(self.spectrum.wavelength)/self.numchunks
        for i,sp in enumerate(self.spls):
                sp.clear()
                prange=range(waveidx1+i*wlen,waveidx1+(i+1)*wlen)
    
                if ((len(self.fitpars[0])>0)):
    
                    sp.plot(self.wave,self.normflux,linestyle='steps-mid')
                    if self.pixtog==1:
                        sp.plot(self.wave[cfg.fitidx],self.normflux[cfg.fitidx],'gs',markersize=4,mec='green')
                    model=joebvpfit.voigtfunc(self.wave,self.fitpars)
                    res=self.normflux-model
                    sp.plot(self.wave,model,'r')
                    if self.restog==1:
                        sp.plot(self.wave,-res,'.',color='black', ms=2)
                    sp.plot(self.wave,[0]*len(self.wave),color='gray')
    
                    ### label lines we are trying to fit
                    if self.labeltog==1:
                        for j in range(len(self.fitpars[0])):
                            labelloc=self.fitpars[0][j]*(1.+self.fitpars[3][j])+self.fitpars[4][j]/c*self.fitpars[0][j]*(1.+self.fitpars[3][j])
                            label = ' {:.1f}_\nz{:.4f}'.format(self.fitpars[0][j], self.fitpars[3][j])
                            sp.text(labelloc, cfg.label_ypos, label, rotation=90, withdash=True, ha='center', va='bottom', clip_on=True, fontsize=cfg.label_fontsize)
                
    
                sp.plot(self.wave,self.normsig,linestyle='steps-mid',color='red', lw=0.5)
                sp.plot(self.wave,-self.normsig,linestyle='steps-mid',color='red', lw=0.5)
                sp.set_ylim(cfg.ylim)
                sp.set_xlim(self.wave[prange[0]],self.wave[prange[-1]])
                sp.set_xlabel('wavelength (A)', fontsize=cfg.xy_fontsize, labelpad=cfg.x_labelpad)
                sp.set_ylabel('normalized flux', fontsize=cfg.xy_fontsize, labelpad=cfg.y_labelpad)
        self.changefig(self.fig)

    def changefig(self, item):
        #text = str(item.text())
        self.rmmpl()
        self.addmpl(self.fig)

    def changesidefig(self, item):
        #text = str(item.text())
        self.rmsidempl()
        self.addsidempl(self.sidefig)
        
        
    def on_click(self, event):
        self.lastclick=event.xdata
        self.sideplot(self.lastclick)
        
    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.canvas.mpl_connect('button_press_event',self.on_click)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,self.mplwindow,coordinates=True)
        self.mplvl.addWidget(self.toolbar)

    def addsidempl(self, sidefig):
        self.sidecanvas = FigureCanvas(sidefig)
        self.sidecanvas.setParent(self.sideMplWindow)
        if len(self.sidefig.axes) == 0:
            self.sidemplvl = QVBoxLayout()
        if len(self.sidefig.axes) != 0:
            self.sidemplvl.addWidget(self.sidecanvas)
        if len(self.sidefig.axes) == 0:
            self.sideMplWindow.setLayout(self.sidemplvl)
        self.sidecanvas.draw()

        
    def rmmpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()

    def rmsidempl(self, ):
        self.sidemplvl.removeWidget(self.sidecanvas)
        self.sidecanvas.close()
        #self.sideMplWindow.setLayout(None)
        #self.mplvl.removeWidget(self.toolbar)
        #self.toolbar.close()

def go(specfilename,parfilename):
    import sys
    import numpy as np
    from astropy.io import fits as pf
    from astropy.io import ascii
    from linetools.spectra.io import readspec

    app = QtGui.QApplication(sys.argv)
    main = Main(specfilename, parfilename)
    main.show()
    app.exec_()
    #sys.exit(app.exec_())
    #app.quit()

if __name__ == '__main__':
        import sys
        specfilename=str(sys.argv[1])
        parfilename=str(sys.argv[5])
        go(specfilename,parfilename)
