# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 22:50:45 2016

@author: burchett
"""

from PyQt4.uic import loadUiType
from PyQt4.QtGui import *
from PyQt4.QtCore import *
import joebvp.atomicdata as atomicdata
import joebvp.joebgoodies as jbg
from joebvp import cfg
from joebvp import joebvpfit


from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
    
Ui_MainWindow, QMainWindow = loadUiType('mainvpwindow.ui')

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



    def data(self,index,role):

        if role == Qt.EditRole:
            row = index.row()
            column = index.column()
            if column == 1: toret=atomicdata.lam2ion(self.fitpars[0][row])
            elif column == 2: toret=round(self.fitpars[3][row],5)
            elif column == 3: toret=round(self.fitpars[1][row],5)
            elif column == 4: toret=round(self.fiterrors[1][row],5)
            elif column == 5: toret=round(self.fitpars[2][row],5)
            elif column == 6: toret=round(self.fiterrors[2][row],5)
            elif column == 7: toret=round(self.fitpars[4][row],5)
            elif column == 8: toret=round(self.fiterrors[4][row],5)
            else: toret=QVariant(round(self.fitpars[column][row],3))
            return toret

        if index.isValid() and role == Qt.DisplayRole:
            row = index.row()
            column = index.column()
            if column == 1: toret=atomicdata.lam2ion(self.fitpars[0][row])
            elif column == 2: toret=round(self.fitpars[3][row],5)
            elif column == 3: toret=round(self.fitpars[1][row],5)
            elif column == 4: toret=round(self.fiterrors[1][row],5)
            elif column == 5: toret=round(self.fitpars[2][row],5)
            elif column == 6: toret=round(self.fiterrors[2][row],5)
            elif column == 7: toret=round(self.fitpars[4][row],5)
            elif column == 8: toret=round(self.fiterrors[4][row],5)
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


    def updatedata(self,fitpars,fiterrors,parinfo):
        self.fitpars=fitpars
        self.fiterrors=fiterrors
        self.parinfo=parinfo
        self.dataChanged.emit(self.index(0,0),self.index(self.rowCount(self)-1,self.columnCount(self)-1))

    def headerData(self,section,orientation,role):
        if role == Qt.DisplayRole:
            if orientation==Qt.Horizontal:
                return self.headers[section]

class Main(QMainWindow, Ui_MainWindow):
    def __init__(self,spec,parfilename=None,wave1=None,wave2=None,numchunks=8):
        super(Main,self).__init__()
        self.setupUi(self)
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

        self.spectrum=spec
        self.wave=spec.wavelength.value
        self.normflux=spec.flux/spec.co
        self.normsig=spec.sig/spec.co

        if not parfilename==None:
            self.linelist = self.initialpars(parfilename)

        self.fitButton.clicked.connect(self.fitlines)
        self.boxLineLabel.clicked.connect(self.toglabels)
        self.boxFitpix.clicked.connect(self.togfitpix)
        self.boxResiduals.clicked.connect(self.togresiduals)
        self.loadParsButton.clicked.connect(self.showParFileDialog)


        fig=Figure()
        self.fig=fig
        self.initplot(fig)



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


    def initplot(self,fig,numchunks=8):
        wlen=len(spec.wavelength)/numchunks
        self.spls=[]
        if self.wave1==None:  waveidx1=0  # Default to plotting entire spectrum for now
        else: waveidx1=jbg.closest(wave,wave1)
        sg=jbg.subplotgrid(numchunks)
        for i in range(numchunks):
            self.spls.append(fig.add_subplot(sg[i][0],sg[i][1],sg[i][2]))
            #print sg[i][0],sg[i][1]
            #self.spls.append(plt.subplot(gs[sg[i][0],sg[i][1]]))
            pixs=range(waveidx1+i*wlen,waveidx1+(i+1)*wlen)
            self.spls[i].plot(self.wave[pixs],self.normflux[pixs],linestyle='steps')
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
        self.datamodel = LineParTableModel(self.fitpars,self.fiterrors,self.parinfo)
        self.tableView.setModel(self.datamodel)

    def fitlines(self):
        print 'Fitting line profile(s)...'
        self.fitpars,self.fiterrors=joebvpfit.joebvpfit(self.wave[cfg.fitidx],self.normflux[cfg.fitidx],self.normsig[cfg.fitidx],self.datamodel.fitpars,self.datamodel.parinfo)
        self.datamodel.updatedata(self.fitpars,self.fiterrors,self.parinfo,)
        self.updateplot()


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

    def showParFileDialog(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open line parameter file','.')
        self.initialpars(str(fname))

    def updateplot(self):
        if self.wave1==None:  waveidx1=0  # Default to plotting entire spectrum for now
        else: waveidx1=jbg.closest(self.wave,self.wave1)
        wlen=len(spec.wavelength)/self.numchunks
        for i,sp in enumerate(self.spls):
                sp.clear()
                prange=range(waveidx1+i*wlen,waveidx1+(i+1)*wlen)
    
                if ((len(self.fitpars[0])>0)):
    
                    sp.plot(self.wave,self.normflux,linestyle='steps')
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
                
    
                sp.plot(self.wave,self.normsig,linestyle='steps',color='red', lw=0.5)
                sp.plot(self.wave,-self.normsig,linestyle='steps',color='red', lw=0.5)
                sp.set_ylim(cfg.ylim)
                sp.set_xlim(self.wave[prange[0]],self.wave[prange[-1]])
                sp.set_xlabel('wavelength (A)', fontsize=cfg.xy_fontsize, labelpad=cfg.x_labelpad)
                sp.set_ylabel('normalized flux', fontsize=cfg.xy_fontsize, labelpad=cfg.y_labelpad)
        self.changefig(self.fig)

    def changefig(self, item):
        #text = str(item.text())
        self.rmmpl()
        self.addmpl(self.fig)
        
        
    def addfig(self, name, fig):
        self.fig_dict[name] = fig
        #self.mplfigs.addItem(name)
        
    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas,self.mplwindow,coordinates=True)
        self.mplvl.addWidget(self.toolbar)
        
    def rmmpl(self,):
        self.mplvl.removeWidget(self.canvas)
        self.canvas.close()
        self.mplvl.removeWidget(self.toolbar)
        self.toolbar.close()



if __name__ == '__main__':
        import sys
        from PyQt4 import QtGui
        import numpy as np
        from astropy.io import fits as pf
        from astropy.io import ascii
        from linetools.spectra.io import readspec


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

        ### Read in spectrum and list of lines to fit
        spec=readspec(filename)

        app = QtGui.QApplication(sys.argv)
        main = Main(spec,parfilename)
        #main.addfig('One plot',fig1)
        #main.addfig('Two Plots',fig2)
        #main.addfig('Pcolormesh',fig3)
        main.show()
        sys.exit(app.exec_())



