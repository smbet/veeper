import numpy as np
import joebvpfit
import makevoigt
from linetools.spectra.io import readspec
from linetools.spectra.xspectrum1d import XSpectrum1D
from astropy.table import Table,vstack
from astropy.io import ascii
import joebvp.cfg as cfg

def compose_model(spec,filelist,outfile):
    '''
    Generates full-model spectrum (XSpectrum1D) from a set of joebvp output files and writes out a file with the spectrum.

    Parameters
    ----------
    spec : string or XSpectrum1D
        The spectrum to be fitted with the input lines

    filelist : list of strings or str
        This should be a list containing the names of VP input files or a string referring to a file simply
        listing the input files.
        See joebvpfit.readpars for details of file format

    outfile : str
        Name of model spectrum output file

    '''

    ### Deal with alternate input types
    if isinstance(spec,str):
        specobj = readspec(spec)
    else:
        specobj = spec

    if isinstance(filelist, str):
        lstarr=np.genfromtxt(filelist,dtype=None)
        listofiles=lstarr.tolist()
    else:
        listofiles=filelist

    ### Load essentials from spectrum
    wave = specobj.wavelength.value
    normsig=specobj.sig.value/specobj.co.value
    cfg.wave = wave

    ### Loop through and concatenate all the parameter arrays
    tabs=[]
    for i, ff in enumerate(listofiles):
        tabs.append(ascii.read(ff))
    bigpartable=vstack(tabs)
    ascii.write(bigpartable,output='compiledVPinputs.dat',delimiter='|')  # write out compiled table
    fitpars, fiterrors, parinfo = joebvpfit.readpars('compiledVPinputs.dat')  # read this back in to load params
    cfg.fitidx = joebvpfit.fitpix(wave, fitpars)  # set the pixels in the line regions
    model = joebvpfit.voigtfunc(wave,fitpars)  # generate the Voigt profiles of all lines

    ### Instantiate XSpectrum1D object and write it out
    outspec=XSpectrum1D.from_tuple((wave,model,normsig))
    outspec.write_to_fits(outfile)

