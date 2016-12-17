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

    ### Load essentials from spectrum
    wave = specobj.wavelength.value
    normsig=specobj.sig.value/specobj.co.value
    cfg.wave = wave

    ### Concatenate all the parameter tables
    concatenate_line_tables(filelist,outtablefile='compiledVPinputs.dat')

    ### Make the model!
    fitpars, fiterrors, parinfo = joebvpfit.readpars('compiledVPinputs.dat')  # read this back in to load params
    cfg.fitidx = joebvpfit.fitpix(wave, fitpars)  # set the pixels in the line regions
    model = joebvpfit.voigtfunc(wave,fitpars)  # generate the Voigt profiles of all lines

    ### Instantiate XSpectrum1D object and write it out
    outspec=XSpectrum1D.from_tuple((wave,model,normsig))
    outspec.write_to_fits(outfile)

def concatenate_line_tables(filelist,outtablefile='compiledVPinputs.dat'):
    '''
    Compiles the output from several fitting runs into a single table

    Parameters
    ----------
    filelist : list of strings or str
        This should be a list containing the names of VP input files or a string referring to a file simply
        listing the input files.
        See joebvpfit.readpars for details of file format

    outtablefile : str
        Name of compiled model parameter output file

    '''

    if isinstance(filelist, str):
        lstarr=np.genfromtxt(filelist,dtype=None)
        listofiles=lstarr.tolist()
    else:
        listofiles=filelist

    tabs = []
    for i, ff in enumerate(listofiles):
        tabs.append(ascii.read(ff))
    bigpartable = vstack(tabs)
    ascii.write(bigpartable, output=outtablefile, delimiter='|')  # write out compiled table

def abslines_from_VPfile(parfile,ra=None,dec=None):
    '''
    Takes a joebvp parameter file and builds a list of linetools AbsLines from the measurements therein.

    Parameters
    ----------
    parfile : str
        Name of the parameter file in the joebvp format
    ra : float
        Right Ascension of the QSO in decimal degrees
    dec : float
        Declination of the QSO in decimal degress

    Returns
    -------
    abslinelist: list
        List of AbsLine objects
    '''
    from linetools.spectralline import AbsLine
    import astropy.units as u
    fitpars, fiterrors, parinfo = joebvpfit.readpars(parfile)
    linetab = ascii.read(parfile)
    abslinelist = []
    for i,row in enumerate(linetab):
        line=AbsLine(row['restwave']*u.AA, z=row['zsys'])
        vlims=[row['vlim1'],row['vlim2']]*u.km/u.s
        line.limits.set(vlims)
        line.attrib['logN'] = row['col']
        line.attrib['sig_N'] = row['sigcol']
        line.attrib['b'] = row['bval']
        line.attrib['sig_b'] = row['sigbval']
        abslinelist.append(line)
    return abslinelist



