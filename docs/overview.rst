.. highlight:: rest

.. _joebvp:

******
joebvp
******

.. index:: joebvp

Overview
========

joebvp is a not-so-cleverly named collection of tools useful for generating and
fitting Voigt profiles of absorption lines, particularly in QSO spectra.  The
project was initiated and is maintained by Joe Burchett, although the
core fitting functionality of this software was based on IDL code by Joe
Meiring, which was in turn based partly on code written by Dan Welty.  Much
gratitude is due to these authors.

This software is optimized for handling large numbers of absorption lines/
systems simultaneously though these lines maybe lie at large separations
from each other across the spectrum.  However, simple systems are also very
easily fit using the code and user interface.  Also, the instrumental line
spread functions are served up by the *linetools* package, another open source
software project actively being developed.

The real strength of joebvp is in Voigt profile fitting, which can be accessed
in two modes: through a graphical user interface (GUI) or a 'batch mode'.

Using joebvp
============

As stated above, `joebvp` may be used for Voigt profile fitting using either a
GUI or a non-interactive batch mode.  Both modes of operation require a
continuum-fitted spectrum file as output by the `fit_continuum()` method of
`linetools.spectra` objects and a line parameter input file, which is
described below.  Also, both are accessed via the `VPmeasure` module, with each
called as follows.

First, assuming the `joebvp` package is contained within a directory set in
your PYTHONPATH environment variable, import the `VPmeasure` module::

    >from joebvp import VPmeasure

Then to access the GUI mode::

    >VPmeasure.go(specfile, paramfile)

where both `specfile` and `paramfile` should be strings.

Batch mode::

    >VPmeasure.batch_fit(specfile, parfilelist)

`parfilelist` is simply a file containing a list of input parameter file names,
one per line.


Line parameter input file
+++++++++++++++++++++++++

The starting parameter values for the lines to fit are defined in an input file
that is read in using the `astropy.io.ascii.read()` function.  The columns are
shown below as pipe-delimited, but the file need not be as long as the format
is readable by `astropy.io.ascii.read()`.  Each row in the table should define
one line, corresponding to a single atomic transition.

The columns of the `paramfile`::

    specfile|restwave|zsys|col|bval|vel|nflag|bflag|vflag|vlim1|vlim2|wobs1|wobs2|pix1|pix2|trans

Column definitions:

    specfile
        Name of the spectrum file used for this fit

    restwave
        Restframe wavelength of the transition in angstroms

    zsys
        Redshift of the absorption system that includes this line transition

    col
        Initial guess for the column density

    bval
        Initial guess for the Doppler b parameter

    vel
        Initial guess for the velocity centroid of the line with respect to zsys

    nflag, bflag, vflag
        Indicators to set behavior during fit for the column density, Doppler b,
        and velocity centroid.  A value of 0 for any parameter will allow it to
        freely vary during the fit.  A value of 1 will fix the parameter at the
        value set by the corresponding `col`, `bval`, or `vel` columns.  Values
        other than 0 or 1 tell the fitting routine to tie this parameter to
        that of another line, which will have the same flag value for the same
        parameter.

        For example, doublets such as O VI should likely have their column
        densities, b values, and velocity centroids tied together during the
        fit.  Therefore, the `nflag`, `bflag`, and `vflag` values should be
        set to the same value >1 in the rows for both the 1032 and 1038 \AA
        lines.

    vlim1, vlim2
        Lower and upper extents of the velocity range containing this line w/r
        to sys. These values will determine the range of pixels used to fit this
        line.  The fitter will minimize chi-squared evaluated for the union
        of pixel ranges set for all lines in the parameter input file.

    wobs1, wobs2
        Mininum and maximum observed wavelengths containing absorption from this line

    pix1, pix2
        First and last pixel in `specfile` showing absorption from this line

    trans
        Name of the transition, such as 'H I', 'Si II', 'C IV', etc.
