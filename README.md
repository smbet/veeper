# joebvp

A Voigt profile fitter for CGM absorption lines. To run veeper on a previously-
identified absorption feature in a QSO spectrum, you will need to import
VPmeasure from joebvp and then call the batch_fit function.

e.g.:

  from joebvp import VPmeasure

  VPmeasure.batch_fit('specfile', ['componentlist'], filepath='./dirwherecomponentlistis/')

More detailed instructions, particularly for going from a pyigm igmguesses output to
a completed set of veeper-fitted components, can be found in pyigm_to_veeper_cookbook.txt
under docs.
