"""Export a class to run the Aero-S portion of an aeroelastic sensitivity
analysis.
"""

from copy import deepcopy

from pyaeroopt.interface.aeros import Aeros, AerosInputFile
from pyaeroopt.test.ARW2.pyao.aeros_blk import *

class AerosAeroelasticSensitivity(Aeros):
    """A class to run the aeros part of an aeroelastic sensitivity calculation.
    """

    def __init__(self, **kwargs):
        super(AerosAeroelasticSensitivity, self).__init__(**kwargs)

    def create_input_file(self, p, desc_ext=None, db=None):
        """Creates an AerosInputFile for the aeroelastic sensitivity problem.

        Parameters:
        --------------
        p: 1D array
            A vector containing the design variables to analyze the problem at
        """

       # Define the prefix for this problem's input and output files
        def append_prefix(file_extension):
            """Appends 'aeros_sens' prefix to filename strings."""
            prefix = 'aeros_sens'
            return "{0:s}.{1:s}".format(prefix, file_extension)

        # Create the input file object
        fname = append_prefix('inp')
        log = append_prefix('log')
        self.infile = AerosInputFile(fname,
                                     deepcopy([controlSens, static, qstatic,
                                               aero, renum, sens, grav,
                                               outputSens, sdes, top, cframe,
                                               disp, dim, eframe, comp, mat,
                                               att, group]),
                                     log)

        self.infile.AERO.MPP = None
        self.infile.AERO.A6 = ' '
        self.infile.AERO.READMODES = None
