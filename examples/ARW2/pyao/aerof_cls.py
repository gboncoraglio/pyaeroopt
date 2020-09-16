"""Export a class to run the Aero-F portion of an aeroelastic sensitivity
analysis.
"""


from copy import deepcopy

from pyaeroopt.interface.aerof import Aerof, AerofInputFile
from pyaeroopt.test.ARW2.pyao.aerof_blk import *
from pyaeroopt.test.ARW2.workflow import ARW2

class AerofAeroelasticSensitivity(Aerof):
    """
    p : numpy array
        Shape parameters
    desc : list
        [type]
    desc_ext : list
        [multsoln]
    """
    def __init__(self, **kwargs):
        super(AerofAeroelasticSensitivity, self).__init__(**kwargs)


    def create_input_file(self, p, desc_ext=None, db=None):
        """Create an AerofInputFile for an aeroelastic sensitivity problem.

        Parameters:
        ------------
        p:  1D array
            An array of design variables to analyze the system at.
        """

        # Define the prefix for this problem's input and output files
        def append_prefix(file_extension):
            """Appends 'aerof_sens' prefix to filename strings."""
            prefix = 'aerof_sens'
            return "{0:s}.{1:s}".format(prefix, file_extension)

        # Create the input file object
        fname = append_prefix('inp')
        log = append_prefix('log')
        self.infile = AerofInputFile(fname,
                                     deepcopy([prob, inpu, outp, sensAnal,
                                               equa, refState, bounCond, spac,
                                               time, meshMoti]),
                                     log)

        # Fill in any input file fields that differ from the default
        self.infile.Input.GeometryPrefix = ARW2.geom_pre

        self.infile.Problem.Type = 'SteadyAeroelasticSensitivityAnalysis'
        self.infile.Input.Matcher = 'binaries/ARW2.match'
        self.infile.Input.InitialWallDisplacement  = 'Fdata/fluidposition.idisp'

        self.infile.Output.Postpro.Prefix                  = 'Fresults/'
        self.infile.Output.Postpro.LiftandDragSensitivity  = "../out/ARW2.sensitivity.lift.drag"
        self.infile.Output.Postpro.Pressure                = "ARW2.pressure"
        self.infile.Output.Postpro.Displacement            = "ARW2.disp"
        self.infile.Output.Postpro.Frequency               = 2000

        self.infile.SensitivityAnalysis.SensitivityFSI = 'On'
        self.infile.SensitivityAnalysis.AdaptiveEpsFSI = 'On'
        self.infile.SensitivityAnalysis.LinearSolver.Eps = 1.0e-7
