"""Define objects that set the design point for a simulation.

Exports the following three subclasses of CodeInterface:
    * StructureDeformer: a class that runs Sdesign to move the structural
            mesh and compute shape derivatives.
    * FluidSkinDeformer: a class that runs SDesign to move the portion of the
            fluid mesh that lies against the wing surface and compute shape
            derivatives
    * MaterialFileChanger: a class that writes a material file that reflects
            the change in the three structural thickness variables in this
            problem.
"""

import os
from copy import deepcopy

from pyaeroopt.interface.sdesign import Sdesign, SdesignInputFile
from pyaeroopt.interface.interface import CodeInterface
from pyaeroopt.test.ARW2.pyao.sdesign_blk import *
from pyaeroopt.test.ARW2.workflow import ARW2

class StructureDeformer(Sdesign):

    """
    A class to deform the structural mesh and compute its shape derivatives

    Usage:
        StructureDeformer() creates the mesh deformer object
        StructureDeformer.execute(p=[...]) deforms meshes and sets up material files
            for a problem with the parameter vector p
    """

    def __init__(self, **kwargs):
        super(StructureDeformer, self).__init__(**kwargs)
        self.infile = None

    def create_input_file(self, p, desc_ext=None, db=None):
        """Create an SdesignInputfile to move the structural mesh."""
        # Define the prefix for this problem's input and output files
        def append_prefix(file_extension):
            """Appends 'ARW2.deformStructure' prefix to filename strings."""
            prefix = 'ARW2.deformStructure'
            return "{0:s}.{1:s}".format(prefix, file_extension)

        # Create the input file object
        fname = append_prefix('inp')
        log = append_prefix('log')

        # absvar is a function that returns an input block and takes in the
        # parameter vector p.  All other input blocks are static
        self.infile = SdesignInputFile(fname,
                                       deepcopy([outstruct, define, nodes,
                                                 edges, patch, volume, dsgvar,
                                                 absvar(p), link,
                                                 femeshstruct]),
                                       log)

class FluidSkinDeformer(Sdesign):

    """
    A class to deform the fluid mesh that lies on the ARW2's skin.

    Usage:
        FluidSkinDeformer() creates the mesh deformer object
        FluidSkinDeformer.execute(p=[...]) deforms meshes for a problem with
            the parameter vector p
    """

    def __init__(self, **kwargs):
        super(FluidSkinDeformer, self).__init__(**kwargs)
        self.infile = None

    def create_input_file(self, p, desc_ext=None, db=None):
        """Create an SdesignInputFile to move the fluid mesh."""
        ## Define the prefix for this problem's input and output files
        def append_prefix(file_extension):
            """Appends 'ARW2.deformFluidSkin' prefix to filename strings."""
            prefix = 'ARW2.deformFluidSkin'
            return "{0:s}.{1:s}".format(prefix, file_extension)

        # Create the input file object
        fname = append_prefix('inp')
        log = append_prefix('log')

        # absvar is a function that returns an input block but takes in the
        # parameter vector p.  All other input blocks are static.
        self.infile = SdesignInputFile(fname,
                                       deepcopy([outfluid, define, nodes,
                                                 edges, patch, volume, dsgvar,
                                                 absvar(p), link,
                                                 femeshfluid]),
                                       log)

    def execute(self, p, desc_ext=None, hpc=None, make_call=True):
        """
        Check that instance of problem has not been run (if database specified),
        create input file, call executable, and add instance to database (if
        applicable)
        """

        # Execute sdesign to find the deformation of the fluid mesh on the wing
        super(FluidSkinDeformer, self).execute(p, desc_ext, hpc, make_call)

        # Use sower to prepare the wing skin deformation for Aero-F to process
        ARW2.sower_fluid_mesh_motion('Sdesign/fluidposition',
                                     'Fdata/fluidposition.idisp')

class MaterialFileChanger(CodeInterface):

    """
    A class to change the material file for a new design point.

    Usage:
        MaterialFileChanger() creates the mesh deformer object
        MaterialFileChanger.execute(p=[...]) sets up material files
            for a problem with the parameter vector p
    """

    def __init__(self, **kwargs):
        super(MaterialFileChanger, self).__init__(**kwargs)

    def execute(self, p, desc_ext=None, hpc=None, make_call=True):
        """Write a material file that reflects the new design point."""

        # Write the material file in this object's execution directory
        fname = "materialfile"
        current_dir = os.getcwd()
        os.chdir(self.exec_dir)
        f = open(fname, 'w')
        materialvector = [0.1*x for x in p[3:]]
        f.write(materialfile(materialvector))
        f.close()
        os.chdir(current_dir)

        # Write the group file in this object's execution directory
        fname = "groupfile"
        current_dir = os.getcwd()
        os.chdir(self.exec_dir)
        f = open(fname, 'w')
        f.write('GROUP\n')
        for attribute, group in groupDict.items():
            f.write(str(attribute) + '    ' + str(group) + '\n')
        f.write('*')
        f.close()
        os.chdir(current_dir)
