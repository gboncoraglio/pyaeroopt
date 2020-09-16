#TODO: Use db (database) to CREATE filename!
import numpy as np
from copy import deepcopy

from pyaeroopt.interface.aerof import Aerof, AerofInputFile
from pyaeroopt.interface import CodeInterface


class Aerofs(CodeInterface):
    """
    A class that allows an Aerof object and an Aeros object to be bundled for
    joint execution

    Usage: Aerofs(Aerof(input params), Aeros(input params))

    The first argument is an Aerof object, the second argument is an Aeros object
    """
    def __init__(self, aerof, aeros, **kwargs):

        super(Aerofs, self).__init__(**kwargs)

        self.aerof = aerof
        self.aeros = aeros
        self.db = None

    def writeInputFile(self):
        self.aerof.infile.write()
        self.aeros.infile.write()

    def create_input_file(self, p, desc_ext=None, db=None):
        self.aerof.create_input_file(p, desc_ext=desc_ext, db=db)
        self.aeros.create_input_file(p, desc_ext=desc_ext, db=db)
        self.infile = [self.aerof.infile, self.aeros.infile]
        self.bin    = [self.aerof.bin,    self.aeros.bin]





