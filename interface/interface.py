"""This module exports classes to interface with external codes.

Exported Classes:
---------------------
    * CodeInterface: an object to interface with a command-line
                     input-file-based codes.
"""

import os
from pyaeroopt.util.hpc_util import execute_code, batch_slurm_pbs

class CodeInterface(object):

    """
    A class to facilitate interfacing with external codes.

    This class can be used to interface with any code that is
        * callable from the command line
        * takes an input file (provided as a command line argument) for input
        * returns its output in an output file

    Data Members
    ------------
    db : SQLAlchemy object
        Database object
    bin : str
        Path to executable
    desc: python object (string, list, tuple, dict, etc)
        Object specifying class instance
    exec_dir: str
        Path (either absolute or relative to the directory that Python will
        launch from) to the directory that the code should be run in. (E.g.
        sdesign must be run in the same directory as its input file)
    """

    def __init__(self, **kwargs):

        # Extract anticipated input
        self.db = kwargs.pop('db')   if 'db'   in kwargs else None
        self.bin = kwargs.pop('bin')  if 'bin'  in kwargs else None
        self.desc = kwargs.pop('desc') if 'desc' in kwargs else None
        self.exec_dir = kwargs.pop('exec_dir') if 'exec_dir' in kwargs else './'
        self.infile = None

        # Extract remaining input
        for kwarg in kwargs:
            setattr(self, kwarg, kwargs[kwarg])

    def set_database(self, db):
        self.db = db

    def initialize_database(self, fname):
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        engine = create_engine("sqlite:///{0:s}".format(fname), echo=True)
        self.db = sessionmaker(bind=engine)

    def create_input_file(self, p, desc_ext=None, db=None):
        """Problem-specific: must be defined by user."""
        pass

    def check_database(self, desc_ext):
        """Check whether this instance of the problem is in the database."""
        exist = False
        if self.db is not None:
            exist = self.db.does_exist(self.infile)
        if exist:
            # If instance already exists, update database with external desc
            self.db.update_entry(self.infile, desc_ext)
        return exist

    def writeInputFile(self):
        current_dir = os.getcwd()
        os.chdir(self.exec_dir)
        self.infile.write()
        os.chdir(current_dir)

    def execute(self, p=None, desc_ext=None, hpc=None, make_call=True, name="myJob"):
        """
        Check that instance of problem has not been run (if database specified),
        create input file, call executable, and add instance to database (if
        applicable)
        """

        # Create input file and check if instance contained in database
        try:
            self.create_input_file(p, desc_ext, self.db)
        except TypeError:
            self.create_input_file(p)
        exist = self.check_database(desc_ext)

        if not exist:
            # If instance does not exist, write input file, call executable,
            # and add to database
            self.writeInputFile()

            # for coupled simulations, self.infile and self.log will be arrays
            if hasattr(self.infile, "__iter__"):
                infileArray = [x.fname for x in self.infile]
                logFile     = self.infile[0].log + "-" + self.infile[1].log
            else:
                infileArray = self.infile.fname
                logFile = self.infile.log
            if hpc.batch:
                batch_slurm_pbs(name, hpc, self.exec_dir, self.bin, infileArray, logFile)
            else:
                exec_str = hpc.execute_str(self.bin, infileArray)
                current_dir = os.getcwd()
                os.chdir(self.exec_dir)
                execute_code(exec_str, logFile, make_call)
                os.chdir(current_dir)
            if self.db is not None: self.db.add_entry(p, self.desc, self.infile)
