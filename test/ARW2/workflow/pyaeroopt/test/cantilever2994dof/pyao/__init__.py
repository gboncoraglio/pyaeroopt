from pyaeroopt.interface import Frg, Hpc
from pyaeroopt.test.cantilever2994dof import src_dir

# FRG utility class (default binary locations, $SOWER, $XP2EXO, etc)
frg = Frg(top=src_dir+'top') 

# HPC utility class
hpc = Hpc(machine='independence', mpi='mpiexec', # use modules on independence
          batch=False, bg=False)                 # so mpiexec is sufficient to
                                                 # use appropriate version of
                                                 # MPI
