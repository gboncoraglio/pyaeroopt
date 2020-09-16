"""Declare problem configuration constants.

This module declares constants that do not define the problem itself,
(which is defined in factory.py), but rather the way in which the problem
will be run (on which machine, in which folder, etc).
"""

from pyaeroopt.interface import Frg, Hpc

ROM_file = 'ROMs_ready.mat'
ROM_interpolation_folder = '/home/aspenser/ROM_Interpolation'


# FRG utility class (default binary locations, $SOWER, $XP2EXO, etc)
ARW2 = Frg(top='AeroF-Files/ARW2.top',
           feminp="AeroS-Files/ARW2.input",
           geom_pre='binaries/ARW2')

# HPC utility class
hpc = Hpc(machine='independence', mpi='mpiexec',
          batch=False,
          nproc=32,
          bg=False)

hpc1 = Hpc(machine='independence',
           nompi=True,
           batch=False,
           bg=False)

hpcCoupled = Hpc(machine='independence', mpi='mpiexec',
                 batch=False,
                 nproc=[32, 1],   # [aerof processors, aeros processors]
                 bg=False)
