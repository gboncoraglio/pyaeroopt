import sys, os

from pyaeroopt.interface.aeros import AerosInputBlock

# Problem Instance
from pyaeroopt.test.cantilever2994dof.pyao import frg

control=AerosInputBlock('CONTROL', None, 1, '"nodeset"', '"elemset"')
static=AerosInputBlock('STATIC', 'sparse')
dynam=AerosInputBlock('DYNAMICS', ['time', 0.0, 0.1, 1.0],
                                  ['newmark'],
                                  ['mech', 0.0, 0.0, 0.0, 0.0])

nonlin=AerosInputBlock('NONLINEAR',
                       ['maxit', 100],
                       ['nltol', 1.0e-8],
                       ['dlambda', 1.00, 1.00])

robc=AerosInputBlock('ROBC',
                     ['snapfi', None],
                     ['mnorma', 1])

readmode=AerosInputBlock('READMODE', ['READMODE', None, 1])

output_hdm=AerosInputBlock('OUTPUT',
                          ['gdisplac', 20, 16, None, 1],
                          ['stressvm', 20, 16, None, 1],
                          ['statevct', 20, 16, None, 1])

output_svd=AerosInputBlock('OUTPUT',
                           ['robdataf', None, 1])

output_rom=AerosInputBlock('OUTPUT',
                          ['gdisplac', 20, 16, None, 1],
                          ['stressvm', 20, 16, None, 1])

topo=AerosInputBlock('INCLUDE_TOPO', ['INCLUDE', frg.top])
disp=AerosInputBlock('INCLUDE_DISP', ['INCLUDE', None])
forc=AerosInputBlock('INCLUDE_FORC', ['INCLUDE', None])
