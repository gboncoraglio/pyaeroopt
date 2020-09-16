import numpy as np

from pyaeroopt.test.cantilever2994dof.pyao import frg, hpc
from pyaeroopt.test.cantilever2994dof.pyao.factory import aeros_hdm
from pyaeroopt.test.cantilever2994dof.pymt.gen_def_grad import gen_def_grad

F = gen_def_grad()

for i in range(F.shape[2]):
    Fi = F[:, :, i]
    aeros_hdm.execute(Fi.flatten('F'), None, hpc)

#n = 4
#v = np.linspace(-0.5, 0.5, n)
#for i in range(n):
#    for j in range(n):
#        for k in range(n):
#            for l in range(n):
#                for p in range(n):
#                    for q in range(n):
#                        x = np.array([v[i], v[j], v[k], v[l], v[p], v[q]])
#                        aeros_hdm.execute(x, None, hpc)
