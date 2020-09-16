import numpy as np

from pyaeroopt.util.aeros_util import read_aeros_mesh
from pyaeroopt.test.cantilever2994dof import src_dir
from pyaeroopt.test.cantilever2994dof.pyao import frg

def extract_boundary_nodes(t):

    perm = [[0, 1, 2], [1, 2, 3], [0, 1, 3], [0, 2, 3]]
    S = set([])
    for i, ti in enumerate(t):
        for j in range(4):
            fij = set(ti[perm[j]])

            bndy = True
            for k, tk in enumerate(t):
                if k == i: continue
                if fij.issubset(set(tk)):
                    bndy=False
                    break
            if bndy:
                S = set.union(S, fij)
    return np.array(list(S))

if __name__ == '__main__':

    p, t_ptr, t, etype, dbc_loc, dbc_val, fext = read_aeros_mesh(frg.top)
    bndy = extract_boundary_nodes(t.reshape((len(t_ptr)-1, 4))-1)
    p_bndy = np.hstack((bndy[:, None]+1, p[:, bndy].T))
    np.savetxt(src_dir+'bndy.nodeset', p_bndy)
