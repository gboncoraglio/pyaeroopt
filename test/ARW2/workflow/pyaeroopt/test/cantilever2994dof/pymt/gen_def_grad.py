import sys, os, time, copy
sys.path.append('/home/mzahr/mjz_codes/pymortestbed/dev/')
sys.path.append('/home/mzahr/mjz_codes/pymortestbed/dev/'+
                'pymortestbed/app/struct_dynam_3d/src/_cython')

import numpy as np
import scipy.sparse

from pymortestbed.linalg.linsys      import LinearSolver, LinearSolverCholmod
from pymortestbed.linalg.decomp      import Svd
from pymortestbed.optimization.nlsys import NonlinearSolver
from pymortestbed.optimization.opt   import OPT

from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleVolume
from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleDefGrad
from pymortestbed.app.struct_dynam_3d.solid_fem import FilterSensitivity
import pymortestbed.app.struct_dynam_3d.src._cython.fem_struct as fem_struct
from pymortestbed.app.struct_dynam_3d.util import write_vtk

from fem_util import setup, construct_full_vec
from hdm_util import hdm_solv_fixmat 

from sklearn import cluster

def gen_def_grad():
    ## Linear/nonlinear solvers (sparse, dense)
    #linsolv = LinearSolverCholmod()
    linsolv = LinearSolver(scipy.sparse.linalg.spsolve)
    nlsolv  = NonlinearSolver(linsolv, 10, 1.0e-5, lambda x: np.max(np.abs(x)))
    
    # Geometry/Material/Force parameters
    which = 'Trestle'
    xlim = [0.0, 1.0]
    ylim = [0.0, 1.0]
    zlim = [0.0, 5.0]
    nelx = 10
    nely = 10
    nelz = 20
    forc =  -150000000.0
    mat00 = [7870.0, 9.695e10, 7.617e10]
    
    ( nsd, nel, nv, p, t_ptr, t, etype, ndof, ndbc, dbc_loc, dbc_val,
      nfbc, fbc_loc, fext0, fext, mat0, claw, mat_ptr, mat, nmat_per_el,
      eshape, nquad ) = setup(which.lower(), nelx, nely, nelz,
                              xlim, ylim, zlim, forc, mat00)
    
    # Make FEM structure
    msh  = fem_struct.MeshFem(p, t_ptr, t, etype, dbc_loc, dbc_val, dbc_val,
                              dbc_val)
    msh.create_sparsity(nmat_per_el, False, False, True)
    phys = fem_struct.PhysicsFem(claw, mat_ptr, mat)
    disc = fem_struct.DiscretizationFem(eshape, etype[0], nquad, nsd)
    
    u, success = hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, None, 1)
    u_nodal = construct_full_vec(msh.dbc_loc, u, msh.dbc_val)
    write_vtk('macro', msh, None, u_nodal, None, None, None, None, True)
    
    print success
    F = AssembleDefGrad(u, msh, disc)
    print F.shape
    return F

def cluster_save(F, nc=100):

    Fc = F.reshape((9, F.shape[2]), order='F')
    c, l, e = cluster.k_means(Fc.T, nc)
    np.savetxt('out/defgrad.nclust{0:d}.centers'.format(nc), c.T )
    for i in range(nc):
        np.savetxt('out/defgrad.nclust{0:d}.cluster{1:d}'.format(nc, i),
                   Fc[:, l==i])

def ccenter_load(nc=100):
        F=np.loadtxt('out/defgrad.nclust{0:d}.centers'.format(nc))
        return cluster_return(F, nc)

def cluster_load(nc=100):
    F = [[]]*nc
    for i in range(nc):
        F[i] = np.loadtxt('out/defgrad.nclust{0:d}.cluster{1:d}'.format(nc, i))
    return [cluster_return(f, nc) for f in F]

def cluster_return(F, nc):
    Fc = np.reshape(F.T, (3, 3, nc), order='F')
    return Fc

if __name__ == '__main__':
    F = gen_def_grad()
    cluster_save(F, 100)
