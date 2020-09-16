import sys, os, time, copy
sys.path.append('/Users/mzahr/Documents/codes/pymortestbed/dev/')
sys.path.append('/Users/mzahr/Documents/codes/pymortestbed/dev/'+
                'pymortestbed/app/struct_dynam_3d/src/_cython')

import numpy as np
import scipy.sparse.linalg

from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleVolume
from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleVolumeMatrix
from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleVolFracProjection
import pymortestbed.app.struct_dynam_3d.src._cython.fem_struct as fem_struct

from fem_util import *
from topopt_util import *

# HDM residual/jacobian
def hdm_resjac_fixmat(u, msh, phys, disc, fext, alpha, return_jac=None):

    Fint = AssembleVolume(u, msh, phys, disc, 'intforce')
    Fext = fext # + AssembleVolume(u, msh, phys, disc, 'body')
    if return_jac is None:
        return ( Fint-alpha*Fext )
    elif return_jac == 'mat':
        K = AssembleVolumeMatrix(u, msh, phys, disc, 'intforce_state')
        return ( Fint-alpha*Fext, K )
    elif return_jac == 'matvec':
        Kv = scipy.sparse.linalg.LinearOperator( (msh.ndof, msh.ndof),
                  lambda v: AssembleVolume(u, msh, phys, disc, 'intforce_state',
                                           None, None, v.reshape((-1, 1))))
        return ( Fint-alpha*Fext, Kv )

def hdm_resjac(u, p, msh, claw, mat_ptr, mat, disc, fext, alpha, simp=1.0,
               proj=True, beta=0.0, return_jac=None):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)
    # HDM res/jac, fixed material
    return hdm_resjac_fixmat(u, msh, phys, disc, fext, alpha, return_jac)

# HDM solve
def hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0=None, nstep=1):

    U0 = None

    # Initialize and check for best starting value
    U, res, zero_best = find_best_start(msh.ndof,
                      lambda U: hdm_resjac_fixmat(U, msh, phys, disc, fext, 1.0,
                                                  None), U0, nlsolv.norm)

    # Return if solution in U0
    if nlsolv.am_i_converged( res ):
        return U, True

    # Do not load step if not starting from 0
    if not zero_best:
        nstep = 1

    # Load stepping
    for i in range(nstep):
        alpha = float(i+1)/float(nstep)
        U, success = nlsolv.newton_raphson(lambda u: hdm_resjac_fixmat(
                                     u, msh, phys, disc, fext, alpha, 'mat'), U)
    return U, success

def hdm_solv(p, msh, claw, mat_ptr, mat, disc, fext, nlsolv, U0=None, nstep=1,
             simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    # HDM solve, fixed material
    return hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0=None, nstep=1)


# HDM adjoint
def hdm_adjoint_fixmat_nosolv(u, func, msh, phys, disc, fext, rho_el, proj=True,
                              beta=0.0):

    # Derivative of residual w.r.t. state
    drdu = AssembleVolumeMatrix(u, msh, phys, disc, 'intforce_state')

    # Evaluate derivatives of functional
    f, dfdu, dfdr = func(u, rho_el, msh, phys, disc, fext, True)
    dfdp = projection_chain_rule(dfdr, msh, p, proj, beta)

    # Adjoint method: df = dfdp - (inv(drdu')*dfdu')'*drdp
    lam = linsolv.solve(dfdu, drdu.transpose()).reshape((-1,1))
    return lam

def hdm_adjoint_nosolv(u, p, func, msh, claw, mat_ptr, mat, disc, fext,
                       simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    return hdm_adjoint_fixmat_nosolv(u, func, msh, phys, disc, fext, rho_el,
                                     proj, beta)

def hdm_adjoint(p, func, msh, claw, mat_ptr, mat, disc, fext, nlsolv, linsolv,
                U0=None, nstep=1, simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    # Solve for displacement, update solution
    u, success = hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0, nstep)
    U0[:, -1] = u

    return hdm_adjoint_fixmat_nosolv(u, func, msh, phys, disc, fext, rho_el,
                                     proj, beta)

# HDM functional/gradient
def hdm_func_nosolv_fixmat(u, p, func, msh, phys, disc, fext):

    # Evaluate functional
    f = func(u, None, msh, phys, disc, fext, False)
    return f

def hdm_func_nosolv(u, p, func, msh, claw, mat_ptr, mat, disc, fext,
                    simp=1.0, proj=True, beta=0.0, inf=1.0e20):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys= create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj, beta)[0]

    # Evaluate functional
    return hdm_func_nosolv_fixmat(u, p, func, msh, phys, disc, fext, inf)
    
def hdm_func_fixmat(p, func, msh, phys, disc, fext, nlsolv, U0=None, nstep=1,
                    inf=1.0e20):

    # Solve for displacement, update solution
    u, success = hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0, nstep)
    U0[:, -1] = u

    # Evaluate functional
    if success:
        return hdm_func_nosolv_fixmat(u, p, func, msh, phys, disc, fext)
    else:
        return inf

def hdm_func(p, func, msh, claw, mat_ptr, mat, disc, fext, nlsolv,
                   U0=None, nstep=1, simp=1.0, proj=True, beta=0.0, inf=1.0e20):
    
    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys= create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj, beta)[0]

    # Evaluate functional
    return hdm_func_fixmat(p, func, msh, phys, disc, fext, nlsolv, U0, nstep)

def hdm_grad_nosolv_noadj_fixmat(u, lam, p, func, msh, phys, disc, fext, mat,
                                 rho_el, simp=1.0, proj=True, beta=0.0):
    # Evaluate functional
    f, dfdu, dfdr = func(u, rho_el, msh, phys, disc, fext, True)
    dfdp = projection_chain_rule(dfdr, msh, p, proj, beta)

    # Adjoint method: df = dfdp - (inv(drdu')*dfdu')'*drdp
    df0 = AssembleVolume(u, msh, phys, disc, 'intforce_mat', None, lam)
    df = dfdp - projection_chain_rule_frommat(df0, msh, mat, p, simp,
                                              proj, beta)
    return df

def hdm_grad_nosolv_fixmat(u, p, func, msh, phys, disc, fext, mat, rho_el,
                           linsolv, simp=1.0, proj=True, beta=0.0):

    # Derivative of residual w.r.t. state
    drdu = AssembleVolumeMatrix(u, msh, phys, disc, 'intforce_state')

    # Evaluate derivatives of functional
    f, dfdu, dfdr = func(u, rho_el, msh, phys, disc, fext, True)

    # Adjoint method: df = dfdp - (inv(drdu')*dfdu')'*drdp
    lam = linsolv.solve(dfdu, drdu.transpose()).reshape((-1,1))
    return  hdm_grad_nosolv_noadj_fixmat(u, lam, p, func, msh, phys, disc, fext,
                                         mat, rho_el, simp, proj, beta)

def hdm_grad_noadj_fixmat(lam, p, func, msh, phys, disc, fext, mat, rho_el,
                          nlsolv, U0=None, nstep=1, simp=1.0, proj=True,
                          beta=0.0):

    # Solve for displacement, update solution
    u, success = hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0, nstep)
    U0[:, -1] = u

    # Adjoint method: df = dfdp - (inv(drdu')*dfdu')'*drdp
    return  hdm_grad_nosolv_noadj_fixmat(u, lam, p, func, msh, phys, disc, fext,
                                         mat, rho_el, simp, proj, beta)
    
def hdm_grad_nosolv_noadj(u, lam, p, func, msh, claw, mat_ptr, mat, disc, fext,
                          simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    # Evaluate gradient
    return hdm_grad_nosolv_noadj_fixmat(u, lam, p, func, msh, phys, disc, fext,
                                        mat, rho_el, simp, proj, beta)

def hdm_grad_nosolv(u, p, func, msh, claw, mat_ptr, mat, disc, fext, linsolv,
                    simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    # Evaluate gradient
    return hdm_grad_nosolv_fixmat(u, p, msh, phys, disc, fext, mat, rho_el,
                                  linsolv, simp, proj, beta)

def hdm_grad_noadj(lam, p, func, msh, claw, mat_ptr, mat, disc, fext, nlsolv,
                   U0=None, nstep=1, simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    # Solve for displacement, update solution
    u, success = hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0, nstep)
    U0[:, -1] = u

    # Adjoint method: df = dfdp - (inv(drdu')*dfdu')'*drdp
    return  hdm_grad_nosolv_noadj_fixmat(u, lam, p, func, msh, phys, disc, fext,
                                         mat, rho_el, simp, proj, beta)

def hdm_grad(p, func, msh, claw, mat_ptr, mat, disc, fext, nlsolv, linsolv,
             U0=None, nstep=1, simp=1.0, proj=True, beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, simp, proj,
                                  beta)

    # Solve for displacement, update solution
    u, success = hdm_solv_fixmat(msh, phys, disc, fext, nlsolv, U0, nstep)
    U0[:, -1] = u

    # Derivative of residual w.r.t. state
    drdu = AssembleVolumeMatrix(u, msh, phys, disc, 'intforce_state')

    # Evaluate derivatives of functional
    f, dfdu, dfdr = func(u, rho_el, msh, phys, disc, fext, True)
    dfdp = projection_chain_rule(dfdr, msh, p, proj, beta)

    # Adjoint method: df = dfdp - (inv(drdu')*dfdu')'*drdp
    lam = linsolv.solve(dfdu, drdu.transpose()).reshape((-1,1))
    return  hdm_grad_nosolv_noadj_fixmat(u, lam, p, func, msh, phys, disc, fext,
                                         mat, rho_el, simp, proj, beta)

# HDM functional/gradient (state-indep)
def hdm_func_nostate(p, func, msh, claw, mat_ptr, mat, disc, fext, proj=True,
                     beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, 1.0, proj,
                                  beta)
    return hdm_func_nostate_fixmat(func, msh, phys, disc, fext, rho_el)

def hdm_func_nostate_fixmat(func, msh, phys, disc, fext, rho_el):

    # Evaluate functional
    f = func(rho_el, msh, phys, disc, fext, False)
    return f

def hdm_grad_nostate(p, func, msh, claw, mat_ptr, mat, disc, fext, proj=True,
                     beta=0.0):

    # Physics structure based on value of parameter p
    matn = np.empty_like(mat, dtype=float, order='F')
    phys, rho_el = create_physics(p, msh, claw, mat_ptr, mat, matn, 1.0, proj,
                                  beta)

    return hdm_grad_nostate_fixmat(p, func, msh, phys, disc, fext, rho_el, proj,
                                   beta)
    
def hdm_grad_nostate_fixmat(p, func, msh, phys, disc, fext, rho_el, proj=True,
                            beta=0.0):

    # Evaluate derivatives of functional
    dfdr = func(rho_el, msh, phys, disc, fext, True)[1]
    dfdp = projection_chain_rule(dfdr, msh, p, proj, beta)
    return dfdp
