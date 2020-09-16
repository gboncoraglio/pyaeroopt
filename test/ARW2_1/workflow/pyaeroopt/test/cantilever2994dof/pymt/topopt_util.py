import sys, copy
import numpy as np

sys.path.append('/Users/mzahr/Documents/codes/pymortestbed/dev/')
sys.path.append('/Users/mzahr/Documents/codes/pymortestbed/dev/'+
                'pymortestbed/app/struct_dynam_3d/src/_cython')
from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleVolume
from pymortestbed.app.struct_dynam_3d.solid_fem import AssembleVolFracProjection
import pymortestbed.app.struct_dynam_3d.src._cython.fem_struct as fem_struct
from pymortestbed.optimization.opt import OPT

def solve_topopt(nvar, p0, pl, pu, obj, grad, constr, dconstr, cl, cu,
                 solver='pyopt:SNOPT', options=None, minf=-1.0e20, inf=1.0e20):

    opt = OPT(minf=-1.0e20, inf=1.0e20)
    opt.addVariables(nvar, p0, pl, pu)
    opt.addObjective(obj, grad)
    for i in range( len(constr) ):
        opt.addNonlinConstraints(1, constr[i], dconstr[i],
                                 cl[i]*np.ones(1), cu[i]*np.ones(1))
    xStar, fStar = opt.optimize('pyopt:SNOPT', 'analytical', options=options)
    return xStar    

def solve_topopt_mc(nvar, p0, pl, pu, obj, grad, constr, dconstr, cl, cu,
                 solver='pyopt:SNOPT', options=None, minf=-1.0e20, inf=1.0e20):

    opt = OPT(minf=-1.0e20, inf=1.0e20)
    opt.addVariables(nvar, p0, pl, pu)
    opt.addObjective(obj, grad)
    opt.addNonlinConstraints(2, constr, dconstr, cl, cu)
    xStar, fStar = opt.optimize('pyopt:SNOPT', 'analytical', options=options)
    return xStar

def create_physics(p, msh, claw, mat_ptr, mat, new_mat, simp=1.0, proj=False,
                   beta=0.0):
    # new_mat is a buffer that must exist in the calling function since
    # phys only stores a POINTER to it; if it does not, assemble will either
    # crash (seg fault) or return zero vectors/matrices (assume matprop = 0)

    rho_el  = get_elem_volfrac(p, msh, proj, beta)
    new_mat[:] = mat*np.vstack((rho_el**simp, rho_el**simp, rho_el**simp)
                                                                  ).flatten('F')
    phys = fem_struct.PhysicsFem(claw, mat_ptr, new_mat)
    return phys, rho_el

def find_best_start(ndof, res_func, x0, norm):

    u = np.zeros(ndof, dtype=float, order='F')
    zero_best = True
    res = norm( res_func(u) )
    if ndof < 10:
        return u, res, zero_best
    if x0 is not None:
        x0.reshape((-1, 1))
        for i in range( x0.shape[1] ):
            R = res_func(x0[:, i])
            nres = norm(R)
            if nres < res:
                u = x0[:, i]
                res = nres
                zero_best = False
    return u, res, zero_best

def get_elem_volfrac(mu, msh, proj, beta):

    if proj:
        rho_el = AssembleVolFracProjection(msh, mu, beta)
    else:
        rho_el = mu
    return rho_el

def projection_chain_rule_frommat(dmat, msh, mat, mu, simp, proj, beta):

    rho_el = get_elem_volfrac(mu, msh, proj, beta)
    dmat_new = mat*np.vstack((simp*(rho_el**(simp-1)),
                              simp*(rho_el**(simp-1)),
                              simp*(rho_el**(simp-1)))).flatten('F')
    drho = ( dmat[::3]*dmat_new[::3] + dmat[1::3]*dmat_new[1::3] +
             dmat[2::3]*dmat_new[2::3] )
    dmu = projection_chain_rule(drho, msh, mu, proj, beta)
    return dmu

def projection_chain_rule(drho, msh, mu, proj, beta):

    if proj:
        dmu = AssembleVolFracProjection(msh, mu, beta, drho.reshape((-1, 1)))[1]
    else:
        dmu = drho
    return dmu

def elbasis_from_pbasis(p_basis, msh, proj=True, beta=0.0):

    npart = p_basis.shape[1]
    if proj:
        el_basis = np.zeros((msh.nel, npart), dtype=float, order='F')
        for n in range(npart):
            el_basis[:, n] = AssembleVolFracProjection(msh, p_basis[:, n],
                                                       beta)
    else:
        el_basis = copy.copy(p_basis)
    return el_basis

def compliance(u, rho, msh, phys, disc, fext, return_derivs=False):

    c = np.dot(fext, u)
    if not return_derivs: return c

    dcdu = fext[:]
    dcdr = np.zeros(rho.size, dtype=float, order='F')
    return c, dcdu, dcdr

def rcompliance(u_r, p_r, ubar_Fext, Vt_Fext, return_derivs=False):

    cr = ubar_Fext+np.dot(u_r, Vt_Fext)
    if not return_derivs: return cr

    dcrdur = Vt_Fext[:]
    dcrdrr = np.zeros(p_r.size, dtype=float, order='F')
    return cr, dcrdur, dcrdrr

def volume(rho, msh, phys, disc, fext, return_derivs=False):

    if rho is None: rho = np.ones((msh.nel, 1), dtype=float, order='F')
    if rho.ndim == 1: rho = rho.reshape((-1, 1))
    V = AssembleVolume(np.zeros(msh.ndof, dtype=float, order='F'),
                       msh, phys, disc, 'volume', None, rho)
    if not return_derivs: return V

    dVdr = AssembleVolume(np.zeros(msh.ndof, dtype=float, order='F'),
                          msh, phys, disc, 'volume')
    return V, dVdr 

def rvolume(p_r, Vr, return_derivs=False):

    if not return_derivs: return np.dot(p_r, Vr)
    return np.dot(p_r, Vr), Vr[np.newaxis, :]

def volume_vec(msh, phys, disc):

    V = AssembleVolume(np.zeros(msh.ndof, dtype=float, order='F'), 
                       msh, phys, disc, 'volume')
    return ( V )

def rvolume_vec(msh, phys, disc, p_basis, proj, beta=0.0):

    el_basis =  elbasis_from_pbasis(p_basis, msh, proj, beta)
    Vr = AssembleVolume(np.zeros(msh.ndof, dtype=float, order='F'), 
                        msh, phys, disc, 'volume', None, el_basis)
    return ( Vr )

def create_param_basis(part, nopt, npart):

    k = 0
    param_basis = np.zeros((nopt, npart), dtype=float, order='F')
    for ip in part:
        param_basis[ip, k] = 1.0
        k += 1
    return param_basis

def partition_param_vec_hexregion(npartx, nparty, npartz, msh,
                                  nelx, nely, nelz, proj=False):

    # Convert element numbers in each direction to global element numbers
    nnx = nelx+1; nny = nely+1; nnz = nelz+1;
    if proj:
        stride = lambda nx, ny, nz: nz+nnz*(ny+nny*nx)
    else:
        stride = lambda ex, ey, ez: ex+nelx*(ez+nelz*ey)
    nx = nnx if proj else nelx
    ny = nny if proj else nely
    nz = nnz if proj else nelz

    # Partition
    npart = npartx*nparty*npartz
    part = []
    startx = 0
    for i in range(npartx):
        if proj:
            nperx = nnx/npartx
            if i < np.mod(nnx, npartx): nperx+=1
        else:
            nperx = nelx/npartx
            if i < np.mod(nelx, npartx): nperx+=1
        ex = np.arange(startx, startx+nperx)
        startx += nperx

        starty = 0
        for j in range(nparty):
            if proj:
                npery = nny/nparty
                if i < np.mod(nny, nparty): npery+=1
            else:
                npery = nely/nparty
                if i < np.mod(nely, nparty): npery+=1
            ey = np.arange(starty, starty+npery)
            starty += npery

            startz = 0
            for k in range(npartz):
                if proj:
                    nperz = nnz/npartz
                    if i < np.mod(nnz, npartz): nperz+=1
                else:
                    nperz = nelz/npartz
                    if i < np.mod(nelz, npartz): nperz+=1
                ez = np.arange(startz, startz+nperz)
                startz += nperz

                cpart = np.zeros((ex.size, ey.size, ez.size), dtype=int,
                                 order='F')
                for ii in range(ex.size):
                    for jj in range(ey.size):
                        for kk in range(ez.size):
                            cpart[ii, jj, kk] = stride(ex[ii], ey[jj], ez[kk])
                part.append( cpart.flatten('F') )
    return ( part )

def partition_param_vec_cluster(param, minstruct, npart, msh, proj=True):

    # Partition into on/off
    part = []
    off = np.array(list(set(np.where(param == 0.0)[0]) - set(minstruct)),
                   dtype=int, order='F')
    on  = np.array(list(set(np.where(param == 1.0)[0]) - set(minstruct)),
                   dtype=int, order='F')
    if proj:
        off_loc = msh.p[:, off]
        on_loc  = msh.p[:, on]
    else:
        off_loc = msh.el_cent[:, off]
        on_loc  = msh.el_cent[:, on]

    # Number in each partition
    npart_off = 0 if len(off) == 0 else int(npart/2)
    npart_on  = npart - npart_off

    # Cluster
    from sklearn import cluster
    if npart_off > 0:
        c_off, l_off, inertia_off = cluster.k_means(off_loc.transpose(),
                                                    npart_off)
    if npart_on  > 0:
        c_on , l_on , inertia_on  = cluster.k_means(on_loc.transpose() ,
                                                    npart_on)

    # Partition
    part = [ ]
    if minstruct.size > 0: part.append(minstruct)
    for i in range(npart_off):
        part.append( off[ np.where(l_off==i)[0] ] )
    for i in range(npart_on):
        part.append( on[ np.where(l_on==i)[0] ] )
    return part, npart_off, npart_on

def partition_param_vec_simple(p, npart):

    part = []
    off = np.where(p == 0.0)[0]
    on  = np.where(p == 1.0)[0]

    # Number in each partition
    if len(off) == 0 and len(on) > 0:
        npart_off = 0; npart_on = npart;
    elif len(on) == 0 and len(off) > 0:
        npart_on = 0 ; npart_off = npart;
    else:
        npart_on  = int(npart/2)
        npart_off = int(npart/2)

    # Partition off/on
    nper_off = off.size / npart_off if npart_off > 0 else 0
    for i in range(npart_off-1):
        #part.append( off[i::nper_off] )
        part.append( off[i*nper_off:(i+1)*nper_off] )
    if npart_off > 0: part.append( off[(i+1)*nper_off:] )

    nper_on = on.size / npart_on if npart_on > 0 else 0
    for i in range(npart_on-1):
        #part.append( on[i::nper_on] )
        part.append( on[i*nper_on:(i+1)*nper_on] )
    if npart_on > 0: part.append( on[(i+1)*nper_on:] )

    return part, npart_off, npart_on
