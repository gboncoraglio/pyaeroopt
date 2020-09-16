import sys, os, time, copy
import numpy as np

sys.path.append('/Users/mzahr/Documents/codes/pymortestbed/dev/')
sys.path.append('/Users/mzahr/Documents/codes/pymortestbed/dev/'+
                'pymortestbed/app/struct_dynam_3d/src/_cython')
from pymortestbed.app.struct_dynam_3d.util import read_aeros
from pymortestbed.app.struct_dynam_3d.util import construct_hex_mesh, write_vtk
from pymortestbed.app.struct_dynam_3d.util import remove_inside_circle
from pymortestbed.app.struct_dynam_3d.util import remove_inside_shape

### Mesh, Physics, Discretization construction functions ###
def mesh_properties(nx, ny, nz, shape=None, param=None):
    p, t_ptr, t, etype = construct_hex_mesh(nx, ny, nz)
    if shape is not None and param is not None:
        p, t_ptr, t, etype = remove_inside_shape(p, t_ptr, t, etype, shape,
                                                 param)
    nsd = p.shape[0]
    nv  = p.shape[1]
    nel = t_ptr.size-1
    return nsd, nel, nv, p, t_ptr, t, etype

def dirichlet_bcs(p, where, val):
    dbc_loc = np.zeros(p.shape, dtype=bool, order='F')
    for j in range(p.shape[1]):
        if where(p[:,j]):
            dbc_loc[:, j] = True
    ndbc = np.sum(dbc_loc)
    dbc_val = val*np.ones(ndbc, dtype=float, order='F')
    return ndbc, dbc_loc, dbc_val

def force_bcs(p, ndof, where, comp, val, dbc_loc):
    fbc_loc = np.zeros(p.shape, dtype=bool, order='F')
    for j in range(p.shape[1]):
        if where(p[:, j]):
            fbc_loc[comp, j] = True
    nfbc = np.sum(fbc_loc)
    fext0 = np.zeros(p.shape, dtype=float, order='F')
    fext0[fbc_loc] = val/np.float(nfbc)
    fext = extract_dof_from_nodes(ndof, dbc_loc, fext0)
    return nfbc, fbc_loc, fext0, fext

def physics_properties(nel, mat0):
    claw = np.zeros(nel, dtype=int) # StVK everywhere
    nmat = mat0.size
    mat_ptr = np.arange(0, nmat*nel+1, nmat)
    mat = np.tile( mat0, (1, nel) ).flatten()
    nmat_per_el = nmat*np.ones(nel, dtype=int, order='F')
    return claw, mat_ptr, mat, nmat_per_el

def discretization_properties(eshape=1, nquad=6):
    return eshape, nquad

### Unstructured setup ###
def setup_uns(fname, which, forc, mat00, dbc, fbc=None):

    # Read AEROS/XPOST
    p, t_ptr, t, etype, dbc_loc, dbc_val, fext0 = read_aeros(fname)
    if which == 'lax0':
        ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext = \
                                                    lax_setup(p, forc, dbc, fbc)
    # Set values
    nsd = 3
    nel = t_ptr.size-1
    nv  = p.shape[1]

    # Physics/Discretization
    mat0 = np.array(mat00, dtype=float, order='F')
    claw, mat_ptr, mat, nmat_per_el = physics_properties(nel, mat0)
    eshape, nquad      = discretization_properties(0, 4) # TET

    return ( nsd, nel, nv, p, t_ptr, t, etype, ndof, ndbc, dbc_loc, dbc_val,
             nfbc, fbc_loc, fext0, fext, mat0, claw, mat_ptr, mat, nmat_per_el,
             eshape, nquad )

def lax_setup(p, forc, dbc_cyl, fbc):

    # Dirichlet boundary conditions (fixed at x = xlim[0])
    dbc_center = dbc_cyl[0]
    dbc_radius = dbc_cyl[1]
    dbc_lim    = dbc_cyl[2]
    dbcfcn = ( lambda v: v[2] > dbc_lim and
              (v[0]-dbc_center[0])**2+(v[1]-dbc_center[1])**2 < dbc_radius**2 )
    ndbc, dbc_loc, dbc_val = dirichlet_bcs(p, dbcfcn, 0.0)
    ndof = p.size - ndbc


    # Force boundary conditions (forc at x = xlim[1], y = 0.5*(zlim[0]+zlim[1]))
    nfbc, fbc_loc, fext0, fext = force_bcs(p, ndof, lambda v: (v[0] < fbc[0]),
                                           0, forc, dbc_loc)
    return ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext

# Setup (general)
def setup(which, nelx, nely, nelz, xlim, ylim, zlim, forc, mat00, shape=None,
          param=None):

    # Location of nodes
    nx = np.linspace(xlim[0], xlim[1], nelx+1)
    ny = np.linspace(ylim[0], ylim[1], nely+1)
    nz = np.linspace(zlim[0], zlim[1], nelz+1)

    # Create mesh quantities
    nsd, nel, nv, p, t_ptr, t, etype = mesh_properties(nx, ny, nz, shape, param)

    if which == 'cantilever':
        ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext = \
                   cantilever_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz)
    elif which == 'trestle':
        ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext = \
                   trestle_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz)
    elif which == 'michell':
        ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext = \
               michell_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz,
                             shape, param)
    elif which == 'lacross0':
        ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext = \
               lacross0_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz,
                              shape, param)
    # Physics/Discretization
    mat0 = np.array(mat00, dtype=float, order='F')
    claw, mat_ptr, mat, nmat_per_el = physics_properties(nel, mat0)
    eshape, nquad      = discretization_properties()

    return ( nsd, nel, nv, p, t_ptr, t, etype, ndof, ndbc, dbc_loc, dbc_val,
             nfbc, fbc_loc, fext0, fext, mat0, claw, mat_ptr, mat, nmat_per_el,
            eshape, nquad )

### Setup specific (cantilever, michell, trestle, lacrossblock0) ###
def cantilever_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz):
    """ Cantilever setup: Requires EVEN number of elements in z-direction! """

    # Dirichlet boundary conditions (fixed at x = xlim[0])
    ndbc, dbc_loc, dbc_val = dirichlet_bcs(p, lambda v: v[0] == xlim[0], 0.0)
    ndof = p.size - ndbc

    # Force boundary conditions (forc at x = xlim[1], y = 0.5*(zlim[0]+zlim[1]))
    nfbc, fbc_loc, fext0, fext = force_bcs(p, ndof,
                  lambda v: (v[0] == xlim[1] and v[2] == 0.5*(zlim[0]+zlim[1])),
                  2, forc, dbc_loc)
    return ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext

def michell_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz, c, r):
    """ Michell setup: Requires EVEN number of elements in z-direction! """

    dx = xlim[1] - xlim[0]; dy = ylim[1] - ylim[0]; dz = zlim[1] - zlim[0];
    hx = dx / float(nelx) ; hy = dy / float(nely) ; hz = dz / float(nelz) ;
    h = min([hx, hy, hz])

    # Dirichlet boundary conditions (fixed at x = xlim[0])
    ndbc, dbc_loc, dbc_val = dirichlet_bcs(p,
                                   lambda v: np.sum((v-c)**2) < (r+h/2)**2, 0.0)
    ndof = p.size - ndbc

    # Force boundary conditions (forc at x = xlim[1], y = 0.5*(zlim[0]+zlim[1]))
    nfbc, fbc_loc, fext0, fext = force_bcs(p, ndof,
                  lambda v: (v[0] == xlim[1] and v[2] == 0.5*(zlim[0]+zlim[1])),
                  2, forc, dbc_loc)
    return ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext

def trestle_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz):
    """ Trestle setup """

    dx = xlim[1] - xlim[0]; dy = ylim[1] - ylim[0]; dz = zlim[1] - zlim[0];
    hx = dx / float(nelx) ; hy = dy / float(nely) ; hz = dz / float(nelz) ;

    # Dirichlet boundary conditions (fixed at x = xlim[0])
    xarea1 = xarea2 = xlim
    yarea1 = yarea2 = ylim
    #xarea1 = [xlim[0] - 2*hx, xlim[0] + 2*hx]
    #xarea2 = [xlim[1] - 2*hx, xlim[1] + 2*hx]
    #yarea1 = [ylim[0] - 2*hy, ylim[0] + 2*hy]
    #yarea2 = [ylim[1] - 2*hy, ylim[1] + 2*hy]

    dbc_func = lambda v: ( (( v[0] < xarea1[1] and v[0] > xarea1[0] and
                              v[1] < yarea1[1] and v[1] > yarea1[0] ) or
                            ( v[0] < xarea2[1] and v[0] > xarea2[0] and
                              v[1] < yarea1[1] and v[1] > yarea1[0] ) or
                            ( v[0] < xarea1[1] and v[0] > xarea1[0] and
                              v[1] < yarea2[1] and v[1] > yarea2[0] ) or
                            ( v[0] < xarea2[1] and v[0] > xarea2[0] and
                              v[1] < yarea2[1] and v[1] > yarea2[0] )) and
                              v[2] == zlim[0] )
    ndbc, dbc_loc, dbc_val = dirichlet_bcs(p, dbc_func, 0.0)
    ndof = p.size - ndbc

    # Force boundary conditions (forc at x = xlim[1], y = 0.5*(zlim[0]+zlim[1]))
    xload = [xlim[0] + 0.5*dx - 2*hx, xlim[0] + 0.5*dx + 2*hx]
    yload = [ylim[0] + 0.5*dy - 2*hy, ylim[0] + 0.5*dy + 2*hy]
    fbc_func = lambda v: ( (v[0] < xload[1] and v[0] > xload[0]) and
                           (v[1] < yload[1] and v[1] > yload[0]) and
                            v[2] == zlim[1] )
    nfbc0, fbc_loc0, fext00, fext0 = force_bcs(p, ndof, fbc_func, 2, forc,
                                               dbc_loc)

    fbc_func = lambda v: ( v[0] == xlim[0] and
                           v[2] > zlim[0] + 0.75*(zlim[1]-zlim[0]) )
    nfbc1, fbc_loc1, fext01, fext1 = force_bcs(p, ndof, fbc_func, 1, forc,
                                               dbc_loc)

    fbc_func = lambda v: ( v[0] == xlim[1] and
                           v[2] > zlim[0] + 0.75*(zlim[1]-zlim[0]) )
    nfbc2, fbc_loc2, fext02, fext2 = force_bcs(p, ndof, fbc_func, 1, -forc,
                                               dbc_loc)

    fbc_func = lambda v: ( v[1] == ylim[0] and
                           v[2] > zlim[0] + 0.75*(zlim[1]-zlim[0]) )
    nfbc3, fbc_loc3, fext03, fext3 = force_bcs(p, ndof, fbc_func, 2, -forc,
                                               dbc_loc)

    fbc_func = lambda v: ( v[1] == ylim[1] and
                           v[2] > zlim[0] + 0.75*(zlim[1]-zlim[0]) )
    nfbc4, fbc_loc4, fext04, fext4 = force_bcs(p, ndof, fbc_func, 2, forc,
                                               dbc_loc)

    nfbc = nfbc0 + nfbc1 + nfbc2 + nfbc3 + nfbc4
    fbc_loc = fbc_loc0 + fbc_loc1 + fbc_loc2 + fbc_loc3 + fbc_loc4
    fext = fext0 + fext1 + fext2 + fext3 + fext4
    fext0 = fext00 + fext01 + fext02 + fext03 + fext04
    return ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext

def lacrossblock0_setup(p, xlim, ylim, zlim, forc, nelx, nely, nelz, c, r):
    """ Lacross block setup """

    dx = xlim[1] - xlim[0]; dy = ylim[1] - ylim[0]; dz = zlim[1] - zlim[0];
    hx = dx / float(nelx) ; hy = dy / float(nely) ; hz = dz / float(nelz) ;

    # Dirichlet boundary conditions
    dbcfcn = lambda v: (v[1] < ylim[0]+0.1*dy) and ((v[0]-(xlim[0]+0.5*dx))**2
                        + (v[2]-(zlim[0]+0.5*dz))**2 < (0.25*min(dx, dz))**2)
    ndbc, dbc_loc, dbc_val = dirichlet_bcs(p, dbcfcn, 0.0)
    ndof = p.size - ndbc

    # Force boundary conditions
    nfbc, fbc_loc, fext0, fext = force_bcs(p, ndof,
                  lambda v: (v[0] == xlim[0] and v[1] > ylim[1]-0.1*dy and
                                                 v[2] > zlim[1]-0.2*dz),
                  2, forc, dbc_loc)
    return ndof, ndbc, dbc_loc, dbc_val, nfbc, fbc_loc, fext0, fext


### Helper functions ###
def extract_dof_from_nodes(ndof, dbc_loc, Fnodes):
    nsd, nv = dbc_loc.shape
    Fdof = np.zeros(ndof, dtype=float, order='F')
    cnt = -1
    for j in range(nv):
        for i in range(nsd):
            if not dbc_loc[i, j]:
                cnt +=1
                Fdof[cnt] = Fnodes[i, j]
    return ( Fdof )

def construct_full_vec(dbc_loc, dof, dbc):
    nsd, nv = dbc_loc.shape
    nodal = np.zeros((nsd, nv), dtype=float, order='F')
    dbccnt = -1
    dofcnt = -1
    for j in range(nv):
        for i in range(nsd):
            if not dbc_loc[i, j]:
                dofcnt += 1
                nodal[i, j] = dof[dofcnt]
            else:
                dbccnt += 1
                nodal[i, j] = dbc[dbccnt]
    return ( nodal )
