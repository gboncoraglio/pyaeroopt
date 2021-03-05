import sys, os, pickle, time
sys.path.insert(0,'/home/users/jbho/codes/blender-2.77a/2.77/python/lib/python3.5/site-packages')# ensures that it uses the numpy for this blender, hardcoded hack
import numpy as np 
#This is compatible with Aerof-2 which uses the abs displacement

#This envs "py35blender" has mpi4py library installed
sys.path.append('/home/users/jbho/codes/anaconda3/envs/py35blender/lib/python3.5/site-packages')
# MPI - if import fails, run in serial
try:
    from mpi4py import MPI
    parallel = True
    der_parallel = True #parallelises derivatives by component during finite difference, rather than parallelizing the computation of each component
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncpu = comm.Get_size()
    if rank == 0:
        print("ncpus: ",ncpu)
except ImportError:
    parallel = False
    der_parallel = False 

# pyAeroOpt modules
from pyaeroopt.util.frg_util import write_vmo,write_vmo_abs, write_der
from pyaeroopt.util.frg_util import cat_split, cat_split_gen
from pyaeroopt.util.frg_util import read_nodeset, read_top, is_nodeset
from pyaeroopt.util.misc     import extents

# Built-in Blender modules
import bpy
from mathutils import Vector, Euler

# George Anderson's custom modules (modified and extended by MJZ)
from pyaeroopt.util.blender.mesh import extract_mesh, create_mesh
from pyaeroopt.util.blender.mesh import create_mesh_object
from pyaeroopt.util.blender.mesh import extract_object_displacement
from pyaeroopt.util.blender.object import delete_object, delete_all_objects
from pyaeroopt.util.blender.modifier import create_lattice, add_modifier

# Extract input using pickle
pkf  = open('tmpFileForPassingClassesToBlenderViaPickle.pkl', 'rb')
pIn = pickle.load(pkf,encoding='latin1')
pkf.close()

# Deal pickled structure into appropriate variables
x            = pIn[0]
bObj         = pIn[1]
ptcloud_file = pIn[2]
vmo_file     = pIn[3]
der_file     = pIn[4]
eps          = pIn[5]
xpost_file   = pIn[6]

# Parallel extension
if parallel:
    print('Parallel support = True')
else:
    print('Parallel support = False')
sys.stdout.flush()
if parallel and ncpu > 1:
    vmo_file_base = vmo_file
    der_file_base = der_file

    parallelExt ='.'+str(ncpu)+'parts.part'+str(rank).zfill(len(str(ncpu)))
    ptcloud_file0 = ptcloud_file
    ptcloud_file += parallelExt
    vmo_file     += parallelExt
    if der_file is not None: der_file += parallelExt

if parallel and rank == 0:
    print(x)
    print(bObj)
    print(ptcloud_file)
    print(vmo_file)
    print(der_file)
    print(eps)
    print(xpost_file)

# Read ptCloud and convert to list of tuples
t0 = time.time()
if is_nodeset(ptcloud_file):
    nodes, elems = read_nodeset(ptcloud_file), []
else:
    nodes, elems = read_top(ptcloud_file)
    nodes = nodes[0][-1]
    elems = [tuple([int(e-1) for e in el]) for el in elems[0][-1]]

if not parallel or rank == 0: print('TIME FOR NODESET/TOP READ = {0:e}'.format(time.time()-t0))
sys.stdout.flush()
ptcloud = [tuple(pt) for pt in nodes]

# Count nodes
if parallel:
    nnodes = comm.allreduce(nodes.shape[0], op=MPI.SUM)
    if rank == 0: print("NUMBER OF NODES EACH = {0:d}".format(nodes.shape[0]))
else:
    nnodes = nodes.shape[0]
if not parallel or rank == 0: print("NUMBER OF NODES = {0:d}".format(nnodes))
sys.stdout.flush()

# Ensure scene is empty
delete_all_objects()

# Make a Blender object out of all of the modifiers, the nodes of the modifier,
# and then link the modifiers (to enable sequence of modifiers, i.e. Skeleton
# to control Cage to control Lattice to control ptCloud)
bObj.blender_make_deform_object_from_self()
bObj.blender_link_deform_objects()

# Make Blender object out of ptCloud and set it as the deformee of bObj
ob = create_mesh_object('mesh', ptcloud, [], elems)
bObj.blender_add_deformee(ob)

# Invoke deformation, extract deformation, and write to file
# Extract displacement from blender object and write to VMO file
t0 = time.time()
disp = bObj.blender_deform(x)
if not parallel or rank == 0: print('TIME FOR BLENDER LATTICE DEFORM = {0:e}'.format(time.time()-t0))
sys.stdout.flush()

t0 = time.time()
if not parallel or (parallel and ncpu == 1):
    # write_vmo(vmo_file, disp)
    write_vmo_abs(vmo_file, disp,nodes)
else:
    # write_vmo(vmo_file, disp, 0 if rank == 0 else None, rank == 0,
    #           different_size = nnodes)
    write_vmo_abs(vmo_file, disp,nodes, 0 if rank == 0 else None, rank == 0,
              different_size = nnodes)
    comm.Barrier()
    if rank == 0:
        cat_split(vmo_file_base, ncpu, cleanup=True)
    comm.Barrier()
if not parallel or rank == 0: print('TIME FOR VMO and CAT WRITE = {0:e}'.format(time.time()-t0))
sys.stdout.flush()

# Make xpost file
if xpost_file is not None:
    if not parallel or (parallel and rank == 0):
        bObj.write_xpost(xpost_file, x)

# Remove lattice meshObj; only needed for writing xpost files
for mod in bObj.modifiers.get_from_id():
    for def_obj in mod.obj.list_of_deform:
        if def_obj.__class__.__name__ == 'Lattice':
            delete_object( def_obj.blend_objs.mesh_obj )

# Derivatives via finite difference
der_file_gen = lambda it: '{0:s}.shpder{1:d}'.format(der_file_base, it)
der_file_par_gen = lambda it: '{0:s}.shpder{1:d}{2:s}'.format(der_file_base, it,
                                                              parallelExt)
if der_file is not None:
  tderiv0 = time.time()
  if not der_parallel: 
    for i, xi in enumerate(x):
        ei = np.zeros(x.size, dtype=float)
        ei[i] = eps

        # Deformation at x + h and x - h
        dp = bObj.blender_deform(x+ei)
        dm = bObj.blender_deform(x-ei)

        # Write derivative
        if not parallel or (parallel and ncpu == 1):
            write_der(der_file, (0.5/eps)*(dp-dm), [i])
        else:
            # write_vmo used since writing individual file for each derivative
            # which will be 'cat'(ed) later. write_der assumes single file
            # written directly for all derivatives
            write_vmo(der_file_par_gen(i), (0.5/eps)*(dp-dm),
                      i if rank == 0 else None,
                      rank == 0 and i == 0,
                      different_size = nnodes)
            comm.Barrier()
            if rank == 0:
                cat_split(der_file_gen(i), ncpu, cleanup=True)
            comm.Barrier()
    if parallel and ncpu > 1:
        if rank == 0:
            cat_split_gen(der_file_base, der_file_gen, len(x), cleanup=True)
  else:
    #read entire nodeset into each process
    if is_nodeset(ptcloud_file0):
        nodes, elems = read_nodeset(ptcloud_file0), []
    else:
        nodes, elems = read_top(ptcloud_file0)
        elems = [tuple([int(e-1) for e in el]) for el in elems[0][-1]]

    ptcloud = [tuple(pt) for pt in nodes]

    # clearscene
    delete_all_objects()

    # Make a Blender object out of all of the modifiers, the nodes of the modifier,
    # and then link the modifiers (to enable sequence of modifiers, i.e. Skeleton
    # to control Cage to control Lattice to control ptCloud)
    bObj.blender_make_deform_object_from_self()
    bObj.blender_link_deform_objects()

    # Make Blender object out of ptCloud and set it as the deformee of bObj
    ob = create_mesh_object('mesh', ptcloud, [], elems)
    bObj.blender_add_deformee(ob)

    
    n = len(x)
    n_percpu = int(n / ncpu)
    nextra = n - n_percpu * ncpu
    if rank < nextra:
      n_percpu = n_percpu + 1 
      starting_i = rank * n_percpu
    else:
      starting_i = rank * n_percpu + nextra
    for i in range(starting_i, starting_i + n_percpu):
      t0 = time.time()
      xi = x[i]
      ei = np.zeros(x.size, dtype=float)
      ei[i] = eps

      # Deformation at x + h and x - h
      dp = bObj.blender_deform(x+ei)
      dm = bObj.blender_deform(x-ei)
      # write_vmo used since writing individual file for each component of the derivative
      # which will be 'cat'(ed) later. write_der assumes single file
      write_vmo(der_file_gen(i), (0.5/eps)*(dp-dm), i, i == 0)
    comm.Barrier()
    if rank == 0:
      t0 = time.time()
      cat_split_gen(der_file_base, der_file_gen, n, cleanup=True)
  if not parallel or rank == 0: print('TIME FOR DER TOTAL = {0:e}'.format(time.time()-tderiv0))
  sys.stdout.flush()
