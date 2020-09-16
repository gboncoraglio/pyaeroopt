import sys, os, pickle, time
import numpy as np
# import mkl
# mkl.set_num_threads(1)

# MPI - if import fails, run in serial
#This envs "conda-python-blender" has mpi4py library installed
start = time.time()
sys.path.append('/home/gbonco/anaconda3/envs/conda-python-blender/lib/python3.5/site-packages')
# sys.path.insert(0,'/home/gbonco/anaconda3/envs/conda-python-blender/lib/python3.5/site-packages')
print('TIME 0a {}'.format(time.time()-start))

import imp

try:

    # (name,fp, pathname, description) = ("mpi4py",None, "/home/gbonco/anaconda3/envs/conda-python-blender/lib/python3.5/site-packages/mpi4py", ('', '', 5))
    # imp.load_module(name, fp, pathname, description)

    from mpi4py import MPI
    parallel = True
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    ncpu = comm.Get_size()
    if rank == 0:
        print("ncpus: ",ncpu)
except ImportError:
    parallel = False
print('RANK {},  TIME 0b {}'.format(rank,time.time()-start))
# pyAeroOpt modules
from pyaeroopt.util.frg_util import write_vmo, write_der
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
print('RANK {},  TIME 1 {}'.format(rank,time.time()-start))
# Extract input using pickle
pkf  = open('tmpFileForPassingClassesToBlenderViaPickle.pkl', 'rb')
pIn = pickle.load(pkf)
pkf.close()
print('RANK {},  TIME 2 {}'.format(rank,time.time()-start))
# Deal pickled structure into appropriate variables
x            = pIn[0]
bObj         = pIn[1]
ptcloud_file = pIn[2]
vmo_file     = pIn[3]
der_file     = pIn[4]
eps          = pIn[5]
xpost_file   = pIn[6]
print('RANK {},  TIME 3 {}'.format(rank,time.time()-start))
# Parallel extension
# if parallel:
#     print('Parallel support = True')
# else:
#     print('Parallel support = False')
if parallel and ncpu > 1:
    vmo_file_base = vmo_file
    der_file_base = der_file

    parallelExt ='.'+str(ncpu)+'parts.part'+str(rank).zfill(len(str(ncpu)))
    ptcloud_file += parallelExt
    vmo_file     += parallelExt
    if der_file is not None: der_file += parallelExt
print('RANK {},  TIME 4 {}'.format(rank,time.time()-start))
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
print('RANK {},  TIME 5 {}'.format(rank,time.time()-start))

print('TIME FOR NODESET/TOP READ = {0:e}'.format(time.time()-t0))
ptcloud = [tuple(pt) for pt in nodes]
print('RANK {},  TIME 6 {}'.format(rank,time.time()-start))

# Count nodes
if parallel:
    nnodes = comm.allreduce(nodes.shape[0], op=MPI.SUM)
    print("NUMBER OF NODES EACH = {0:d}".format(nodes.shape[0]))
else:
    nnodes = nodes.shape[0]
print("NUMBER OF NODES = {0:d}".format(nnodes))
print('RANK {},  TIME 7 {}'.format(rank,time.time()-start))

# Ensure scene is empty
delete_all_objects()
print('RANK {},  TIME 8 {}'.format(rank,time.time()-start))

# Make a Blender object out of all of the modifiers, the nodes of the modifier,
# and then link the modifiers (to enable sequence of modifiers, i.e. Skeleton
# to control Cage to control Lattice to control ptCloud)
bObj.blender_make_deform_object_from_self()
bObj.blender_link_deform_objects()
print('RANK {},  TIME 9 {}'.format(rank,time.time()-start))

# Make Blender object out of ptCloud and set it as the deformee of bObj
ob = create_mesh_object('mesh', ptcloud, [], elems)
bObj.blender_add_deformee(ob)
print('RANK {},  TIME 10 {}'.format(rank,time.time()-start))

# Invoke deformation, extract deformation, and write to file
# Extract displacement from blender object and write to VMO file
t0 = time.time()
disp = bObj.blender_deform(x)
print('RANK {},  TIME 11 {}'.format(rank,time.time()-start))
print('TIME FOR BLENDER LATTICE DEFORM = {0:e}'.format(time.time()-t0))

t0 = time.time()
if not parallel or (parallel and ncpu == 1):
    write_vmo(vmo_file, disp)
else:
    write_vmo(vmo_file, disp, 0 if rank == 0 else None, rank == 0,
              different_size = nnodes)
    comm.Barrier()
    if rank == 0:
        cat_split(vmo_file_base, ncpu, cleanup=True)
    comm.Barrier()
print('RANK {},  TIME 12 {}'.format(rank,time.time()-start))
print('TIME FOR VMO and CAT WRITE = {0:e}'.format(time.time()-t0))

# Make xpost file
if xpost_file is not None:
    if not parallel or (parallel and rank == 0):
        bObj.write_xpost(xpost_file, x)

print('RANK {},  TIME 13 {}'.format(rank,time.time()-start))
# Remove lattice meshObj; only needed for writing xpost files
for mod in bObj.modifiers.get_from_id():
    for def_obj in mod.obj.list_of_deform:
        if def_obj.__class__.__name__ == 'Lattice':
            delete_object( def_obj.blend_objs.mesh_obj )

# Derivatives via finite difference
t0 = time.time()
der_file_gen = lambda it: '{0:s}.shpder{1:d}'.format(der_file_base, it)
der_file_par_gen = lambda it: '{0:s}.shpder{1:d}{2:s}'.format(der_file_base, it,
                                                              parallelExt)
if der_file is not None:
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
print('TIME FOR DER = {0:e}'.format(time.time()-t0))
print('RANK {},  TIME 14 {}'.format(rank,time.time()-start))