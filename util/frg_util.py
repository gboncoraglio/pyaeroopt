import os, copy, subprocess, sys
import numpy as np

from pyaeroopt.util.hpc_util import execute_code
from pyaeroopt.util.misc import count_lines, split_line_robust, is_numeric

def is_nodeset(fname):
    f = open(fname)
    line = f.readline().rstrip('\n')
    f.close()
    return line == 'FENODES'
    # return line == 'NODE'

def top_stats(fname, only_nodes=False):
    pass

def nodeset_stats(fname):
    return ['nodes', count_lines(fname)-2]

def xpost_stats(fname):

    nline = count_lines(fname)
    out = {'desc':None, 'nnode':None, 'ndof':None, 'nstep':None}

    f = open(fname,'r')
    for k, line in enumerate(f):
        if k == 0:
            out['desc']  = line.rstrip('\n')
        if k == 1:
            out['nnode'] = int(line.rstrip('\n'))
        if k == 3:
            out['ndof']  = len( [x for x in line.split(' ') if x] )
            break
    f.close()

    out['nstep'] = (nline-2)/(out['nnode']+1)
    return out

def read_top(fname, only_nodes=False):

    # Read all lines from file
    f = open(fname,'r')
    lines = f.readlines()
    f.close()

    nodes, elems = [], []
    for k, line in enumerate(lines):
        # Determine section currently reading from; create appropriate data
        # structures to store statistics and values of mesh.
        if 'Nodes' in line.lstrip(' ').lstrip('\t')[:5]:
            inNodes = True
            inElems = False

            ne = 0

            # Seek ahead to find next instance of Nodes or Elems to determine
            # size of node set
            for nn, lline in enumerate(lines[k+1:]):
                if ('Nodes' in lline.lstrip(' ').lstrip('\t')[:5] or
                    'Elements' in lline.lstrip(' ').lstrip('\t')[:8]):
                    break
            if lline is lines[-1]: nn+=1

            lineStrip = line.lstrip(' ').lstrip('\t')
            lineStrip = lineStrip.rstrip(' ').rstrip('\t').rstrip('\n')
            nodes.append([x for x in lineStrip.split(' ') if x])
            nodes[-1].append( np.zeros(nn,dtype=int) )
            nodes[-1].append( np.zeros((nn,3),dtype=float) )
            continue

        elif 'Elements' in line.lstrip(' ').lstrip('\t')[:8]:
            inNodes = False
            if only_nodes: continue

            inElems = True

            ne = 0

            # Use next line to determine number of nodes per element
            # ASSUMES ALL ELEMENTS IN ELEMENT SET ARE OF SAME TYPE
            nen = len( split_line_robust(lines[k+1]) ) - 2

            # Seek ahead to find next instance of Nodes or Elems to determine
            # size of element set
            for nn, lline in enumerate(lines[k+1:]):
                if ('Nodes' in lline.lstrip(' ').lstrip('\t')[:5] or
                    'Elements' in lline.lstrip(' ').lstrip('\t')[:8]):
                    break
            if lline is lines[-1]: nn+=1

            lineStrip = line.lstrip(' ').lstrip('\t')
            lineStrip = lineStrip.rstrip(' ').rstrip('\n').rstrip('\t')
            sep = '\t' if '\t' in lineStrip else ' '
            elems.append([x for x in lineStrip.split(sep) if x and
                                            x != 'using' and not is_numeric(x)])
            elems[-1].append( np.zeros(nn,dtype=int) )
            elems[-1].append( np.zeros(nn,dtype=int) )
            elems[-1].append( np.zeros((nn,nen),dtype=int) )

            continue

        # Split/strip line (lines might contain ' ' and \t as deliminaters)
        lineSplit = split_line_robust(line)

        # Add node number and coordinates to current node
        if inNodes:
            nodes[-1][-2][ne] = int(lineSplit[0])
            nodes[-1][-1][ne,:] = np.array([float(x) for x in lineSplit[1:]]
                                                               ).reshape(1,-1)

        # Add element number, type, and nodes to current element
        if inElems:
            elems[-1][-3][ne] = int(lineSplit[0])
            elems[-1][-2][ne] = int(lineSplit[1])
            elems[-1][-1][ne,:] = np.array([int(x) for x in lineSplit[2:]]
                                                               ).reshape(1,-1)
        ne+=1
    return ( nodes, elems )

def read_nodeset(fname):
    return np.loadtxt(fname, dtype='float', skiprows=2, comments='END')[:, 1:]

def read_vmo(fname):
    return np.loadtxt(fname, dtype='float', skiprows=3)

def read_xpost(fname):

    # Read all lines from file
    f = open(fname,'r')
    lines = f.readlines()
    f.close()

    nnode = int( lines[1].rstrip('\n') )
    ndof  = len( [x for x in lines[3].rstrip('\n').split(' ') if x] )
    nstep = ( len(lines) - 2 ) / ( nnode + 1)

    # Extract "time" information from lines
    time = np.array([float(x.rstrip('\n')) for x in lines[2:-1:nnode+1]])

    # Extract value information from lines
    val = np.zeros((nnode, ndof, int(nstep) ))
    for t, tim in enumerate(time):
        for k, line in enumerate(lines[3+t*(nnode+1):3+(t+1)*nnode+t]):
            val[k,:,t] = np.array( [float(x.rstrip('\n')) for
                                                    x in line.split(' ') if x] )
    return ( time, val )

def read_vmo(fname, coords=[0, 1, 2]):
    return np.loadtxt(fname, dtype='float', skiprows=3)[:, coords]

def read_der(fname, coords=[0, 1, 2]):
    return read_xpost(fname)[1][:, coords, :]

def write_top(fname, nodes, elems):

    f = open(fname,'w')

    # Write nodes
    for nodeset in nodes:
        f.write(nodeset[0]+' '+nodeset[1]+'\n')
        for k, nnum in enumerate(nodeset[2]):
            # f.write(str(nnum)+'  '+str(nodeset[3][k,0])+' '+
            #                  str(nodeset[3][k,1])+' '+str(nodeset[3][k,2])+'\n')
            f.write("{0:d} {1:.16e} {2:.16e} {3:.16e} \n".format(nnum, nodeset[3][k,0], nodeset[3][k,1],nodeset[3][k,2])) #for more accuracy

    # Write elements
    for elemset in elems:
        f.write(elemset[0]+' '+elemset[1]+' using '+elemset[2]+'\n')
        for k, enum in enumerate(elemset[3]):
            f.write(str(enum)+'    '+str(elemset[4][k]))
            # print(elemset[5][k,:])
            for e in elemset[5][k,:]:
                f.write('    '+str(e))
            f.write('\n')

    f.close()

'''
********************************************
#change in position of arguments
#change in first and last lines of the file
********************************************
'''
def write_nodeset(fname, nodes , node_nums=None):
    f = open(fname,'w')
    # f.write('FENODES\n')
    f.write('NODE\n')
    f.write('*\n')
    for k, node in enumerate(nodes):
        if node_nums is None:
            f.write(str(k+1)+' ')
        else:
            f.write(str(node_nums[k])+' ')
        f.write(str(node[0])+' '+str(node[1])+' '+str(node[2])+'\n')
    # f.write('END')
    f.close()


def write_nodeset_1(fname, nodes , node_nums=None):
    f = open(fname,'w')
    f.write('FENODES\n')
    # f.write('NODE\n')
    f.write('*\n')
    for k, node in enumerate(nodes):
        if node_nums is None:
            f.write(str(k+1)+' ')
        else:
            f.write(str(node_nums[k])+' ')
        f.write(str(node[0])+' '+str(node[1])+' '+str(node[2])+'\n')
    f.write('END')
    f.close()

def write_xpost(fname, name, tags, vals, header=None):

    # Shape of data to write
    if len(vals.shape) < 3: vals = vals[:, :, None]
    nnode, ndof, nstep = vals.shape

    # nstep nneds to match tags and vals
    if len(tags) != nstep:
        raise ValueError('len(tags) and vals.shape[2] must be equal')

    # Open file and write header
    f = open(fname,'w')
    f.write('Scalar ') if vals.shape[1]==1 else f.write('Vector ')
    f.write(header+' under load for '+name+'\n')
    f.write(str(nnode)+'\n')
    for t, time in enumerate(tags):
        f.write(str(time)+'\n')
        for v, val in enumerate(vals[:,:,t]):
            f.write('  '.join([str(vv) for vv in val])+'\n')
    f.close()

def combine_xpost(fname, files):

    # Statistics from xpost file
    stats = xpost_stats(files[0])
    nnode = stats['nnode']
    ndof  = stats['ndof']
    desc  = stats['desc']
    nfile = len(files)

    # Read/concatenate from each file
    time = np.zeros(nfile, dtype=float)
    val  = np.zeros((nnode, ndof, nfile), dtype=float)
    for k, file in enumerate(files):
        _, tmp = read_xpost(file)
        time[k]    = k
        val[:, :, k] = tmp[:, :, -1]

    # Write xpost file
    s = desc.split(' ')
    write_xpost(fname, s[-1], time, val, s[1])

def top2nodeset(top_name, nodeset_name):
    nodes, _ = read_top(top_name, True)
    write_nodeset( nodeset_name,nodes[0][-1], nodes[0][-2])
    # write_nodeset(nodes[0][-1], nodeset_name, nodes[0][-2])

def top2nodeset_1(top_name, nodeset_name):
    nodes, _ = read_top(top_name, True)
    write_nodeset_1( nodeset_name,nodes[0][-1], nodes[0][-2])

def flip_coord(fname, fname_flipped, swap=[0, 2, 1], ftype='top'):
    pass

def write_vmo(fname, dat, step_num=0, with_header=True, different_size=None):

    f = open(fname,'w')
    if with_header:
        f.write('Vector MODE under Attributes for nset\n')
        if different_size is not None:
            f.write(str(different_size)+'\n')
        else:
            f.write(str(dat.shape[0])+'\n')
    if step_num is not None:
        f.write('  {0:d}\n'.format(step_num))
    for d in dat:
        # f.write(str(d[0])+'  '+str(d[1])+'  '+str(d[2])+'\n')
        f.write("{0:.16e} {1:.16e} {2:.16e}\n".format(d[0], d[1], d[2]))
    f.close()

def write_vmo_abs(fname, dat,nodes_abs, step_num=0, with_header=True, different_size=None):

    f = open(fname,'w')
    if with_header:
        f.write('Vector MODE under Attributes for nset\n')
        if different_size is not None:
            f.write(str(different_size)+'\n')
        else:
            f.write(str(dat.shape[0])+'\n')
    if step_num is not None:
        f.write('  {0:d}\n'.format(step_num))
    # for d in dat:
    for i in range(dat.shape[0]):
        f.write(str(dat[i][0]+nodes_abs[i][0])+'  '+str(dat[i][1]+nodes_abs[i][1])+'  '+str(dat[i][2]+nodes_abs[i][2])+'\n')
    f.close()


def write_der(fname, dat, vars=None):

    # Write all derivatives to file
    if vars is None and len(dat.shape) == 3:
        vars = np.arange(dat.shape[2])

    # Make sure var iterable
    if not hasattr(vars,'__getitem__'):
        vars = [vars]

    # If writing variable 0, open file and write header
    if vars[0] == 0:
        f = open(fname,'w')
        f.write('Vector MODE under Attributes for FluidNodes\n')
        f.write(str(dat.shape[0])+'\n')
    else:
        f = open(fname,'a')

    # Write derivatives
    for var in vars:
        f.write('  '+str(var)+'\n')
        for n in np.arange(dat.shape[0]):
            if len(dat.shape) == 2:
                f.write(str(dat[n,0])+'  '+str(dat[n,1])+'  '+
                                                            str(dat[n,2])+'\n')
            elif len(dat.shape) == 3:
                f.write(str(dat[n,0,var])+'  '+str(dat[n,1,var])+'  '+
                                                        str(dat[n,2,var])+'\n')
    f.close()


def read_multsoln(fname):
    multsoln = []
    f = open(fname, 'r')
    for line in f:
        multsoln.append(line.rstrip('\n'))
    f.close()
    return multsoln

def write_multsoln(fname, multsoln):
    f = open(fname, 'w')
    f.write('{0:d}'.format(len(multsoln)))
    for i, imultsoln in enumerate(multsoln):
        f.write('{0:s}'.format(imultsoln))
    f.close()

def flip_coord_top(fname_in, fname_out, swap=[0, 2, 1],scale = 1.0):

    f_in  = open(fname_in, 'r')
    f_out = open(fname_out, 'w')
    for line in f_in:
        if 'Elements' in line:
            underNodes = False
        elif 'Nodes' in line:
            underNodes = True
            f_out.write(line)
            continue

        if underNodes:
            if swap == [0,1,2] or swap == [2,0,1] or swap ==  [1,2,0]: sign= 1.0
            else:                                                      sign=-1.0

            coords = [x for x in line.rstrip('\n').split() if x]
            if scale != 1.0:
              f_out.write(coords[0]+' '+str(scale*float(coords[swap[0]+1]))+
                                  ' '+str(scale*float(coords[swap[1]+1]))+
                                  ' '+str(sign*scale*float(coords[swap[2]+1]))+'\n')
            else:
              f_out.write(coords[0]+' '+coords[swap[0]+1]+
                                  ' '+coords[swap[1]+1]+
                                  ' '+str(sign*float(coords[swap[2]+1]))+'\n')
        else:
            f_out.write(line)
    f_in.close()
    f_out.close()

def displace_top_with_vmo(fname, top, vmo, scale = 1.0, scale_disp = 1.0):
    nodes, elems   = read_top(top)
    disp           = read_vmo(vmo)
    nodes[-1][-1] += scale_disp * disp
    nodes[-1][-1] *= scale
    write_top(fname, nodes, elems)

def join_tops(fname, top1, top2):#for joining two top files each with only one nodeset
#keeps nodeset and elemsets separate
    nodes1, elems1   = read_top(top1)#assumes it only has one nodeset
    nodeoffset = np.max(nodes1[-1][2])
    elemoffset = np.max(elems1[-1][3])
    nodes2, elems2   = read_top(top2)#assumes it only has one nodeset
    nodes2[-1][2] += nodeoffset
    elems2[-1][3] += elemoffset
    elems2[-1][-1] += nodeoffset #node numbering in the elements
    nodes1[-1][2] = np.hstack((nodes1[-1][2],nodes2[-1][2])) #join nodesets together 
    nodes1[-1][3] = np.vstack((nodes1[-1][3],nodes2[-1][3])) #join nodeset together
    elems1[-1][1] += "_1"
    elems2[-1][1] += "_2"
    elems2[-1][2] = elems1[-1][2]
    write_top(fname, nodes1, [elems1[-1], elems2[-1]])
def join_derivs(fname, der1, der2):
    t1, val1 = read_xpost(der1)
    nnodes1, dim1, nder1 = val1.shape
    t2, val2 = read_xpost(der2)
    nnodes2, dim2, nder2 = val2.shape
    if dim1 != dim2:
      print("Dimensions of der1 and der2 do not match in join join_derivs")
      exit()

    dim = dim1
    der = np.zeros((nnodes1 + nnodes2, dim, nder1 + nder2))
    for i in range(nder1):
      der[:nnodes1,:,i] = val1[:,:,i]
    for i in range(nder2):
      der[nnodes1:,:,nder1+i] = val2[:,:,i]
    write_xpost(fname, "Attributes", range(nder1 + nder2), der, header = 'MODE')
def concatenate_xposts(fname, xpost1, xpost2, concatenateTimes = True):#for concatenating multiple xposts for the same nodeset
#keeps nodeset and elemsets separate
    t1, val1 = read_xpost(xpost1)
    nnodes1, dim1, n1 = val1.shape
    t2, val2 = read_xpost(xpost2)
    nnodes2, dim2, n2 = val2.shape
    if nnodes1 != nnodes2 or dim1 != dim2:
      print("Cannot concatenate two xposts with different number of nodes and/or dimensions")
      exit()
    if concatenateTimes:
      t = np.concatenate((t1, t2 + t1[-1] + 1))
    else:
      t = np.concatenate((t1, t2))
    write_xpost(fname, "Attributes", t, np.concatenate((val1,val2),axis=2), header = 'MODE')
def rewrite_top(fname, top, elemmask):#rewrites nodeset with mask
    nodes, elems = read_top(top)
    mask = np.genfromtxt(elemmask)
    elemsnew = [copy.copy(elems[0]), copy.copy(elems[0])]
    for i in range(2):
      elemsnew[i][1] += str(i)
      for j in [3,4,5]:
        elemsnew[i][j] = elems[0][j][mask==i]
    write_top(fname, nodes, elemsnew)
#from scipy.spatial.transform import Rotation as R
def rotate_top_multiple_axes(fname, fname_derivs, top, axes, origin, thetas, scale = 1.0, top2 = None):
    naxis = len(axes)
    nodes, elems   = read_top(top)
    nnodes, dim = nodes[-1][-1].shape
    nextra = 0
    if top2 is not None:
      nodes2, elems2 = read_top(top2)
      nextra = nodes2[-1][-1].shape[0]
    der = np.zeros((nnodes + nextra, dim,naxis))
    for i in range(naxis):
      axis = axes[i] / np.sqrt(np.sum(np.power(axes[i],2)))
      #rotation
      K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
      t = np.radians(thetas[i])
      R = np.eye(3) + np.sin(t)*K + (1 - np.cos(t))*np.matmul(K, K)
      nodes[-1][-1] = np.matmul(nodes[-1][-1],R.transpose()) -  R.dot(origin) + origin
      #for k in range(i):#rotate derivatives as well
      #  der[:,:,k] = np.matmul(der[:,:,k],R.transpose()) -  R.dot(origin) + origin
      for j in range(nnodes):
        der[j,:,i] = np.cross(axis,nodes[-1][-1][j] - origin) / scale 
      write_top(fname, nodes, elems)
    write_xpost(fname_derivs, "Attributes", range(naxis), der, header='MODE')
def rotate_top(fname, top, axis, origin, theta):
    nodes, elems   = read_top(top)
    axis = axis / np.sqrt(np.sum(np.power(axis,2)))
    #rotation
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    t = np.radians(theta)
    R = np.eye(3) + np.sin(t)*K + (1 - np.cos(t))*np.matmul(K, K)
    nodes[-1][-1] = np.matmul(nodes[-1][-1],R.transpose()) -  R.dot(origin) + origin
    write_top(fname, nodes, elems)
def dihedral_deform(fname, top, s):
    nodes, elems   = read_top(top)
    nnodes = nodes[-1][-1].shape[0]
    for i in range(nnodes):
      absy = np.fabs(nodes[-1][-1][i][1])
      if absy > 1.2:
        nodes[-1][-1][i][2] += 0.001 * s * (absy - 1.2) ** 2
    write_top(fname, nodes, elems)
def rotate_derivs(fname, top, axis, origin, scale = 1.0, top2 = None):
    #top2 is a second top file with which top is joined which means the xpost needs to be padded with zeros
    nodes, elems = read_top(top)
    axis = axis / np.sqrt(np.sum(np.power(axis,2)))
    nnodes, dim = nodes[-1][-1].shape
    nextra = 0
    if top2 is not None:
      nodes2, elems2 = read_top(top2)
      nextra = nodes2[-1][-1].shape[0]
    der = np.zeros((nnodes + nextra, dim))
    for i in range(nnodes):
      der[i] = np.cross(axis,nodes[-1][-1][i] - origin) / scale
    write_xpost(fname, "Attributes", [0], der, header='MODE')
def displace_top_along_axis(fname, top, axes,s):
    #axes is some collection of displacement vectors
    nodes, elems   = read_top(top)
    displacement = s[0] * axes[0]
    if len(s) > 1:
      for i in range(1, len(s)):
        displacement += s[i] * axes[i]
    nodes[-1][-1] = nodes[-1][-1] + displacement
    write_top(fname, nodes, elems)
def displace_top_along_axis_derivs(fname, top, axes, top2 = None):
    #top2 is a second top file with which top is joined which means the xpost needs to be padded with zeros
    nodes, elems = read_top(top)
    nnodes, dim = nodes[-1][-1].shape
    nshapevar = len(axes)
    nextra = 0
    if top2 is not None:
      nodes2, elems2 = read_top(top2)
      nextra = nodes2[-1][-1].shape[0]
    der = np.zeros((nnodes + nextra, dim, nshapevar))
    for i in range(nshapevar):
      der[0:nnodes,:,i] = axes[i]
    write_xpost(fname, "Attributes", range(nshapevar), der, header='MODE')
def pad_derivs(fname, deriv_xpost, top2):
    #pads derivx_xpost with n zeros where n is the size of top2
    t, val = read_xpost(deriv_xpost)
    nnodes, dim, nders = val.shape
    nodes2, elems2 = read_top(top2)
    nextra = nodes2[-1][-1].shape[0]
    der = np.zeros((nnodes + nextra, dim, nders))
    for ider in range(nders):
        der[:nnodes,:,ider] = val[:,:,ider]
    write_xpost(fname, "Attributes", range(nders), der, header='MODE')
def scale_top(fname, top, scale,offsety = False):
    nodes, elems   = read_top(top)
    if scale is not None:
      nodes[-1][-1] *= scale
    nnodes = nodes[-1][-1].shape[0]
    if offsety:
      for i in range(nnodes):
        if nodes[-1][-1][i][1] < 0.001:
          nodes[-1][-1][i][1] = -0.001
    write_top(fname, nodes, elems)

def scale_top_custom(fname, top, y0, ycut):
    nodes, elems   = read_top(top)
    nnodes = nodes[-1][-1].shape[0]
    stretch = (ycut - y0) / (3.010 - y0)
    for i in range(nnodes):
      y = nodes[-1][-1][i][1]
      if y < y0:
        nodes[-1][-1][i][1] = stretch * (y - y0) + y0;
    write_top(fname, nodes, elems)

def update_top_with_vmo(fname, top, vmo):
    nodes, elems   = read_top(top)
    position       = read_vmo(vmo)
    nodes[-1][-1] = position
    write_top(fname, nodes, elems)

def displace_nodeset_with_vmo(fname, nodeset, vmo):
    nodes = read_nodeset(nodeset)+read_vmo(vmo)
    write_nodeset(fname, nodes)

def displaceABS_nodeset_with_vmo(fname, nodeset, vmo):

    nodes = read_nodeset(nodeset) + read_vmo(vmo)
    with open(vmo) as myfile:
        head = [next(myfile) for x in range(3)]

    f = open(fname,'w')
    for i in range(3):
        f.write(head[i])

    for k, node in enumerate(nodes):
        f.write(str(node[0])+' '+str(node[1])+' '+str(node[2])+'\n')
    f.close()


def displace_nodeset_with_vmo_numbered(fname, nodeset, vmo):
    nodes = read_nodeset(nodeset)+read_vmo(vmo)
    number_nodes = (np.loadtxt(nodeset, dtype='float', skiprows=2, comments='END')[:, 0]).astype(int)
    write_nodeset(fname,nodes,number_nodes)

def displace_nodeset_with_vmo1(fname, nodeset, vmo):
    nodes = read_nodeset(nodeset)+read_vmo(vmo)
    write_nodeset_1(fname, nodes)

def part_mesh(top, ndec, log=None, make_call=True, partmesh=None):

    # Default partnmesh executable
    if partmesh is None: partmesh = os.path.expandvars('$PARTMESH')

    # Build execution string, execute
    if make_call:
        top_old = copy.copy(top)
        top = top.split('/')[-1]+'.copy'
        subprocess.call('cp {0:s} {1:s}'.format(top_old, top), shell=True)
    exec_str = "{0:s} {1:s} {2:d}".format(partmesh, top, ndec)
    execute_code(exec_str, log, make_call)
    if make_call:
        subprocess.call('rm {0:s}'.format(top), shell=True)
        subprocess.call('mv {0:s}.dec.{2:d} {1:s}.dec.{2:d}'.format(top,
                                                     top_old, ndec), shell=True)
        top = copy.copy(top_old)
    return "{0:s}.dec.{1:d}".format(top, ndec)

def run_matcher(top, feminput, geom_prefix, log=None, make_call=True, matcher=None):

    # Default matcher executable
    if matcher is None: matcher = os.path.expandvars('$MATCHER')

    # Build execution string, execute
    exec_str = "{0:s} {1:s} {2:s} -output {3:s}".format(matcher, top, feminput,
                                                        geom_prefix)
    execute_code(exec_str, log, make_call)
    return "{0:s}.match.fluid".format(geom_prefix)

def run_matcher_special(top, feminput, geom_prefix,special, log=None, make_call=True, matcher=None):

    # Default matcher executable
    if matcher is None: matcher = os.path.expandvars('$MATCHER')

    # Build execution string, execute
    exec_str = "{0:s} {1:s} {2:s} {3:s} -output {4:s}".format(matcher, top, feminput,special,
                                                        geom_prefix)

    execute_code(exec_str, log, make_call)
    return "{0:s}.match.fluid".format(geom_prefix)


def sower_fluid_top(top, dec, cpus, nclust, geom_prefix, log=None,
                    make_call=True, sower=None, match=None):

    # Default sower executable
    if sower is None: sower = os.path.expandvars('$SOWER')

    # Build execution string, execute
    exec_str = "{0:s} -fluid -mesh {1:s} -dec {2:s}".format(sower, top, dec)
    for cpu in cpus:
        exec_str += " -cpu {0:d}".format(cpu)
    if match is not None:
        exec_str += " -match {0:s}".format(match)
    exec_str += " -cluster {0:d} -output {1:s}".format(nclust, geom_prefix)
    execute_code(exec_str, log, make_call)

def sower_fluid_extract_surf(msh, con, surf_top, bccode=-3, log=None,
                             make_call=True, sower=None):

    # Default sower executable
    if sower is None: sower = os.path.expandvars('$SOWER')

    # Build execution string, execute
    exec_str = ("{0:s} -fluid -merge -con {1:s} -mesh {2:s} "+
                "-skin {3:s} -bc {4:d}").format(sower, con, msh, surf_top,
                                                bccode)
    execute_code(exec_str, log, make_call)
    subprocess.call('mv {0:s}.xpost {0:s}'.format(surf_top), shell=True)

def sower_fluid_mesh_motion(mm_file, msh, con, out, bccode=-3, log=None,
                            make_call=True, sower=None):

    # Default sower executable
    if sower is None: sower = os.path.expandvars('$SOWER')

    # Build execution string, execute
    exec_str = ("{0:s} -fluid -split -con {1:s} -mesh {2:s} "+
                "-result {3:s} -ascii -bc {4:d} -out {5:s}").format(
                                          sower, con, msh, mm_file, bccode, out)
    execute_code(exec_str, log, make_call)

def sower_fluid_split(file2split, msh, con, out, from_ascii=True, log=None,
                      make_call=True, sower=None):

    # Default sower executable
    if sower is None: sower = os.path.expandvars('$SOWER')

    # Build execution string, execute
    exec_str = ("{0:s} -fluid -split -con {1:s} -mesh {2:s} "+
                "-result {3:s}").format(sower, con, msh, file2split)
    if from_ascii: exec_str = "{0:s} -ascii".format(exec_str)
    exec_str = "{0:s} -out {1:s}".format(exec_str, out)
    execute_code(exec_str, log, make_call)

def sower_fluid_merge(res_file, msh, con, out, name, from_bin=False, log=None,
                      make_call=True, sower=None):

    # Default sower executable
    if sower is None: sower = os.path.expandvars('$SOWER')

    # Build execution string, execute
    exec_str = "{0:s} -fluid -merge -con {1:s}".format(sower, con)
    exec_str = "{0:s} -mesh {1:s} -result {2:s}".format(exec_str, msh, res_file)
    exec_str = "{0:s} -name {1:s} -out {2:s}".format(exec_str, name, out)
    exec_str = "{0:s} -precision 15".format(exec_str)
    if from_bin: exec_str = "{0:s} -binary".format(exec_str)
    execute_code(exec_str, log, make_call)

def sower_fluid_merge_geom(msh, con, dec, from_bin=False, log=None,
                      make_call=True, sower=None):

    # Default sower executable
    if sower is None: sower = os.path.expandvars('$SOWER')

    # Build execution string, execute
    exec_str = "{0:s} -fluid -merge -con {1:s}".format(sower, con)
    exec_str = "{0:s} -mesh {1:s} -dec {2:s}".format(exec_str, msh, dec)
    if from_bin: exec_str = "{0:s} -binary".format(exec_str)
    execute_code(exec_str, log, make_call)

def run_cd2tet_fromtop(top, out, log=None, make_call=True, cd2tet=None):

    # Default cd2tet executable
    if cd2tet is None: cd2tet = os.path.expandvars('$CD2TET')

    # Build execution string, execute
    exec_str = "{0:s} -mesh {1:s} -output {2:s}".format(cd2tet, top, out)
    execute_code(exec_str, log, make_call)

def run_xp2exo(top, exo_out, xpost_in=[], log=None, make_call=True,
               xp2exo=None):

    # Default xp2exo executable
    if xp2exo is None: xp2exo = os.path.expandvars('$XP2EXO')

    # Build execution string, execute
    exec_str = "{0:s} {1:s} {2:s}".format(xp2exo, top, exo_out)
    for res in xpost_in:
        exec_str = "{0:s} {1:s}".format(exec_str, res)
    execute_code(exec_str, log, make_call)

def meshtools_plane(res_file, msh, con, out, plane, log=None, make_call=True,
                    meshtools=None):
    pass

def split_nodeset(fname, npartition):
    X = read_nodeset(fname)
    stride = np.ceil( float(X.shape[0])/float(npartition) )
    for k in range(npartition):
        fname_partk = '{0:s}.{1:d}parts.part{2:s}'.format(fname, npartition,
                                             str(k).zfill(len(str(npartition))))
        if k == 0:
            write_nodeset_1(fname_partk, X[:int(stride), :] )
        elif k == npartition-1:
            write_nodeset_1(fname_partk, X[int(k*stride):, :])
        else:
            write_nodeset_1(fname_partk, X[int(k*stride) :int((k+1)*stride), :])

def cat_split(fname_root, npartition, cleanup=False):
    s = ''
    for k in range(npartition):
        s += ' {0:s}.{1:d}parts.part{2:s}'.format(fname_root, npartition,
                                             str(k).zfill(len(str(npartition))))
    cat_str = 'cat {0:s} > {1:s}'.format(s, fname_root)
    execute_code(cat_str, None, True)
    if cleanup:
        rm_str  = 'rm {0:s}'.format(s)
        execute_code(rm_str, None, True)

def cat_split_gen(fname_root, fname_gen, n, cleanup=False):
    s = ''
    for k in range(n):
        s += ' {0:s}'.format(fname_gen(k))
    cat_str = 'cat {0:s} > {1:s}'.format(s, fname_root)
    execute_code(cat_str, None, True)
    if cleanup:
        rm_str  = 'rm {0:s}'.format(s)
        execute_code(rm_str, None, True)
