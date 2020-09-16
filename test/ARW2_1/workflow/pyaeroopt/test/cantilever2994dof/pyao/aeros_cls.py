
import numpy as np

from pyaeroopt.interface.aeros import Aeros, AerosInputFile

from pyaeroopt.test.cantilever2994dof import home_dir, src_dir
from pyaeroopt.test.cantilever2994dof.pyao.aeros_blk import *

def deformation_gradient(p):
    #F = np.array([[p[0], p[1], p[2]], [p[1], p[3], p[4]], [p[2], p[4], p[5]]])
    F = p.reshape((3, 3), order='F')
    return F

class AerosHdm(Aeros):
    """
    p : numpy array of size (3, 3)
        Deformation gradient
    desc : None
    desc_ext : None
    """
    def __init__(self, **kwargs):
        super(AerosHdm, self).__init__(**kwargs)

    def create_input_file(self, p, desc_ext=None, db=None):

        # Prefix for filenames
        p_str  = '_'.join(["{0:12.10f}".format(pi) for pi in p])
        prefix = "hdm_F{0:s}".format(p_str)
        prefix = prefix.replace('-', 'm')

        fname    = '{0:s}input/{1:s}.input'.format(home_dir, prefix)
        log      = '{0:s}log/{1:s}.log'.format(home_dir, prefix)
        idisp    = '{0:s}input/{1:s}.idisp'.format(home_dir, prefix)
        gdisplac = '"{0:s}output/{1:s}.gdisplac"'.format(home_dir, prefix)
        stressvm = '"{0:s}output/{1:s}.stressvm"'.format(home_dir, prefix)
        statvec   = '"{0:s}output/{1:s}.statvec"'.format(home_dir, prefix)

        output_hdm.gdisplac[2] = gdisplac
        output_hdm.stressvm[2] = stressvm
        disp.INCLUDE       = idisp

        self.infile = AerosInputFile(fname, [static, nonlin, output_hdm, topo,
                                             disp, forc], log)

        # Write displacement file
        F = deformation_gradient(p)
        pbndy = np.loadtxt(src_dir+'bndy.nodeset')
        f = open(idisp, 'w')
        f.write('DISPLACEMENTS\n')
        for i, ipbndy in enumerate(pbndy):
            pnew = np.dot(F, ipbndy[1:])
            for j in range(3):
                f.write("{0:d} {1:d} {2:16.12e}\n".format(int(ipbndy[0]), j+1,
                                                          pnew[j]))
        f.close()

class AerosSvd(Aeros):
    """
    p : numpy array of size (3, 3)
        Deformation gradient
    desc : None
    desc_ext : None
    """
    def __init__(self, **kwargs):
        super(AerosHdm, self).__init__(**kwargs)

    def create_input_file(self, p, desc_ext=None, db=None):

        # Prefix for filenames
        p_str  = '_'.join(["{0:12.10f}".format(pi) for pi in p])
        prefix = "svd_F{0:s}".format(p_str)
        prefix = prefix.replace('-', 'm')

        fname    = '{0:s}input/{1:s}.input'.format(home_dir, prefix)
        log      = '{0:s}log/{1:s}.log'.format(home_dir, prefix)
        statvec   = '"{0:s}output/{1:s}.statvec"'.format(home_dir, prefix)
        robdataf  = '"{0:s}output/{1:s}.robdataf"'.format(home_dir, prefix)

        robc.snapfi = statvec.replace('svd', 'hdm')
        output_svd.robdataf[2] = robdataf

        self.infile = AerosInputFile(fname, [static, nonlin, robc, output_svd,
                                             topo, disp, forc], log)
        f.close()

class AerosRom(Aeros):
    """
    p : numpy array of size (3, 3)
        Deformation gradient
    desc : None
    desc_ext : None
    """
    def __init__(self, **kwargs):
        super(AerosHdm, self).__init__(**kwargs)

    def create_input_file(self, p, desc_ext=None, db=None):

        # Prefix for filenames
        p_str  = '_'.join(["{0:12.10f}".format(pi) for pi in p])
        prefix = "rom_F{0:s}".format(p_str)
        prefix = prefix.replace('-', 'm')

        fname    = '{0:s}input/{1:s}.input'.format(home_dir, prefix)
        log      = '{0:s}log/{1:s}.log'.format(home_dir, prefix)
        gdisplac = '"{0:s}output/{1:s}.gdisplac"'.format(home_dir, prefix)
        stressvm = '"{0:s}output/{1:s}.stressvm"'.format(home_dir, prefix)
        robdataf  = '"{0:s}output/{1:s}.robdataf"'.format(home_dir, prefix)

        readmode.READMODE[0] = robdataf.replace('rom', 'hdm')
        output_rom.gdisplac[2] = gdisplac
        output_rom.stressvm[2] = stressvm

        self.infile = AerosInputFile(fname, [static, nonlin, readmode,
                                             output_rom, topo, disp, forc], log)
        f.close()
