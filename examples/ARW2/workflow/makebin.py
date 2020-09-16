"""Make AEROF/S binaries from ASCII geometry files (.top, fem.input)."""

import sys

from pyaeroopt.test.ARW2.workflow import ARW2, hpc


def make_binaries(n_cpus=None):
    """ Makes the binaries for the ARW2 problem"""
    if n_cpus is None:
        n_cpus = hpc.nproc

    # Decomposition
    ndec = n_cpus
    nclust = n_cpus
    cpus = [n_cpus]

    # Make binaries
    ARW2.part_mesh(ndec, log='./binaries/binLog/partnmesh.log')
    ARW2.run_matcher(log='./binaries/binLog/match.log')
    ARW2.sower_fluid_top(cpus, nclust,
                         log='./binaries/binLog/sower.fluid.top.log')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        make_binaries(int(sys.argv[1]))
    else:
        make_binaries()
