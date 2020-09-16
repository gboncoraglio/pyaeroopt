"""Make AEROF/S binaries from ASCII geometry files (.top, fem.input)."""

from pyaeroopt.test.ARW2.workflow import ARW2


def make_binaries():
    """ Makes the binaries for the ARW2 problem"""

    # Decomposition
    ndec = 32
    nclust = 32
    cpus = [32]

    # Make binaries
    ARW2.part_mesh(ndec, log='./binaries/binLog/partnmesh.log')
    ARW2.run_matcher(log='./binaries/binLog/match.log')
    ARW2.sower_fluid_top(cpus, nclust,
                         log='./binaries/binLog/sower.fluid.top.log')


if __name__ == '__main__':
    make_binaries()
