################################################################################
# Import pyAeroOpt into the script.                                            #
################################################################################
import pyaeroopt


################################################################################
# Create an Frg object to describe the problem geometry (.top file, FEM input  #
# file, binary prefixes, etc.) and an Hpc object to describe the computational #
# environment (cluster name, number of cores, etc).  Then create the binaries  #
# for the problem using methods of the Frg object.                             #
################################################################################
frg = pyaeroopt.interface.Frg(top="../mesh/top",
                              geom_pre="../mesh/binary/naca0012")  
hpc = pyaeroopt.interface.Hpc(machine="independence", 
                              mpi="mpiexec",
                              batch=False,
                              bg=False,
                              nproc=12)
frg.part_mesh(hpc.ppn)                               # part the mesh
frg.sower_fluid_top([hpc.ppn], hpc.ppn)              # sower the top file 


###############################################################################
# Create objects to describe the simulations that need to be run.  For this   #
# problem we need to run a steady Aero-F simulation and calculate the lift,   #  
# so we will use the Aerof class from pyAeroOpt. While the entire input file  #
# could be passed to the Aerof constructor, it is usually easier to create a  #
# subclass for each simulation, so we will do this.  The input file           #
# constructor takes in a list of AerofInputBlocks, which we define one at a   #
# time.  Note that any input file block can be defined as a function of the   #
# optimization variable(s).                                                   #
###############################################################################
from pyaeroopt.interface import Aerof, AerofInputFile, AerofInputBlock
class NACA_Steady(Aerof):

    # Note that the input file is a function of optimization variable(s)
    def create_input_file(self, p):
        """Define the input file for this problem."""
        mach   = p[0]      # this simulation is parameterized by Mach number
        alpha  = p[1]      # and angle of attack

	# Define the input file one block at a time. Note that the inlet
        # block is a function of the optimization variable
        prob = AerofInputBlock('Problem', ['Type', 'Steady'],
                                          ['Mode', 'Dimensional'])
        inp  = AerofInputBlock('Input',   ['GeometryPrefix', frg.geom_pre])
        rest = AerofInputBlock('Restart',
                               ['Solution', 'data/naca0012.sol'])
        post = AerofInputBlock('Postpro',
                               ['LiftandDrag', 'data/liftdrag'],
                               ['Pressure', 'data/pressure'])
        outp = AerofInputBlock('Output',  ['Postpro', post],
                                          ['Restart', rest])
        model= AerofInputBlock('FluidModel[0]', ['Fluid', 'PerfectGas'])        
        equa = AerofInputBlock('Equations', ['Type', 'Euler'],
                                            ['FluidModel[0]', model])
        wall = AerofInputBlock('Wall', ['Type', 'Adiabatic'])
        inlet= AerofInputBlock('Inlet', ['Mach', mach], 
                                        ['Beta', alpha],  # this mesh is oriented so
                                        ['Alpha', 0],     # that AoA is actually Beta
                                        ['Pressure', 12.71],
                                        ['Density', 1.0193e-07])
        bc   = AerofInputBlock('BoundaryConditions', ['Inlet', inlet],
                                                     ['Wall', wall])
        nav  = AerofInputBlock('NavierStokes', ['Reconstruction', 'Constant'])
        space= AerofInputBlock('Space', ['NavierStokes', nav])
        prec = AerofInputBlock('Preconditioner', ['Type', 'Ras'], ['Fill', 0])
        navNewton = AerofInputBlock('NavierStokes', ['MaxIts', 100],
                                    ['KrylovVectors', 100], ['Eps', 0.001],
                                    ['Preconditioner', prec])
        navLin = AerofInputBlock('LinearSolver', ['NavierStokes', navNewton])
        newton = AerofInputBlock('Newton', ['MaxIts', 200], ['Eps', 0.001],
                                           ['LinearSolver', navLin])
        impl = AerofInputBlock('Implicit', ['Newton', newton])
        cfl  = AerofInputBlock('CflLaw', ['Strategy', 'Hybrid'], ['Cfl0', 1.],
                                         ['CflMax', 1000])
        time = AerofInputBlock('Time', ['MaxIts', 1000], ['Eps', 1e-8],
                                       ['Implicit', impl], ['CflLaw', cfl])
        
        fname = "naca.input"
        log   = "naca.log"
        self.infile = AerofInputFile(fname, [prob, inp, outp, equa, 
                                             bc, space, time], log)


###############################################################################
# Set up the optimization problem                                             #
###############################################################################
import numpy, math
def getLD(alpha):
    NACA_Steady().execute([0.5, alpha], hpc=hpc)
    D0, L0 = numpy.loadtxt('data/liftdrag', usecols=[4,5], unpack=True)
    # rotate lift and drag into the proper frame
    Df = D0[-1];   Lf = L0[-1];   alpha_r = math.radians(-alpha)
    D = Df*math.cos(alpha_r) - Lf*math.sin(alpha_r)
    L = Df*math.sin(alpha_r) + Lf*math.cos(alpha_r) 
    print("...At alpha = {:7.8f}, \t L/D is {}".format(alpha, L/D))
    return -L/D

optimizer = pyaeroopt.opt.Optimize(name="NACA 0012")
optimizer.add_variables(1, [5], [-90], [90])
optimizer.add_objective(lambda alpha: getLD(alpha[0]), name="L/D")
optimizer.optimize('scipy:L-BFGS-B', options={"ftol": 0.01})


###############################################################################
# Create a .exo file for the pressure field so we can visualize the results   #
###############################################################################
frg.sower_fluid_merge('data/pressure', 'pressure', 'pressure')
frg.run_xp2exo('pressure.exo', ['pressure.xpost'])

