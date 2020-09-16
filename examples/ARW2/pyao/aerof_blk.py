"""Exports many Aero-F input blocks, to help create Aero-F InputFile objects.

Note: For expediency, some input file fields are exported from this module
with numeric values already filled in, and these values are not necessarily
the same as the aerof defaults.  CHECK ALL VALUES FOR YOUR PROBLEM.
"""

import sys, os

# PYAEROOPT
from pyaeroopt.interface.aerof import AerofInputBlock

## Problem
prob = AerofInputBlock('Problem',
                       ['Type', None],
                       ['Mode', 'Dimensional'],
                       ['Prec', None])

## Input
inpu = AerofInputBlock('Input',
                       ['Prefix', ""],
                       ['GeometryPrefix', None],
                       ['ShapeDerivative', None],
                       ['MultipleSolutions', None],
                       ['OptimalPressure', None],
                       ['OptimalPressureDimensionality', None],
                       ['StateSnapshotData', None],
                       ['SensitivitySnapshotData', None],
                       ['StateSnapshotReferenceSolution', None],
                       ['Solution', None],
                       ['PODData', None],
                       ['Matcher', None],
                       ['InitialWallDisplacement', None])

## Output
post = AerofInputBlock('Postpro',
                       ['Prefix', ''],
                       ['LiftandDrag', None],
                       ['Force', None],
                       ['MatchPressure', None],
                       ['LiftandDragSensitivity', None],
                       ['MatchPressureSensitivity', None],
                       ['PressureSensitivity', None],
                       ['StateVectorSensitivity', None],
                       ['FluxPartialSensitivity', None],
                       ['Pressure', None],
                       ['Displacement', None],
                       ['PODData', None],
                       ['ROM', None],
                       ['AeroelasticEigenvalues', None],
                       ['Frequency', 0])

nonlRom = AerofInputBlock('NonlinearROM',
                          ['Prefix', ''],
                          ['StateVector', None],
                          ['StateVectorOutputFrequencyTime', 200],
                          ['StateVectorOutputFrequencyNewton', 0],
                          ['FluxNormSensitivity', None])

rest = AerofInputBlock('Restart',
                       ['Prefix', 'Fresults/Restart/'],
                       ['Solution', 'sol'],
                       ['Position', None],
                       ['RestartData', None],
                       ['Frequency', 0])

outp = AerofInputBlock('Output',
                       ['Postpro', post],
                       ['Restart', rest])

## Linearized computations
lin = AerofInputBlock('Linearized',
                      ['Amplification', 0.1],
                      ['Eps', 1e-4],
                      ['ToleranceEigenAeroelastic', None],
                      ['MaxItsEigenAeroelastic', None],
                      ['NumStrModes', 6],
                      ['StrModes', None],
                      ['ExcMode', 1],
                      ['Domain', 'Frequency'],
                      ['FreqStep', None],
                      ['NumPOD', 100])


## NonlinearRomFileSystem
nonlRomDire = AerofInputBlock('Directories',
                              ['Prefix', ''],
                              ['TopLevelDirectoryName', None],
                              ['ClusterDirectoryName', None],
                              ['SensitivityClusterDirectoryName', None])

nonlRomFile = AerofInputBlock('Files',
                              ['StatePrefix', None],
                              ['StateBasisPrefix', None],
                              ['SensitivityPrefix', None],
                              ['SensitivityBasisPrefix', None],
                              ['ResidualPrefix', None],
                              ['GNATPrefix', None])

nonlRomFileSyst=AerofInputBlock('NonlinearRomFileSystem',
               ['NumClusters', 1],
               ['Directories', nonlRomDire],
               ['Files', nonlRomFile])

## NonlinearRomOffline
# Clustering
clus=AerofInputBlock('Clustering',
               ['ClusteringAlgorithm', 'KMeansWithBounds'],
               ['PercentOverlap', 1.0],
               ['KMeansMaxIterations', 1],
               ['KMeansMaxAggressiveIterations', 1],
               ['KMeansRandomSeed',10101],
               ['MinClusterSize', 1],
               ['UseExistingClusters', False])

# State snapshots
stSnap=AerofInputBlock('Snapshots',
               ['NormalizeSnaps', None],
               ['IncrementalSnaps', False],
               ['SubtractNearestSnapToCenter', False],
               ['SubtractClusterCenters', None],
               ['SubtractReferenceState', None])

stDataComp=AerofInputBlock('DataCompression',
               ['ComputePOD', True],
               ['PODMethod', 'ScalapackSVD'],
               ['SingularValueTolerance', 1e-16],
               ['MaxBasisSize', 50000],
               ['MinBasisSize', 10],
               ['MaxEnergyRetained', 1.0])

stRob=AerofInputBlock('StateROB',
               ['Snapshots', stSnap],
               ['DataCompression', stDataComp])


# Sensitivity snapshots
sensSnap=AerofInputBlock('Snapshots',
               ['NormalizeSnaps', False])

sensDataComp=AerofInputBlock('DataCompression',
               ['ComputePOD', True],
               ['PODMethod', 'ScalapackSVD'],
               ['SingularValueTolerance', 1e-16],
               ['MaxBasisSize', 50000],
               ['MinBasisSize', 10],
               ['MaxEnergyRetained', 1.0])

sensRob=AerofInputBlock('SensitivityROB',
               ['Snapshots', sensSnap],
               ['DataCompression', sensDataComp])


# Residual snapshots
resSnap=AerofInputBlock('Snapshots',
               ['NormalizeSnaps', None])

resDataComp=AerofInputBlock('DataCompression',
               ['ComputePOD', True],
               ['PODMethod', 'ScalapackSVD'],
               ['SingularValueTolerance', 1e-16],
               ['MaxBasisSize', 10000],
               ['MinBasisSize', 10],
               ['MaxEnergyRetained', 1.0])

resRob=AerofInputBlock('ResidualROB',
               ['Snapshots', resSnap],
               ['DataCompression', resDataComp])

onliBasiUpda=AerofInputBlock('OnlineBasisUpdates',
               ['PreprocessForNoUpdates', 'On'],
               ['PreprocessForFastExactUpdates', 'On'])

# Construct ROB
consROB=AerofInputBlock('ConstructROB',
               ['Clustering', clus],
               ['StateROB', stRob],
               ['SensitivityROB', sensRob],
               ['ResidualROB', resRob],
               ['OnlineBasisUpdates', onliBasiUpda])

# Construct GNAT
consGNAT=AerofInputBlock('ConstructGNAT',
               ['UseUnionOfSampledNodes', None],
               ['MaxDimensionStateROB', None],
               ['MinDimensionStateROB', None],
               ['EnergyStateROB', None],
               ['MaxDimensionResidualROB', None],
               ['MinDimensionResidualROB', None],
               ['EnergyResidualROB', None],
               ['ROBGreedy', 'Residual'],
               ['MaxDimensionROBGreedy', None],
               ['MinDimensionROBGreedy', None],
               ['ROBGreedyFactor', None],
               ['MaxSampledNodes', None],
               ['MinSampledNodes', None],
               ['SampledNodesFactor',None],
               ['NumPseudoInvNodesAtATime', 1000000],
               ['OutputReducedBases', 'True'])

nonlRomOffl=AerofInputBlock('NonlinearRomOffline',
               ['ConstructROB', consROB],
               ['ConstructGNAT', consGNAT])

## NonlinearRomOnline
sensOnli=AerofInputBlock('Sensitivities',
               ['Include', 'On'],
               ['GramSchmidt', 'On'],
               ['MaximumDimension', 50000],
               ['MinimumDimension', 1],
               ['MaximumEnergy', None])

# NonlinearRomOnline
nonlRomOnli=AerofInputBlock('NonlinearRomOnline',
               ['Projection', 'PetrovGalerkin'],
               ['PerformLineSearch', None], # 'Backtracking'],
               ['MinimumDimension', 5],
               ['MaximumDimension', 50000],
               ['MaximumEnergy', None],
               ['FastDistanceComparisons', 'Off'],
               ['BasisUpdates', None],
               ['ProjectSwitchStateOntoAffineSubspace', None],
               ['StoreAllClustersInMemory', False],
               ['Sensitivities', sensOnli])

## SensitivityAnalysis
precSa=AerofInputBlock('Preconditioner',
               ['Type', 'Ras'],
               ['Fill', 0])

lineSolvSa=AerofInputBlock('LinearSolver',
               ['Type', 'Gmres'],
               ['MaxIts', 2000],
               ['KrylovVectors', 2000],
               ['Eps', 1e-10],
               ['Output', '"stdout"'],
               ['Preconditioner', precSa])

sensAnal=AerofInputBlock('SensitivityAnalysis',
               ['Method', 'Direct'],
               ['SensitivityComputation', 'Analytical'],
               ['SensitivityMesh', 'On'],
               ['SensitivityFSI' , None],
               ['AdaptiveEpsFSI' , None],
               ['SparseApproach' , 'On'],
               ['LinearSolver', lineSolvSa])

## MeshMotion
precMm=AerofInputBlock('Preconditioner',
               ['Type', 'Jacobi'],
               ['Fill', 0])

lineSolvMm=AerofInputBlock('LinearSolver',
               ['Type', 'Cg'],
               ['MaxIts', 5000],
               ['KrylovVectors', 1000],
               ['Eps', 1e-9],
               ['Preconditioner', precMm])

newtMm=AerofInputBlock('Newton',
               ['MaxIts', 1],
               ['Eps', 0.01],
               ['LinearSolver', lineSolvMm])

meshMoti=AerofInputBlock('MeshMotion',
               ['Type', 'Basic'],
               ['Element', 'TorsionalSprings'],
               ['Mode', 'NonRecursive'],
               ['NumIncrements', 1],
               ['Newton', newtMm])

## Surfaces
surfData=AerofInputBlock('SurfaceData[2]',
               ['Nx', 0.0],
               ['Ny', 0.0],
               ['Nz', 1.0])

surf=AerofInputBlock('Surfaces',
               ['SurfaceData[2]', surfData])

## Equations
FluidModel=AerofInputBlock('FluidModel[0]',
               ['Fluid', 'PerfectGas'])

ViscosityModel=AerofInputBlock('ViscosityModel',
               ['Type', 'Constant'],
               ['DynamicViscosity', 1e-5])

equa=AerofInputBlock('Equations',
               ['Type', 'Euler'],
               ['FluidModel', FluidModel],
               ['ViscosityModel', ViscosityModel])

## Reference State
refState=AerofInputBlock('ReferenceState',
               ['Length', 1.0],
               ['Mach', 0.8],
               ['Reynolds', 1.0e3],
               ['Temperature', 273.15])


## BoundaryConditions
inle=AerofInputBlock('Inlet',
               ['Mach', 0.8],
               ['Alpha', 0.0],
               ['Beta', 0.0],
               ['Pressure', 12.7162],
               ['Density', 1.0193e-7])

wall=AerofInputBlock('Wall',
               ['Type', 'Adiabatic'])

bounCond=AerofInputBlock('BoundaryConditions',
               ['Inlet', inle],
               ['Wall' , wall])

## Space
naviStokSp=AerofInputBlock('NavierStokes',
               ['Flux', 'Roe'],
               ['Reconstruction', 'Linear'],
               ['AdvectiveOperator', 'FiniteVolume'],
               ['Limiter', 'VanAlbada'],
               ['Gradient', 'LeastSquares'],
               ['Dissipation', 'SecondOrder'],
               ['Beta', 0.6666666666666],
               ['Gamma', 1.0])

spac=AerofInputBlock('Space',
               ['NavierStokes', naviStokSp])

## Time
precTi=AerofInputBlock('Preconditioner',
               ['Type', 'Ras'],
               ['Fill', 0])

naviStokTi=AerofInputBlock('NavierStokes',
               ['Type', 'Gmres'],
               ['EpsFormula', None],
               ['MaxIts', 100],
               ['KrylovVectors', 100],
               ['Eps', 0.001],
               ['Preconditioner', precTi])

lineSolvTi=AerofInputBlock('LinearSolver',
               ['NavierStokes', naviStokTi])

newtTi=AerofInputBlock('Newton',
               ['MaxIts', 100],
               ['Eps', 1e-03],
               ['EpsAbs', None],
               ['FailSafe', 'On'],
               ['LinearSolver', lineSolvTi])

impl=AerofInputBlock('Implicit',
               ['Type', 'ThreePointBackwardDifference'],
               ['MatrixVectorProduct', 'Exact'],
               ['Newton', newtTi])

CflLaw=AerofInputBlock('CflLaw',
               ['Strategy', 'Hybrid'],
               ['Cfl0', 1.0],
               ['CflMax', 1000])

time=AerofInputBlock('Time',
               ['Type', 'Implicit'],
               ['Form', None],
               ['MaxIts', 10000],
               ['Eps', 1e-12],
               ['TypeTimeStep', None],
               ['Implicit', impl],
               ['CflLaw', CflLaw])
