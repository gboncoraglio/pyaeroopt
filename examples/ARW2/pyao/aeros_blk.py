"""Exports many Aero-S input blocks, to help create AerosInputFile objects.

Note: For expediency, some input file fields are exported from this module
with numeric values already filled in, and these values are not necessarily
the same as the aeros defaults.  CHECK ALL VALUES FOR YOUR PROBLEM.
"""

from pyaeroopt.interface.aeros import AerosInputBlock

controlEig = AerosInputBlock('CONTROL',
                             ['ARW2.EIG'],
                             ['1'],
                             ['NodSet'],
                             ['EleSet'])

controlMPP = AerosInputBlock('CONTROL',
                             ['ARW2.DYNAM'],
                             ['1'],
                             ['NodSet'],
                             ['EleSet'])

controlSens = AerosInputBlock('CONTROL',
                              ['ARW2.SENSITIVITY'],
                              ['1'],
                              ['NodSet'],
                              ['EleSet'])

controlSteady = AerosInputBlock('CONTROL',
                              ['ARW2.STEADY'],
                              ['1'],
                              ['NodSet'],
                              ['EleSet'])

static = AerosInputBlock('STATIC', 'skyline')

qstatic = AerosInputBlock('QSTATICS',
                          ['mech', 0.5, 1.0e-6, 2000, 0.0])

renum = AerosInputBlock('RENUM', 'rcm')


eig = AerosInputBlock('EIGEN',
                      ['nsbspv', 14],
                      ['neigpa', 6],
                      ['toleig', '1e-6'],
                      ['toljac', '1e-4'],
                      ['maxiter', 200],
                      ['qrfactorization', 'on'])

dyn = AerosInputBlock('DYNAMICS',
                      ['newmark'],
                      ['mech', 0.25, 0.5],
                      ['time', 0.0, 0.0001, 0.001],
                      ['MODE'])

aero = AerosInputBlock('AERO',
                       ['MPP', 0.1],
                       ['A6', None],
                       ['Matcher', '"binaries/ARW2.match.fem"'],
                       ['READMODES', '"Sresults/ARW2.struct_modes_input"'])



outputEig = AerosInputBlock('OUTPUT',
                            ['geigenpa', '"Sresults/ARW2.struct_modes" 1'],
                            ['xmatrix', '"Sresults/ARW2.xmatrix"'],
                            ['qmatrix', '"Sresults/ARW2.qmatrix"'],
                            ['rmatrix', '"Sresults/ARW2.rmatrix"'])

outputSens = AerosInputBlock('OUTPUT',
                             ['stressvm', 20, 15, '"out/ARW2.stressvm"', 2001, 'upper'],
                             ['vmstshap', 20, 15, '"out/ARW2.vmstshap"', 1, 'upper'],
                             ['vmstthic', 20, 15, '"out/ARW2.vmstthic"', 1, 'upper'],
                             ['weigshap', 20, 15, '"out/ARW2.weigshap"', 1],
                             ['weigthic', 20, 15, '"out/ARW2.weigthic"', 1])

output6 = AerosInputBlock('OUTPUT6',
                          ['geigenpa', '"Sresults/ARW2.struct_modes6" 1'])


# Include statements for ARW2 properties
sdes = AerosInputBlock('INCLUDE_SDESIGN',
                       ['INCLUDE', '"Sdesign/position.nodefile"'])


top = AerosInputBlock('INCLUDE_TOPO',
                      ['INCLUDE', '"AeroS-Files/topologyfile"'])

cframe = AerosInputBlock('INCLUDE_CFRAME',
                         ['INCLUDE', '"AeroS-Files/cframefile"'])

disp = AerosInputBlock('INCLUDE_DISPLACEMENT',
                       ['INCLUDE', '"AeroS-Files/displacementfile"'])

dim = AerosInputBlock('INCLUDE_DIMASS',
                      ['INCLUDE', '"AeroS-Files/dimassfile"'])

eframe = AerosInputBlock('INCLUDE_EFRAME',
                         ['INCLUDE', '"AeroS-Files/eframefile"'])

comp = AerosInputBlock('INCLUDE_COMPSITES',
                       ['INCLUDE', '"AeroS-Files/compositefile"'])

mat = AerosInputBlock('INCLUDE_MATERIAL',
                      ['INCLUDE', '"AeroS-Files/materialfile"'])

att = AerosInputBlock('INCLUDE_ATRRIBUTES',
                      ['INCLUDE', '"AeroS-Files/attributefile"'])

group = AerosInputBlock('INCLUDE_GROUP',
                        ['INCLUDE', '"AeroS-Files/groupfile"'])

grav = AerosInputBlock('GRAVITY',
                       ['', 0, 0, -386.088])

sens = AerosInputBlock('SENSITIVITY',
                       ['method', 'Direct'],
                       ['tolsen', 1e-6],
                       ['relaxationsen', 0.5],
                       #['relaxationsen', None],
                       ['readsensitivity', '"Sdesign/position.der"'],
                       ['thgrli', 1, 2, 3],
                       ['ksparam', 50])
