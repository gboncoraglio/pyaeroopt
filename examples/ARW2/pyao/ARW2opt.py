""" Define ARW2-problem-specific classes."""

import numpy as np

import matlab.engine

from pyaeroopt.test.ARW2.pyao.aerof_cls import AerofAeroelasticSensitivity
from pyaeroopt.test.ARW2.pyao.aeros_cls import AerosAeroelasticSensitivity
from pyaeroopt.test.ARW2.pyao.sdesign_cls import (StructureDeformer,
                                                  FluidSkinDeformer,
                                                  MaterialFileChanger)
from pyaeroopt.interface.aerofs import Aerofs
from pyaeroopt.test.ARW2.workflow import (hpcCoupled, hpc1, ROM_file,
                                          ROM_interpolation_folder)

class ARW2Problem(object):

    """Define the optimization problem (variables, bounds, objective, constraint,
     gradients).
     """

    def __init__(self, load_matlab=True):
        self.designpoint = None
        self.num_interpolations = 0
        self.num_design_points = 0
        self.num_aeroelastic_sims = 0
        self.last_obj_design_point = None
        self.last_obj_value = None
        self.last_obj_gradient = None
        self.last_constraint_design_point = None
        self.last_weight_value = None
        self.last_weight_gradient = None
        self.last_stress_value = None
        self.last_stress_gradient = None
        self.last_damping_value = None
        self.last_damping_gradient = None

        # Launch MATLAB and load ROM interpolation variables into workspace
        if load_matlab:
            self.eng = matlab.engine.start_matlab()
            self.eng.addpath(self.eng.genpath(ROM_interpolation_folder),
                             nargout=0)
            self.eng.load(ROM_file, nargout=0)

    def set_design_point(self, x):
        """Moves meshes and alters input files for a new design point."""
        x = np.array(x, dtype='float')
        if np.array(x != self.designpoint).any():
            # Move the meshes and update the material file
            StructureDeformer(exec_dir='Sdesign/').execute(x, hpc=hpc1)
            FluidSkinDeformer(exec_dir='Sdesign/').execute(x, hpc=hpc1)
            MaterialFileChanger(exec_dir='AeroS-Files').execute(x, hpc=hpc1)

            # Update instance variables and print current design point
            self.designpoint = [float(i) for i in x]
            self.num_design_points += 1
            printstring = '...Setting design point:\t [\t'
            for i in x[:-1]:
                printstring += str(i) + ',\t '
            print(printstring + str(x[-1]) + '\t]\n')

    def evaluate_objective(self, x):
        """ Computes lift and drag and sets the objective and objective gradient
        instance variables.
        """
        x = np.array(x, dtype='float')    # cast x as a numpy float array
        if np.array(x == self.last_obj_design_point).all():
            return dict([('Objective', self.last_obj_value),
                         ('ObjectiveGradient', self.last_obj_gradient)])
        else:
            self.run_simulation(x)
            print('...Current L/D:\t' + str(-self.last_obj_value))
            return dict([('Objective', self.last_obj_value),
                         ('ObjectiveGradient', self.last_obj_gradient)])

    def run_simulation(self, x):
        """ Run an aeroelastic sensitivity calculation at design point x."""
        self.set_design_point(x)
        Aerofs(AerofAeroelasticSensitivity(),
               AerosAeroelasticSensitivity()).execute(x, hpc=hpcCoupled)
        self.set_objective()
        self.last_obj_design_point = x
        self.num_aeroelastic_sims += 1

    def set_objective(self):
        """ Extract lift and drag and gradients from the current output files
            and places the results into this object's instance variables.
        """
        data = np.genfromtxt('out/ARW2.sensitivity.lift.drag', delimiter=' ',
                             skip_header=1)
        lift = data[0, 4]
        drag = data[0, 2]
        self.last_obj_value = -lift/drag
        lift_sens = data[:, 7]
        lift_sens[3:] /= 10
        drag_sens = data[:, 5]
        drag_sens[3:] /= 10
        self.last_obj_gradient = -(lift_sens/drag - lift/(drag**2) * drag_sens)
        return

    def evaluate_constraint(self, x):
        """Evaluates all of the constraints and returns their values in a
        dictionary.

        Returns:
        ----------
        returnDict: dictionary
        Contains the weight in the field 'Weight', the weight-gradient in the field
        'WeightGradient', the stress in the field 'Stress', the stress-gradient
        in the field 'StressGradient', the damping ratios in the field
        'Damping', and the damping ratio gradients in the field
        'DampingGradient'.
        """
        x = np.array(x, dtype='float')
        if np.array(x == self.last_constraint_design_point).all():
            return dict([('Weight', self.last_weight_value),
                         ('WeightGradient', self.last_weight_gradient),
                         ('Stress', self.last_stress_value),
                         ('StressGradient', self.last_stress_gradient),
                         ('Damping', self.last_damping_value),
                         ('DampingGradient', self.last_damping_gradient)])
        else:
            if np.array(x != self.last_obj_design_point).any():
                self.run_simulation(x)
            interp_res = self.interpolate_flutter(x)
            self.set_constraints(interp_res)
            print('...Current Weight:  ' + str(self.last_weight_value) + '\n',
                  '...Current Stress:  ' + str(max(self.last_stress_value)) + '\n',
                  '...Current Damping: ' + str(np.amin(self.last_damping_value)) + '\n')
            return dict([('Weight', self.last_weight_value),
                         ('WeightGradient', self.last_weight_gradient),
                         ('Stress', self.last_stress_value),
                         ('StressGradient', self.last_stress_gradient),
                         ('Damping', self.last_damping_value),
                         ('DampingGradient', self.last_damping_gradient)])

    def set_constraints(self, interp_res):
        """Extracts constraint values from output files and places them into
        instance variables.
        """
        self.set_weight_constraint()
        self.set_stress_constraint()
        self.set_damping_constraint(interp_res)

    def set_weight_constraint(self):
        """Extracts the weight and weight gradients and places them into
        instance variables.
        """
        weight = np.genfromtxt('out/ARW2.weigshap', delimiter=' ',
                               skip_header=0, skip_footer=3)
        weight_shape = np.genfromtxt('out/ARW2.weigshap', delimiter=' ',
                                     skip_header=1)
        weight_thickness = np.genfromtxt('out/ARW2.weigthic', delimiter=' ',
                                         skip_header=1)
        weight_grad = np.append(weight_shape, weight_thickness).reshape([1, 6])
        self.last_weight_value = np.array([float(weight)])
        self.last_weight_gradient = weight_grad

    def set_stress_constraint(self):
        """Extracts the stress and stress gradients and places them into
        instance variables.
        """
        vm_stress = np.genfromtxt('out/ARW2.stressvm', delimiter=' ',
                                  skip_header=460)
        vm_shape = np.genfromtxt('out/ARW2.vmstshap', delimiter=',',
                                 skip_header=1)
        vm_thickness = np.genfromtxt('out/ARW2.vmstthic', delimiter=',',
                                     skip_header=1)
        vm_grad = np.append(vm_shape, vm_thickness, axis=1)
        self.last_stress_value = vm_stress
        self.last_stress_gradient = vm_grad

    def set_damping_constraint(self, interp_res):
        """Sets the damping ratio and gradient instance variables given the
        result of the ROM interpolation as input.
        """
        self.last_damping_value = (np.array(interp_res[0])).reshape([6, 1])[:, 0]
        self.last_damping_gradient = np.transpose(np.array(interp_res[1]))

    def interpolate_flutter(self, x):
        """ Estimates the flutter eigenvalues and their gradients by interpolating
            the database of ROMs loaded into 'StructROM_rotated_c' and
            'AeroelasticROM_rotated_c'.  Uses a MATLAB engine to interpolate.
        """
        p = np.reshape(x, [6, 1]).tolist() # the interpolant requires a column vector
        self.eng.workspace['p'] = matlab.double(p)
        interp_res = self.eng.eval('runFlutterAnalysis(p,'
                                   'StructROM_rotated_c,'
                                   'AeroelasticROM_rotated_c)',
                                   nargout=2)
        self.last_constraint_design_point = x
        self.num_interpolations += 1
        return interp_res
