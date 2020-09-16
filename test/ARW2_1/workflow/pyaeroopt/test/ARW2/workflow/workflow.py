"""Run the optimization problem.

Usage: from the command line, type
        # python workflow.py 
    If you'd like to run the full optimization problem.
    If you'd like only evaluate the objective function, type
        # python workflow objective x1 x2 x3 x4 x5 x6
    where (x1,x2,x3,x4,x5,x6) is the design point to analyze
"""

import warnings
import pdb
import sys

import numpy as np
import os

from pyaeroopt.test.ARW2.pyao.ARW2opt import ARW2Problem
from pyaeroopt.opt import Optimize as Optimize

# Disable warnings about features that are deprecated in Python 3.5
# (PyAeroOpt MUST be run in Python 3.4 so this warning is irrelevant)
warnings.simplefilter(action='ignore', category=FutureWarning)

def setup_optimizer(design_point=None):
    """Run the optimization problem."""

    if design_point is None:
        design_point = [0, 0, 0, 0, 0, 0] # default 

    # Define the optimization problem
    opt_prob = ARW2Problem()

    # Set up the optimizer by passing in methods defined in the opt_prob object
    optimizer = Optimize(name='ARW2')
    optimizer.add_variables(6, design_point,
                            np.ones(6) * (-0.1),
                            np.ones(6) *  (0.1))
    optimizer.add_objective(lambda x: opt_prob.evaluate_objective(x)['Objective'],
                            name='Lift/Drag',
                            grad=lambda x: opt_prob.evaluate_objective(x)['ObjectiveGradient'])
    optimizer.add_nonlin_constraints(456,
                                     lambda x: opt_prob.evaluate_constraint(x)['Stress'],
                                     lambda x: opt_prob.evaluate_constraint(x)['StressGradient'],
                                     optimizer.minf*np.ones(456),
                                     25000*np.ones(456),
                                     name='Stress')
    optimizer.add_nonlin_constraints(1,
                                     lambda x: opt_prob.evaluate_constraint(x)['Weight'],
                                     lambda x: opt_prob.evaluate_constraint(x)['WeightGradient'],
                                     optimizer.minf*np.ones(1),
                                     400*np.ones(1),
                                     name='Weight')
    optimizer.add_nonlin_constraints(6,
                                     lambda x: opt_prob.evaluate_constraint(x)['Damping'],
                                     lambda x: opt_prob.evaluate_constraint(x)['DampingGradient'],
                                     4.2e-4*np.ones(6),
                                     optimizer.inf*np.ones(6),
                                     name='Damping')
    return (opt_prob, optimizer)


def run_optimizer(optimizer):
    """Run the optimizer."""
    optimizer.optimize('pyoptsparse:SNOPT')

def objective_function_only():
    """Evaluate only the objective function at a design point."""
    design_point_string_list = sys.argv[2:8]
    design_point_float_list = [float(x) for x in design_point_string_list]
    opt_prob = ARW2Problem(load_matlab=False)
    opt_prob.evaluate_objective(design_point_float_list)

def no_arguments():
    """If there are no command-line arguments, simply run the optimization problem."""
    (opt_prob, optimizer) = setup_optimizer()
    run_optimizer(optimizer)

if __name__ == '__main__':
    if len(sys.argv) >= 8:
        if sys.argv[1] == 'objective':
            objective_function_only()
    else:
        no_arguments()
