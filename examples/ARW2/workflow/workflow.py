"""Run the optimization problem."""

import warnings
import pdb

import numpy as np

from pyaeroopt.test.ARW2.pyao.ARW2opt import ARW2Problem
from pyaeroopt.opt import Optimize as Optimize

def run_optimizer():
    """Run the optimization problem."""

    # Define the optimization problem
    opt_prob = ARW2Problem()

    # Set up the optimizer by passing in methods defined in the opt_prob object
    optimizer = Optimize(name='ARW2')
    optimizer.add_variables(6, [0, 0, 0, 0, 0, 0],
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

    optimizer.optimize('pyoptsparse:SLSQP')
    #optimizer.optimize('pyoptsparse:SNOPT')

if __name__ == '__main__':
    warnings.simplefilter(action='ignore', category=FutureWarning)
    run_optimizer()
