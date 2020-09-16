#!/usr/bin/env python
"""
Optimization Interface
"""

## @package opt
#  Optimization Interface
#
#  Class for interfacing to various python optimization solvers (pyOpt,
#  scipy.optimize, nlopt currently supported)

import numpy as np
# import pyOpt         #not supported currently
import scipy.optimize
# import pyoptsparse   # removed until it can be integrated more smoothly
# import openopt
# import nlopt


class Optimize(object):
    ## Initialize optimization class with empty/default properties.  Objective
    #  function and constraints can be added via calls to other methods
    #  addObjective, ...
    def __init__(self, **kwargs):

        # Initialize relevant counters, opt variable itertes, \pm \infty
        self.counters = {}
        self.var_iter = []
        self.inf = 1.0e20 if 'inf' not in kwargs else kwargs['inf']
        self.minf = -1.0e20 if 'minf' not in kwargs else kwargs['minf']

        # Initialize variable-related properties
        self.nvar = 0
        self.var_init = np.zeros(self.nvar, dtype='float')
        self.var_bounds = np.zeros((self.nvar, 2), dtype='float')

        # Initialize objective-related properties
        self.objective = None
        self.objective_name = None
        self.gradient = None

        # Initialize linear-constraint-related properties
        self.n_lin_constr = 0
        self.lin_mat = np.zeros((self.n_lin_constr, self.nvar), dtype='float')
        self.lin_bounds = np.zeros((self.n_lin_constr, 2), dtype='float')

        # Initialize nonlinear-constraint-related properties
        self.n_nonlin_constr = 0
        self.nonlin_constr = []
        self.nonlin_constr_fail = []   # note: any([]) = False
        self.nonlin_jac = []
        self.nonlin_jac_fail = []
        self.nonlin_bounds = np.zeros((self.n_nonlin_constr, 2), dtype='float')
        self.n_nonlin_constr_per_group = np.zeros(1, dtype='int')
        self.nonlin_bounds_per_group = dict()
        self.nonlin_constraint_names = dict()

        # Give the problem a name
        self.name = 'Opt. Problem' if 'name' not in kwargs else kwargs['name']

    def add_variables(self, nvar, init, low, up):
        """Add variables and bounds to the Optimize object.

        Parameters:
        ------------
        nvar: int
            the number of variables to add
        init: 1D array of floats
            a list of initial values for the variables
        low: 1D array of floats
            a list of lower bounds for the variables
        up: 1D array of floats
            a list of upper bounds for the variables
        """
        self.nvar += nvar
        self.var_init = np.r_[self.var_init, init]
        self.var_bounds = np.r_[self.var_bounds, np.c_[low, up]]
        # Resize linear matrix
        self.lin_mat.resize((self.n_lin_constr, self.nvar))

    def add_objective(self, obj, name=None, grad=None, hess=None,
                      hess_vec=None, obj_failure=lambda x: False,
                      grad_failure=lambda x: False):
        """Add an objective function, and potentially also its derivatives.

        Parameters:
        ------------
        obj: Python function
            The objective function to be optimized.  Must return a scalar.
        name: string
            the name of the objective function.
        grad: Python function
            The gradient of the objective. Must return a 1D array with
            len(grad_return)=nvars, the number of design variables
        hess: Python function
            The Hessian of the objective function.  Must return a 2D array with
            hess_return.shape=(nvar, nvar)
        hess_vec: Python function
            a function that returns the product of the Hessian and a vector x,
            as a 1D array with len(hess_vec_return)=nvar
        obj_failure: Python function
            a function that returns True if the objective evaluation failed,
            and False otherwise.
        grad_failure: Python function
            a function that returns True if the gradient evaluation failed,
            and false otherwise
        """
        self.objective = obj
        self.objective_name = name if name is not None else 'obj'
        self.obj_failure = obj_failure
        self.gradient = grad
        self.grad_failure = grad_failure
        self.hessian = hess
        self.hessian_vec = hess_vec

    def add_lin_constraints(self, nconstr, mat, low, up):
        """Add linear constraints of the form low <= mat*x <= up.

        Parameters:
        ---------------
        nconstr: int
            the number of variables in the constraint
        mat: 2D array
            the matrix A in the constraint low <= A*x <= up
        low: 1D array of floats
            the lower bound for the linear constraint, with len(low)=nconstr
        up: 1D array of floats
            the upper bound for the linear constraint, with len(up)=nconstr
        """
        self.n_lin_constr += nconstr
        self.lin_mat = np.r_[self.lin_mat, mat]
        self.lin_bounds = np.r_[self.lin_bounds, np.c_[low, up]]

    def add_nonlin_constraints(self, nconstr, constr, jac, low, up, name=None,
                               constr_fail=lambda x: False,
                               jac_fail=lambda x: False):
        """Add nonlinear constraints of the form low <= constr <= up.

        Parameters:
        --------------
        nconstr: int
            the number of variables in the constraint
        constr: Python function
            a function returning a 1D array with len(constr_return)=nconstr,
            which specifies a nonlinear constraint
        jac: Python function
            a callable function returning the derivative of constr as a 2D array
            with jac.shape=(nconstr, nvar)
        low: 1D array of floats
            the lower bound for the constraint, with len(low)=nconstr
        up: 1D array of floats
            the upper bound for the constraint, with len(up)=nconstr
        name: string
            the name of this constraint group
        constr_fail: Python function
            a function returning True if the constraint evaluation failed,
            and False otherwise
        jac_fail: Python function
            a function returning True if the constraint jacobian evalution
            failed, and False otherwise.
        """
        self.n_nonlin_constr += nconstr
        self.nonlin_constr.append(constr)
        self.nonlin_jac.append(jac)
        self.nonlin_bounds = np.r_[self.nonlin_bounds, np.c_[low, up]]
        next_ind = len(self.n_nonlin_constr_per_group) - 1
        self.nonlin_bounds_per_group[next_ind] = np.r_[np.c_[low, up]]
        if name is not None:
            self.nonlin_constraint_names[next_ind] = name
        self.n_nonlin_constr_per_group = np.r_[self.n_nonlin_constr_per_group,
                                               np.array([nconstr])]
        self.nonlin_constr_fail.append(constr_fail)
        self.nonlin_jac_fail.append(jac_fail)

    def eval_nonlin_constraints(self, param, which, ind=None):
        """Evaluate a single nonlinear constraint equation."""
        # Determine constraints to evaluate
        if ind is None: ind = np.arange(self.n_nonlin_constr)
        if type(ind) is int: ind = np.array([ind])

        # Initialize arrays for constraints and/or jacobians
        if which in ['constr', 'both']:
            constrs = np.zeros(len(ind), dtype=float)
        if which in ['jac', 'both']:
            jacs = np.zeros((len(ind), self.nvar), dtype=float)

        ind_per_group = []
        # Determine indices per group
        for g in range(len(self.nonlin_constr)):
            lb = sum(self.n_nonlin_constr_per_group[:g+1])
            ub = sum(self.n_nonlin_constr_per_group[:g+2])
            ind_below = ind[np.logical_and(ind < ub, ind >= lb)]
            ind_per_group.append(ind_below)

        # Evaluate appropriate constraints and/or jacobians
        # cnt = 0
        # for g in range(len(self.nonlin_constr)):
        #     if len(ind_per_group[g]) > 0:
        #         idx = ind_per_group[g]
        #         if which in ['constr','both']:
        #             constrs[cnt:cnt+len(idx)] = self.nonlin_constr[g](param)[idx]
        #         if which in ['jac','both']:
        #             jacs[cnt:cnt+len(idx)]    = self.nonlin_jac[g](param)[idx, :]
        #         cnt += len(idx)
        for cnt, i in enumerate(ind):
            group = (self.n_nonlin_constr_per_group.cumsum() < i+1).nonzero()[0][-1]
            index = i-self.n_nonlin_constr_per_group[:group+1].sum()
            if which in ['constr', 'both']:
                constrs[cnt] = self.nonlin_constr[group](param)[index]
            if which in ['jac', 'both']:
                jacs[cnt] = self.nonlin_jac[group](param)[index]

        if which is 'constr': return constrs
        if which is 'jac': return jacs
        if which is 'both': return (constrs, jacs)

    def optimize(self, solver, sens=None, options=None, callback=None):
        """Solve the problem with the requested optimizer.

        Parameters:
        -------------
        solver: string
            a string specifying the solver to use, of the form 'package:solver'
        sens: string or None
            a string specifying the type of gradient approximation to use.
            'FD' for forward difference, 'CD' for central difference,
            'FDR' for forward difference with relative step size,
            'CDR' for central difference with relative step size,
             and 'CS' for complex step. None indicates that either gradients
             have been provided or the default finite differencing scheme
             for your chosen optimizer should be used.
        options: dictionary
            a dictionary of options to pass to the solver
        callback: Python function
            for scipy, a function that is called after each iteration as
            callback(xk), where xk is the parameter vector after iteration k.

        Returns:
        ----------
        xStar: 1D array
            the optimized parameter vector
        fStar: float
            The optimized objective value
        """
        if 'pyoptsparse' in solver:
            xStar, fStar = self.optimize_pyoptsparse(solver, sens, options,
                                                     callback)
        elif 'scipy' in solver:
            xStar, fStar = self.optimize_scipy(solver, sens, options,
                                               callback)
        elif 'nlopt' in solver:
            xStar, fStar = self.optimize_nlopt(solver, sens, options)
        elif 'openopt' in solver:
            xStar, fStar = self.optimize_openopt(solver, sens, options,
                                                 callback)
        elif 'pyopt' in solver:
            raise ValueError('pyOpt not supported in pyaeroopt34')
        return (xStar, fStar)

    def optimize_pyoptsparse(self, solver, sens, options, callback=None):
        """Solve the problem with pyoptsparse."""

        def objfunc_builder():
            """Builds an objective function with a name chosen at run-time."""
            def objfunc(xdictInside):
                """Give the objective and constraints in pyoptsparse format."""
                funcs = dict()
                x = xdictInside['vars']  # extract the variables as a 1D array
                funcs[self.objective_name] = self.objective(x)

                # add the constraints to the funcs dictionary, one key per group
                n_nonlin_constr_per_group = self.n_nonlin_constr_per_group[1:]
                for ind, numConstr in enumerate(n_nonlin_constr_per_group):
                    if ind in self.nonlin_constraint_names:
                        key = self.nonlin_constraint_names[ind]
                    else:
                        key = 'constraint_group{0:d}'.format(ind)
                    funcs[key] = np.array(self.nonlin_constr[ind](x), ndmin=1)
                obj_fail = self.obj_failure(x)
                constr_fail_list = [f(x) for f in self.nonlin_constr_fail]
                fail = any([constr_fail_list, obj_fail])
                return funcs, fail
            objfunc.__name__ = self.objective_name
            return objfunc

        if sens in (None, 'provided'):
            def sens(xdict, funcs):
                """Specify objective/constraint gradients in pyoptparse format."""
                func_sens = dict()
                x = xdict['vars'] # extract the design variables as a 1D array
                func_sens[self.objective_name, 'vars'] = np.array(self.gradient(x),
                                                                  ndmin=1)

                # add the constraints to the funcs dictionary, one key per group
                n_nonlin_constr_per_group = self.n_nonlin_constr_per_group[1:]
                for ind, numConstr in enumerate(n_nonlin_constr_per_group):
                    if ind in self.nonlin_constraint_names:
                        key = self.nonlin_constraint_names[ind]
                    else:
                        key = 'constraint_group{0:d}'.format(ind)
                    func_sens[key, 'vars'] = np.array(self.nonlin_jac[ind](x),
                                                     ndmin=2)
                obj_grad_fail = self.grad_failure(x)
                constr_jac_fail_list = [f(x) for f in self.nonlin_jac_fail]
                fail = any([constr_jac_fail_list, obj_grad_fail])
                return func_sens, fail

        # Warn the user if they tried to use a callback function
        if callback is not None:
            raise ValueError('pyoptsparse does not support callbacks.')

        # Create an Optimization instance and add the objective function
        objfunc = objfunc_builder()
        print(self.name, objfunc)
        opt_prob = pyoptsparse.Optimization(self.name, objfunc)
        opt_prob.addObj(self.objective_name)

        # Add the design variables
        opt_prob.addVarGroup('vars', self.nvar, lower=self.var_bounds[:, 0],
                             upper=self.var_bounds[:, 1], value=self.var_init)

        # Add the constraints
        # Add the nonlinear constraints
        n_nonlin_constr_per_group = self.n_nonlin_constr_per_group[1:]
        for ind, numConstr in enumerate(n_nonlin_constr_per_group):
            if ind in self.nonlin_constraint_names:
                constraint_name = self.nonlin_constraint_names[ind]
            else:
                constraint_name = 'constraint_group{0:d}'.format(ind)
            opt_prob.addConGroup(constraint_name,
                                 numConstr,
                                 lower=self.nonlin_bounds_per_group[ind][:, 0],
                                 upper=self.nonlin_bounds_per_group[ind][:, 1])

        # Add the linear constraints
        opt_prob.addConGroup('Linear Constraints',
                             self.n_lin_constr,
                             lower=self.lin_bounds[:, 0],
                             upper=self.lin_bounds[:, 1],
                             linear=True,
                             wrt='vars',
                             jac={'vars': self.lin_mat})

        # Solve the optimization problem
        print(opt_prob)
        opt = eval('pyoptsparse.'+ solver.lstrip('pyoptsparse:')+'()')
        sol = opt(opt_prob, sens=sens)
        print(sol)

        # Extract the optimized parameter vector and objective value
        xStar = [None]*self.nvar
        for i in range(self.nvar):
            xStar[i] = sol.variables['vars'][i].value
        fStar = sol.objectives[self.objective_name].value

        return (xStar, fStar)

    def optimize_scipy(self, solver, sens, options, callback=None):
        """Solve the problem with scipy."""

        # Reformat variable bounds
        bnds = [(x[0] if x[0] > self.minf else None,
                 x[1] if x[1] < self.inf  else None) for x in self.var_bounds]

        # Need to define alternate functions to handle scoping issues with
        # lambdas (specifically, when trying to use index in a loop to
        # extract specific constraints to put into scipy format).  This will
        # essentially define new scopes.
        def bndValFactory(k, which):
            if which == 'low':
                return (lambda x: x[k] - self.var_bounds[k, 0])
            if which == 'up':
                return (lambda x: self.var_bounds[k, 1] - x[k])
        def bndGradFactory(k, which):
            ek = np.zeros(self.nvar, dtype=float)
            ek[k] = 1.0
            if which == 'low':
                return lambda x: ek
            if which == 'up':
                return lambda x: -ek
        def linValFactory(k, which):
            if which == 'low':
                return(lambda x: np.dot(self.lin_mat[k, :], x) -
                       self.lin_bounds[k, 0])
            if which == 'up':
                return(lambda x: self.lin_bounds[k, 1] -
                       np.dot(self.lin_mat[k, :], x))
        def linGradFactory(k, which):
            if which == 'low':
                return (lambda x: self.lin_mat[k, :])
            if which == 'up':
                return (lambda x: -self.lin_mat[k, :])
        def nonlinValFactory(k, which):
            if which == 'low':
                return(lambda x: self.eval_nonlin_constraints(x, 'constr', k) -
                       self.nonlin_bounds[k, 0])
            if which == 'up':
                return(lambda x: self.nonlin_bounds[k, 1] -
                       self.eval_nonlin_constraints(x, 'constr', k))
        def nonlinGradFactory(k, which):
            if which == 'low':
                return lambda x: self.eval_nonlin_constraints(x, 'jac', k)
            if which == 'up':
                return lambda x: -self.eval_nonlin_constraints(x, 'jac', k)

        # Reformat constraints (don't distinguish between lin/nonlin)
        constrs = []
        # Add linear equality constraints
        if self.n_lin_constr > 0:
            for k, bnd in enumerate(self.lin_bounds):
                if bnd[0] == bnd[1]:
                    constrs.append({'type':'eq',
                                    'fun': linValFactory(k, 'low'),
                                    'jac': linGradFactory(k, 'low')})

        # Add nonlinear equality constraints
        if self.n_nonlin_constr > 0:
            for k, bnd in enumerate(self.nonlin_bounds):
                if bnd[0] == bnd[1]:
                    if None in self.nonlin_jac:
                        constrs.append({'type':'eq',
                                        'fun': nonlinValFactory(k, 'low')})
                    else:
                        constrs.append({'type':'eq',
                                        'fun': nonlinValFactory(k, 'low'),
                                        'jac': nonlinGradFactory(k, 'low')})

        # Add linear inequality constraints
        if self.n_lin_constr > 0:
            for k, bnd in enumerate(self.lin_bounds):
                if bnd[0] < bnd[1]:
                    if bnd[0] > self.minf:
                        constrs.append({'type':'ineq',
                                        'fun': linValFactory(k, 'low'),
                                        'jac': linGradFactory(k, 'low')})
                    if bnd[1] < self.inf:
                        constrs.append({'type':'ineq',
                                        'fun': linValFactory(k, 'up'),
                                        'jac': linGradFactory(k, 'up')})

        # Add nonlinear inequality constraints
        if self.n_nonlin_constr > 0:
            for k, bnd in enumerate(self.nonlin_bounds):
                if bnd[0] < bnd[1]:
                    if bnd[0] > self.minf:
                        if None in self.nonlin_jac:
                            constrs.append({'type':'ineq',
                                            'fun': nonlinValFactory(k, 'low')})
                        else:
                            constrs.append({'type':'ineq',
                                            'fun': nonlinValFactory(k, 'low'),
                                            'jac': nonlinGradFactory(k, 'low')})
                    if bnd[1] < self.inf:
                        if None in self.nonlin_jac:
                            constrs.append({'type':'ineq',
                                            'fun': nonlinValFactory(k, 'up')})
                        else:
                            constrs.append({'type':'ineq',
                                            'fun': nonlinValFactory(k, 'up'),
                                            'jac': nonlinGradFactory(k, 'up')})

        # If using COBYLA, add simple bounds as general constraints and
        # treat all equality constraints as two inequality constraints
        if 'COBYLA' in solver:
            constrs = []

            # Add linear constraints
            if self.n_lin_constr > 0:
                for k, bnd in enumerate(self.lin_bounds):
                    if bnd[0] > self.minf:
                        constrs.append({'type':'ineq',
                                        'fun': linValFactory(k, 'low'),
                                        'jac': linGradFactory(k, 'low')})
                    if bnd[1] < self.inf:
                        constrs.append({'type':'ineq',
                                        'fun': linValFactory(k, 'up'),
                                        'jac': linGradFactory(k, 'up')})

            # Add nonlinear constraints
            if self.n_nonlin_constr > 0:
                for k, bnd in enumerate(self.nonlin_bounds):
                    if bnd[0] > self.minf:
                        constrs.append({'type':'ineq',
                                        'fun': nonlinValFactory(k, 'low'),
                                        'jac': nonlinGradFactory(k, 'low')})
                    if bnd[1] < self.inf:
                        constrs.append({'type':'ineq',
                                        'fun': nonlinValFactory(k, 'up'),
                                        'jac': nonlinGradFactory(k, 'up')})

            # Add bound constraints
            for k, bnd in enumerate(self.var_bounds):
                ek = np.zeros(self.nvar, dtype=float)
                ek[k] = 1.0
                if bnd[0] > self.minf:
                    constrs.append({'type':'ineq',
                                    'fun': bndValFactory(k, 'low'),
                                    'jac': bndGradFactory(k, 'low')})
                if bnd[1] < self.inf:
                    constrs.append({'type':'ineq',
                                    'fun': bndValFactory(k, 'up'),
                                    'jac': bndGradFactory(k, 'up')})

        summ = scipy.optimize.minimize(self.objective, self.var_init,
                                       method=solver.lstrip('scipy:'),
                                       bounds=bnds, jac=self.gradient,
                                       constraints=constrs, hess=self.hessian,
                                       hessp=self.hessian_vec,
                                       callback=callback, options=options)
        print(summ)
        return (summ.x, summ.fun)



    #TODO: test
    def optimize_pyopt(self, solver, sens, options):
        """Deprecated -- only supported in Python 2.7"""

        # pyOpt requires empty options to be specified as {}, not None
        if options is None: options = {}

        eqConstr = ( 'SLSQP'    not in solver and
                     'CONMIN'   not in solver and
                     'COBYLA'   not in solver and
                     'FILTERSD' not in solver and
                     'SDPEN'    not in solver )

        # Organize constraints
        indLinEq = []; indLinIneq = []; indNonlinEq = []; indNonlinIneq = [];
        for k, bnd in enumerate(self.lin_bounds):
            if bnd[0]==bnd[1] and eqConstr: indLinEq.append(k)
            else:                             indLinIneq.append(k)
        indLinEq=np.array(indLinEq); indLinIneq=np.array(indLinIneq);
        for k, bnd in enumerate(self.nonlin_bounds):
            if bnd[0]==bnd[1] and eqConstr: indNonlinEq.append(k)
            else:                             indNonlinIneq.append(k)
        indNonlinEq=np.array(indNonlinEq);indNonlinIneq=np.array(indNonlinIneq);

        # pyOpt objective
        def objectivePyopt(xIn):

            x = xIn[:self.nvar]

            f = self.objective(x)
            g = np.zeros( 0, dtype = float)
            if len(indLinEq) > 0:
                g = np.r_[ g, np.dot(self.lin_mat[indLinEq,:],x)]
            if len(indNonlinEq) > 0:
                g =np.r_[g, self.eval_nonlin_constraints(x,'constr',indNonlinEq)]
            if len(indLinIneq) > 0:
                g = np.r_[ g, np.dot(self.lin_mat[indLinIneq,:],x)]
            if len(indNonlinIneq) > 0:
                g=np.r_[g, self.eval_nonlin_constraints(x,'constr',indNonlinIneq)]

            fail = 0
            if hasattr(f, '__iter__'):
                f = f[0]
            if f >= self.inf or np.any( g >= self.inf): fail = 1

            return f, g, fail

        # pyOpt gradient
        def gradientPyopt(xIn,f,g):

            x = xIn[:self.nvar]

            df = self.gradient(x)
            dg = np.zeros( (0,self.nvar), dtype = float)
            if len(indLinEq) > 0:
                dg = np.r_[ dg, self.lin_mat[indLinEq,:] ]
            if len(indNonlinEq) > 0:
                dg= np.r_[ dg, self.eval_nonlin_constraints(x,'jac',indNonlinEq) ]
            if len(indLinIneq) > 0:
                dg = np.r_[ dg, self.lin_mat[indLinIneq,:] ]
            if len(indNonlinIneq) > 0:
                dg=np.r_[dg, self.eval_nonlin_constraints(x,'jac',indNonlinIneq)]

            fail = 0
            if hasattr(f, '__iter__'):
                f = f[0]
            if f >= self.inf or np.any( g >= self.inf): fail = 1

            return df.reshape(1, -1), dg, fail

        # Instantiate optimization problem
        optProb = pyOpt.Optimization('pyopt', objectivePyopt)

        # Add objective
        optProb.addObj('objective')

        # Add variables
        optProb.addVarGroup('var',self.nvar,type='c',value=self.var_init,
                            lower=self.var_bounds[:,0],upper=self.var_bounds[:,1])
        # Add constraints
        if len(indLinEq) > 0:
            optProb.addConGroup('lin-equality',len(indLinEq),type='e',
                                               equal=self.lin_bounds[indLinEq,0])
        if len(indNonlinEq) > 0:
            optProb.addConGroup('nonlin-equality',len(indNonlinEq),type='e',
                                         equal=self.nonlin_bounds[indNonlinEq,0])
        if len(indLinIneq) > 0:
            optProb.addConGroup('lin-inequality',len(indLinIneq),type='i',
                                             lower=self.lin_bounds[indLinIneq,0],
                                             upper=self.lin_bounds[indLinIneq,1])
        if len(indNonlinIneq) > 0:
            optProb.addConGroup('nonlin-inequality',len(indNonlinIneq),type='i',
                                       lower=self.nonlin_bounds[indNonlinIneq,0],
                                       upper=self.nonlin_bounds[indNonlinIneq,1])

        # Setup solver
        if 'SNOPT' in solver:
            optimizer = pyOpt.pySNOPT.SNOPT(options=options)
        if 'SLSQP' in solver:
            optimizer = pyOpt.SLSQP(options=options)
        if 'CONMIN' in solver:
            optimizer = pyOpt.CONMIN(options=options)
        if 'ALGENCAN' in solver:
            optimizer = pyOpt.ALGENCAN(options=options)
        if 'ALPSO' in solver:
            optimizer = pyOpt.ALPSO(options=options)
        if 'ALHSO' in solver:
            optimizer = pyOpt.ALHSO(options=options)
        if 'COBYLA' in solver:
            optimizer = pyOpt.COBYLA(options=options)
        if 'FILTERSD' in solver:
            optimizer = pyOpt.FILTERSD(options=options)
        if 'KOPT' in solver:
            optimizer = pyOpt.KOPT(options=options)
        if 'MIDACO' in solver:
            optimizer = pyOpt.MIDACO(options=options)
        if 'KSQP' in solver:
            optimizer = pyOpt.KSQP(options=options)
        if 'SDPEN' in solver:
            optimizer = pyOpt.SDPEN(options=options)
        if 'SOLVOPT' in solver:
            optimizer = pyOpt.SOLVOPT(options=options)

        # Run optimization
        if sens == 'finite-difference':
            optimizer(optProb,sens_type='FD')
        else:
            optimizer(optProb,sens_type=gradientPyopt)

        # Extract solution
        j = len(optProb._solutions)-1
        xStar = np.zeros(self.nvar)
        for k in range(self.nvar):
            xStar[k] = optProb._solutions[j].getVar(k).value

        pdb.set_trace()
        return xStar,objectivePyopt(xStar)[0]

    #TODO: test
    def optimize_nlopt(self,solver,sens,options):
        string = 'nlopt.' + solver.lstrip('nlopt:')
        # optProb = nlopt.opt( getattr( nlopt , solver.lstrip('nlopt:') ),
                             # self.nvar)
        optProb = nlopt.opt( nlopt.LD_SLSQP ,self.nvar)

        def objectiveNlopt(x,grad):

            if grad.size > 0:
                grad[:] = self.gradient( x )

            return ( self.objective( x ) )

        # Organize constraints
        eqConstr = ( 'MMA' not in solver and
                     'COBYLA' not in solver )

        indLinEq = []; indLinIneqUp = []; indLinIneqLow = [];
        indNonlinEq = []; indNonlinIneqUp = []; indNonlinIneqLow = [];
        for k, bnd in enumerate(self.lin_bounds):
            if bnd[0] == bnd[1] and eqConstr: indLinEq.append(k)
            else:
                if bnd[0] > self.minf : indLinIneqLow.append(k)
                if bnd[1] < self.inf  : indLinIneqUp.append(k)
        indLinEq=np.array(indLinEq)
        indLinIneqUp=np.array(indLinIneqUp)
        indLinIneqLow=np.array(indLinIneqLow)

        for k, bnd in enumerate(self.nonlin_bounds):
            if bnd[0] == bnd[1] and eqConstr: indNonlinEq.append(k)
            else:
                if bnd[0] > self.minf : indNonlinIneqLow.append(k)
                if bnd[1] < self.inf  : indNonlinIneqUp.append(k)
        indNonlinEq=np.array(indNonlinEq)
        indNonlinIneqUp=np.array(indNonlinIneqUp)
        indNonlinIneqLow=np.array(indNonlinIneqLow)

        # linear constraints
        linEqConstr = np.zeros( (0,self.nvar) , dtype=float )
        linEqRhs    = np.zeros( 0 , dtype=float )
        if len(indLinEq) > 0:
            linEqConstr = self.lin_mat[indLinEq,:]
            linEqRhs    = self.lin_bounds[indLinEq,0]

        linIneqConstr = np.zeros( (0,self.nvar) , dtype=float )
        linIneqRhs    = np.zeros( 0 , dtype=float )
        if len(indLinIneqLow) > 0:
            linIneqConstr = np.r_[ linIneqConstr,
                                   -self.lin_mat[indLinIneqLow,:] ]
            linIneqRhs    = np.r_[ linIneqRhs,
                                   -self.lin_bounds[indLinIneqLow,0] ]
        if len(indLinIneqUp) > 0:
            linIneqConstr = np.r_[ linIneqConstr,
                                   self.lin_mat[indLinIneqUp,:] ]
            linIneqRhs    = np.r_[ linIneqRhs,
                                   self.lin_bounds[indLinIneqUp,1] ]

        # constraints
        def linEqConstrNlopt(results,x,grad):

            results[:] = np.dot(linEqConstr,x) - linEqRhs
            if grad.size > 0:
                grad[:] = linEqConstr

        def linIneqConstrNlopt(results,x,grad):

            results[:] = np.dot(linIneqConstr,x) - linIneqRhs
            if grad.size > 0:
                grad[:] = linIneqConstr

        def nonlinEqConstrNlopt(results,x,grad):

            results[:] = (self.eval_nonlin_constraints(x,'constr',indNonlinEq) -
                                               self.nonlin_bounds[indNonlinEq,0])
            if grad.size > 0:
                grad[:] = self.eval_nonlin_constraints(x,'jac',indNonlinEq)

        def nonlinIneqConstrNlopt(results,x,grad):

            tmp = np.zeros( 0, dtype=float )
            if len(indNonlinIneqLow) > 0:
                tmp = np.r_[ tmp,
                  -self.eval_nonlin_constraints(x,'constr',indNonlinIneqLow) +
                                          self.nonlin_bounds[indNonlinIneqLow,0]]
            if len(indNonlinIneqUp) > 0:
                tmp = np.r_[ tmp,
                   self.eval_nonlin_constraints(x,'constr',indNonlinIneqUp)  -
                                          self.nonlin_bounds[indNonlinIneqUp,1]]
            results[:] = tmp

            if grad.size > 0:
                tmpGrad = np.zeros( (0,self.nvar), dtype=float)
                if len(indNonlinIneqLow) > 0:
                    tmpGrad = np.r_[ tmpGrad,
                          -self.eval_nonlin_constraints(x,'jac',indNonlinIneqLow)]
                if len(indNonlinIneqUp ) > 0:
                    tmpGrad = np.r_[ tmpGrad,
                           self.eval_nonlin_constraints(x,'jac',indNonlinIneqUp)]
                grad[:] = tmpGrad

        optProb.set_min_objective( objectiveNlopt )
        optProb.set_lower_bounds(self.var_bounds[:,0])
        optProb.set_upper_bounds(self.var_bounds[:,1])

        if linEqConstr.shape[0] > 0:
            optProb.add_equality_mconstraint( linEqConstrNlopt ,
                              1.0e-6*np.ones(linEqConstr.shape[0],dtype=float) )
        if indNonlinEq.shape[0] > 0:
            optProb.add_equality_mconstraint( nonlinEqConstrNlopt ,
                              1.0e-6*np.ones(indNonlinEq.shape[0],dtype=float) )
        if linIneqConstr.shape[0] > 0:
            optProb.add_inequality_mconstraint( linIneqConstrNlopt ,
                            1.0e-6*np.ones(linIneqConstr.shape[0],dtype=float) )
        if indNonlinIneqLow.shape[0] + indNonlinIneqUp.shape[0] > 0:
            optProb.add_inequality_mconstraint( nonlinIneqConstrNlopt ,
                                   1.0e-6*np.ones(indNonlinIneqLow.shape[0]+
                                                  indNonlinIneqUp.shape[0],
                                                  dtype=float) )

        optProb.set_xtol_abs(1.0e-8)
        xopt = optProb.optimize(self.var_init)

        return xopt, self.objective( xopt )

    #TODO: test
    def optimize_openopt(self,solver,sens,options,callback=None):

        # Organize constraints
        indLinEq = []; indLinIneqUp = []; indLinIneqLow = [];
        indNonlinEq = []; indNonlinIneqUp = []; indNonlinIneqLow = [];
        for k, bnd in enumerate(self.lin_bounds):
            if bnd[0] == bnd[1]   : indLinEq.append(k)
            if bnd[0] > self.minf : indLinIneqLow.append(k)
            if bnd[1] < self.inf  : indLinIneqUp.append(k)
        indLinEq=np.array(indLinEq)
        indLinIneqUp=np.array(indLinIneqUp)
        indLinIneqLow=np.array(indLinIneqLow)

        for k, bnd in enumerate(self.nonlin_bounds):
            if bnd[0] == bnd[1]   : indNonlinEq.append(k)
            if bnd[0] > self.minf : indNonlinIneqLow.append(k)
            if bnd[1] < self.inf  : indNonlinIneqUp.append(k)
        indNonlinEq=np.array(indNonlinEq)
        indNonlinIneqUp=np.array(indNonlinIneqUp)
        indNonlinIneqLow=np.array(indNonlinIneqLow)

        # linear constraints
        linEqConstr = np.zeros( (0,self.nvar), dtype=float )
        linEqRhs    = np.zeros( 0 , dtype=float )
        if len(indLinEq) > 0:
            linEqConstr = self.lin_mat[indLinEq,:],
            linEqRhs    = self.lin_bounds[indLinEq,0]

        linIneqConstr = np.zeros( (0,self.nvar), dtype=float )
        linIneqRhs    = np.zeros( 0 , dtype=float )
        if len(indLinIneqLow) > 0:
            linIneqConstr = np.r_[ linIneqConstr,
                                   -self.lin_mat[indLinIneqLow,:] ]
            linIneqRhs    = np.r_[ linIneqRhs,
                                   -self.lin_bounds[indLinIneqLow,0] ]
        if len(indLinIneqUp) > 0:
            linIneqConstr = np.r_[ linIneqConstr,
                                   self.lin_mat[indLinIneqUp,:] ]
            linIneqRhs    = np.r_[ linIneqRhs,
                                   self.lin_bounds[indLinIneqUp,1] ]

        # nonlinear constraints
        def nonlinEqConstrOpenopt(x):

            if len(indNonlinEq) > 0:
                c = (self.eval_nonlin_constraints(x,'constr',indNonlinEq) -
                                               self.nonlin_bounds[indNonlinEq,0])
                return ( c )

        def nonlinEqJacOpenopt(x):

            if len(indNonlinEq) > 0:
                j = self.eval_nonlin_constraints(x,'jac',indNonlinEq)
                return ( j )

        def nonlinIneqConstrOpenopt(x):

            c = np.zeros(0,dtype=float)
            if len(indNonlinIneqLow) > 0:
                c = np.r_[ c,-self.eval_nonlin_constraints(x,'constr',
                       indNonlinIneqLow) +self.nonlin_bounds[indNonlinIneqLow,0]]
            if len(indNonlinIneqUp) > 0:
                c = np.r_[ c, self.eval_nonlin_constraints(x,'constr',
                        indNonlinIneqUp) - self.nonlin_bounds[indNonlinIneqUp,1]]
            return ( c )

        def nonlinIneqJacOpenopt(x):

            j = np.zeros((0,self.nvar),dtype=float)
            if len(indNonlinIneqLow) > 0:
                j = np.r_[j,
                    -self.eval_nonlin_constraints(x,'jac',indNonlinIneqLow) +
                                          self.nonlin_bounds[indNonlinIneqLow,0]]
            if len(indNonlinIneqUp) > 0:
                j = np.r_[j,
                   self.eval_nonlin_constraints(x,'jac',indNonlinIneqUp)  -
                                          self.nonlin_bounds[indNonlinIneqUp,1]]

            return ( j )

        # Turn unused contraints into None
        if linEqConstr.shape[0] == 0:
            linEqConstr = None
            linEqRhs    = None
        if linIneqConstr.shape[0] == 0:
            linIneqConstr = None
            linIneqRhs    = None
        if len(indNonlinEq) == 0:
            h  = None
            dh = None
        if len(indNonlinIneqUp) + len(indNonlinIneqLow) == 0:
            c  = None
            dc = None

        optProb = openopt.NLP(self.objective, self.var_init,
                              df=self.gradient,
                              lb=self.var_bounds[:,0],
                              ub=self.var_bounds[:,1],
                              Aeq=linEqConstr,beq=linEqRhs,
                              A=linIneqConstr,b=linIneqRhs,
                              h=nonlinEqConstrOpenopt,dh=nonlinEqJacOpenopt,
                              c=nonlinIneqConstrOpenopt,dc=nonlinIneqJacOpenopt)

        summ = optProb.solve(solver.lstrip('openopt:'))
        return summ.xf,summ.ff
