import pyoptsparse
def objfunc(xdict):
    x = xdict['xvars']
    funcs = {}
    funcs['obj'] = -x[0]*x[1]*x[2]
    conval = [0]*2
    conval[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
    conval[1] = -x[0] - 2.*x[1] - 2.*x[2]
    funcs['con'] = conval
    fail = False

    return funcs, fail

optProb = pyoptsparse.Optimization('TP037', objfunc)
optProb.addVarGroup('xvars',3, 'c',lower=[0,0,0], upper=[42,42,42], value=10)
optProb.addConGroup('con',2, lower=None, upper=0.0)
optProb.addObj('obj')
print(optProb)
opt = pyoptsparse.SNOPT()
sol = opt(optProb, sens='FD')
print(sol)