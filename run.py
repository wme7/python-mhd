#!/usr/bin/env python


def run_problem(solver, problem, name=None, quiet=True, CFL=0.1, tfinal=0.2):

    from time import time
    from numpy import zeros_like

    if name is None: name = problem.__class__

    P  = problem.initial_model(solver.N, solver.NumComponents)
    t  = 0.0
    dt = 1e-9
    nc = 0

    Ng = solver.get_Ng()
    Nx = P.shape[0]

    min_dx = min([L/N for L,N in zip(solver.L, solver.N)])

    start_time = time()
    while t < tfinal:

        nc += 1
        start = time()

        P[ 0:Ng   ] = P[   Ng  ] # Boundary conditions
        P[Nx-Ng:Nx] = P[Nx-Ng-1]

        solver.advance_state(P,dt)

        t += dt
        step_time = time()-start

        msg_data = (nc, t, dt, 1e6*step_time/P.size, 0)
        msg_text = "N: %05d t: %6.4f dt: %6.4e us/zone: %5.4f failures: %d"
        if not quiet: print msg_text % msg_data

        dt = CFL * min_dx



    print "Python driver finished '%s'... total time: %f" % (name, time() - start_time)
    return P



if __name__ == "__main__":

    from hydro import EulersEquationsSolver
    from numpy import array, zeros
    from hydro.testbench import *
    from hydro import visual

    solver = EulersEquationsSolver(N=[256], L=[1.0])
    problem = SRShockTube1()

    P = run_problem(solver, problem, quiet=False, CFL=0.1, tfinal=0.2)

    from pylab import plot, show, legend

    visual.hyd_shocktube(P)
    show()
