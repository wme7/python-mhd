#!/usr/bin/env python


def run_problem(solver, problem, domain, boundary, name=None, quiet=True, CFL=0.3, tfinal=0.2):
    from time import time
    if name is None: name = problem.__class__
    from numpy import array, float64

    Nq = solver.NumComponents
    Ng = solver.NumGhostCells
    P  = problem.initial_model(domain, Ng, Nq)

    t  = 0.0
    nc = 0
    dt = CFL * min(domain.dx)
    start_time = time()

    solver.new_problem(domain)

    while t < tfinal:

        start = time()

        domain.set_BC(P, Ng, boundary)
        failures = solver.advance_state(P,dt)

        if failures and False:
            from pylab import imshow, show
            print failures
            imshow(solver.get_failure_mask(), interpolation='nearest')
            show()
            exit(1)

        t += dt
        nc += 1
        step_time = time()-start

        msg_data = (nc, t, dt, 1e6*step_time/P.size, failures)
        msg_text = "N: %05d t: %6.4f dt: %6.4e us/zone: %5.4f failures: %d"
        if not quiet: print msg_text % msg_data

    if not quiet:
        print "Python driver finished '%s'... total time: %f" % (name, time() - start_time)

    if len(P.shape) == 2:
        P = P[Ng:-Ng]
    if len(P.shape) == 3:
        P = P[Ng:-Ng,Ng:-Ng]
    if len(P.shape) == 4:
        P = P[Ng:-Ng,Ng:-Ng,Ng:-Ng]

    if domain.is_distributed:
        domain.dump(P)
        return domain.read()
    else:
        return P


if __name__ == "__main__":
    from hydro import visual, EulersEquationsSolver, RMHDEquationsSolver, SRHDEquationsSolver
    from hydro.testbench import *
    from hydro.boundary import OutflowBoundary
    from hydro.domain import SimpleCartesianDomain
    from hydro.parallel import DistributedDomain

    DomainClass = DistributedDomain

    solver = EulersEquationsSolver(scheme='midpoint')
    problem = RMHDCylindricalA()
    domain = DomainClass((32,32), (-0.5,-0.5), (0.5,0.5))
    boundary = OutflowBoundary()

    P = run_problem(solver, problem, domain, boundary,
                    quiet=(domain.rank is not 0), CFL=0.4, tfinal=0.1)

    if domain.rank is 0:
        visual.four_pane_2d(P)
        visual.show()
