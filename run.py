#!/usr/bin/env python


def run_problem(solver, problem, domain, boundary, name=None, quiet=True, CFL=0.3, tfinal=0.2):
    from time import time
    if name is None: name = problem.__class__

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
        solver.advance_state(P,dt)
        t += dt
        nc += 1
        step_time = time()-start

        msg_data = (nc, t, dt, 1e6*step_time/P.size, 0)
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
        from os import listdir, system
        from pickle import load
        from numpy import zeros

        domain.dump(P)

        fnames = [x for x in listdir('.') if x.endswith('.pk')]
        tiles = [load(open(fn)) for fn in fnames]
        ndims = len(tiles[0]['data'].shape)-1
        Nq = tiles[0]['data'].shape[-1]
        N = tuple(domain.global_shape)
        P = zeros(N+(Nq,))

        for i,t in enumerate(tiles):
            if len(N) == 1:
                pi = t['data']
                i0 = t['global_start'][0]
                i1 = t['global_start'][0] + pi.shape[0]
                P[i0:i1,:] = pi

            if len(N) == 2:
                pi = t['data']
                i0,j0 = t['global_start']
                i1,j1 = [g+s for g,s in zip(t['global_start'], pi.shape)]
                P[i0:i1,j0:j1,:] = pi

        return P

    else:
        return P


def eulers():
    from hydro import visual, EulersEquationsSolver
    from hydro.testbench import SRShockTube1, RMHDCylindricalA
    from hydro.boundary import OutflowBoundary
    from hydro.domain import SimpleCartesianDomain
    from hydro.parallel import DistributedDomain

    DomainClass = DistributedDomain

    solver = EulersEquationsSolver(scheme='midpoint')
    problem = RMHDCylindricalA()
    domain = DomainClass((1024,1024), (-1.0,-1.0), (1.0,1.0))
    boundary = OutflowBoundary()

    P = run_problem(solver, problem, domain, boundary,
                    quiet=(domain.rank is not 0), CFL=0.4, tfinal=0.1)

    if domain.rank is 0:
        visual.hyd_four_pane_2d(P)
        visual.show()


if __name__ == "__main__":
    
    eulers()
