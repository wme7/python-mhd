#!/usr/bin/env python



def run_1d_problem(lib, state, problem, Nx=128, CFL=0.5, tfinal=0.2, verbose=False):

    from numpy import zeros

    P = zeros((Nx,8))
    U = zeros((Nx,8))
    L = zeros((Nx,8))

    state.adiabatic_gamma = problem.adiabatic_gamma
    problem.initial_model(P)
    lib.set_state(state)
    lib.initialize(P,Nx,1,1)
    lib.prim_to_cons_array(P,U)

    dx  = 1.0 / Nx
    t   = 0.0
    dt  = CFL * dx

    from time import time

    ttltime = 0.0
    n_cycle = 0
    while t < tfinal:

        start = time()

        e1 = lib.dUdt_1d(U,L)
        e2 = lib.dUdt_1d(U + 0.5*dt*L,L)

        if e1 or e2:
            print "Run crashed! Sorry...", e1, e2
            break

        U += dt*L
        t += dt
        n_cycle += 1

        U[ 0:4,   :] = U[   4,:] # Boundary conditions
        U[Nx-4:Nx,:] = U[Nx-5,:]

        ttltime += time()-start

        if verbose:
            print "t =", t

    lib.cons_to_prim_array(U,P)
    lib.finalize()
    print "Solver averaged %f us/zone" % (ttltime / (Nx*8*n_cycle)*1e6)
    return P




def sr_shocktube():

    from pylab import show
    import rmhd

    problem = SRShockTube2()
    state1  = rmhd.LibraryState(mode_reconstruct=1)
    state2  = rmhd.LibraryState(mode_reconstruct=2)

    P0 = rmhd.riemann.exact_sr_vt(problem,           tfinal=0.2)
    P1 = run_1d_problem(rmhd._lib2, state1, problem, tfinal=0.2, CFL=0.5)
    P2 = run_1d_problem(rmhd._lib2, state2, problem, tfinal=0.2, CFL=0.5)

    rmhd.visual.shocktube(P0, label="exact", linestyle='-', marker='None', lw=2)
    rmhd.visual.shocktube(P1, label="3vel", linestyle='--')
    rmhd.visual.shocktube(P2, label="4vel", linestyle='-.')
    show()




def compare_reconstruct():

    from pylab import show
    import rmhd

    problem = RMHDRotationalWave()

    state0  = rmhd.LibraryState(mode_reconstruct=0)
    state1  = rmhd.LibraryState(mode_reconstruct=1)
    state2  = rmhd.LibraryState(mode_reconstruct=2)

    P0 = run_1d_problem(rmhd._lib, state0, problem, CFL=0.2, tfinal=0.2)
    P1 = run_1d_problem(rmhd._lib, state1, problem, CFL=0.2, tfinal=0.2)
    P2 = run_1d_problem(rmhd._lib, state2, problem, CFL=0.2, tfinal=0.2)

    rmhd.visual.shocktube(P0, label="Piecewise Constant", linestyle='--', mfc='None')
    rmhd.visual.shocktube(P1, label="PLM 3-velocity", linestyle='-.', mfc='None')
    rmhd.visual.shocktube(P2, label="PLM 4-velocity", linestyle='-', marker='None')
    show()


def compare_quartic():

    from pylab import show
    import rmhd

    problem = RMHDShockTube4()

    state0  = rmhd.LibraryState(mode_quartic_solver=0)
    state1  = rmhd.LibraryState(mode_quartic_solver=1)
    state2  = rmhd.LibraryState(mode_quartic_solver=2)

    run_args = {'Nx':1024, 'CFL':0.5, 'tfinal':0.2}

    P0 = run_1d_problem(rmhd._lib, state0, problem, **run_args)
    P1 = run_1d_problem(rmhd._lib, state1, problem, **run_args)
    P2 = run_1d_problem(rmhd._lib, state2, problem, **run_args)

    rmhd.visual.shocktube(P0, label="Exact", linestyle='--', mfc='None')
    rmhd.visual.shocktube(P1, label="Approx1", linestyle='-.', mfc='None')
    rmhd.visual.shocktube(P2, label="Approx2", linestyle='-', marker='None')
    show()


if __name__ == "__main__":

    from rmhd.testbench import *
    #sr_shocktube()
    #compare_reconstruct()
    compare_quartic()
