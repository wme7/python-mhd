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

    while t < tfinal:
        e1 = lib.dUdt_1d(U,L)
        e2 = lib.dUdt_1d(U + 0.5*dt*L,L)

        if e1 or e2:
            print "Run crashed! Sorry...", e1, e2
            break

        U += dt*L
        t += dt

        U[ 0:4,   :] = U[   4,:]
        U[Nx-4:Nx,:] = U[Nx-5,:]

        if verbose:
            print "t =", t

    lib.cons_to_prim_array(U,P)
    lib.finalize()
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




def compare_libs_1_and_2():

    from pylab import show
    import rmhd

    problem = RMHDShockTube1(L={'Pre':1.0})

    state1  = rmhd.LibraryState()
    state2  = rmhd.LibraryState(mode_reconstruct=2)

    P1 = run_1d_problem(rmhd._lib1, state1, problem, CFL=0.5, tfinal=0.4)
    P2 = run_1d_problem(rmhd._lib2, state2, problem, CFL=0.5, tfinal=0.4)

    rmhd.visual.shocktube(P1, label="run1", linestyle='--')
    rmhd.visual.shocktube(P2, label="run2", linestyle='-', mfc='None')
    show()


if __name__ == "__main__":

    from rmhd.testbench import *
    #sr_shocktube()
    compare_libs_1_and_2()

