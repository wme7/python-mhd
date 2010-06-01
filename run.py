#!/usr/bin/env python


def riemann_wave_patter():

    from rmhd import _lib, visual, LibraryState
    from numpy import zeros, array, linspace
    from pylab import show

    problem = RMHDShockTube2(L={'B':[1e-4,2.0,2.0]}, R={'B':[1e-4,-2.0,2.0]})
    Pl, Pr = problem.get_states()

    Pl = array(Pl)
    Pr = array(Pr)

    Nx = 1024
    x = linspace(-1.5,1.5,Nx)
    F, U, P, Ul, Ur = zeros(8),zeros(8),zeros(8),zeros(8),zeros(8)
    P_hllc, P_hll = zeros((Nx,8)),zeros((Nx,8))

    _lib.set_state(LibraryState(cons_to_prim_use_estimate=1))

    for i in range(Nx):

        U = zeros(8)

        _lib.hll_flux(Pl,Pr,U,F,x[i])
        if _lib.cons_to_prim_point(U,P_hll[i,:]):
            print "Warning! HLL generated non-invertible intermediate cons state."

        _lib.hllc_flux(Pl,Pr,U,F,x[i])
        if _lib.cons_to_prim_point(U,P_hllc[i,:]):
            print "Warning! HLLC generated non-invertible intermediate cons state."


    visual.shocktube(P_hll , x=(-1.5,1.5), label="HLL" , linestyle='-.', marker='None', lw=6)
    visual.shocktube(P_hllc, x=(-1.5,1.5), label="HLLC", linestyle='-' , marker='None')
    show()



def run_1d_problem(lib, state, problem, Nx=128, CFL=0.5, tfinal=0.2, verbose=False):

    from numpy import zeros

    P = zeros((Nx,8))
    U = zeros((Nx,8))
    L = zeros((Nx,8))

    state.adiabatic_gamma = problem.adiabatic_gamma
    problem.initial_model(P)
    lib.set_state(state)
    lib.initialize(P,Nx,1,1)
    lib.prim_to_cons_array(P,U,Nx)

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

    lib.cons_to_prim_array(U,P,Nx)
    lib.finalize()
    print "Solver averaged %f us/zone" % (ttltime / (Nx*8*n_cycle)*1e6)
    return P



def run_2d_problem(lib, state, problem, Nx=128, Ny=128, CFL=0.5, tfinal=0.2, verbose=True):

    from numpy import zeros

    P = zeros((Nx,Ny,8))
    U = zeros((Nx,Ny,8))
    L = zeros((Nx,Ny,8))

    state.adiabatic_gamma = 1.4
    problem.initial_model(P)

    lib.set_state(state)
    lib.initialize(P,Nx,Ny,1)
    lib.prim_to_cons_array(P,U,Nx*Ny)

    dx  = 1.0 / Nx
    dy  = 1.0 / Ny
    t   = 0.0
    dt  = CFL * min([dx,dy])

    from time import time

    ttltime = 0.0
    n_cycle = 0
    while t < tfinal:

        start = time()

        e1 = lib.dUdt_2d(U,L)
        e2 = lib.dUdt_2d(U + 0.5*dt*L,L)

        if e1 or e2:
            print "Run crashed! Sorry...", e1, e2
            break
        U += dt*L
        t += dt
        n_cycle += 1

        for i in range(4): # Boundary conditions
            U[   i  ,:] = U[   4,:]
            U[Nx-i-1,:] = U[Nx-5,:]

            U[:,   i  ] = U[:,   4]
            U[:,Ny-i-1] = U[:,Ny-5]

        ttltime += time()-start

        if verbose:
            print "t =", t

    lib.cons_to_prim_array(U,P,Nx*Ny)
    lib.finalize()
    print "Solver averaged %f us/zone" % (ttltime / (Nx*Ny*8*n_cycle)*1e6)
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

    problem = RMHDShockTube2()

    state0  = rmhd.LibraryState(mode_reconstruct=0)
    state1  = rmhd.LibraryState(mode_reconstruct=1)
    state2  = rmhd.LibraryState(mode_reconstruct=2)

    run_args = {'Nx':256, 'CFL':0.5, 'tfinal':0.2}

    P0 = run_1d_problem(rmhd._lib, state0, problem, **run_args)
    P1 = run_1d_problem(rmhd._lib, state1, problem, **run_args)
    P2 = run_1d_problem(rmhd._lib, state2, problem, **run_args)

    rmhd.visual.shocktube(P0, label="Piecewise Constant", linestyle='--', mfc='None')
    rmhd.visual.shocktube(P1, label="PLM 3-velocity", linestyle='-.', mfc='None')
    rmhd.visual.shocktube(P2, label="PLM 4-velocity", linestyle='-', marker='None')
    show()


def compare_riemann_solver():

    from pylab import show
    import rmhd

    problem = RMHDContactWave()

    state0  = rmhd.LibraryState(mode_riemann_solver=0, mode_reconstruct=2)
    state1  = rmhd.LibraryState(mode_riemann_solver=1, mode_reconstruct=2)

    run_args = {'Nx':512, 'CFL':0.1, 'tfinal':0.2}

    P0 = run_1d_problem(rmhd._lib, state0, problem, **run_args)
    P1 = run_1d_problem(rmhd._lib, state1, problem, **run_args)

    rmhd.visual.shocktube(P0, label="HLL", linestyle='--', mfc='None')
    rmhd.visual.shocktube(P1, label="HLLC", linestyle='-', mfc='None')
    show()



def compare_quartic():

    from pylab import show
    import rmhd

    problem = RMHDShockTube4()

    state0  = rmhd.LibraryState(mode_quartic_solver=0)
    state1  = rmhd.LibraryState(mode_quartic_solver=1)
    state2  = rmhd.LibraryState(mode_quartic_solver=2)
    state3  = rmhd.LibraryState(mode_quartic_solver=3)

    run_args = {'Nx':512, 'CFL':0.5, 'tfinal':0.2}

    P0 = run_1d_problem(rmhd._lib, state0, problem, **run_args)
    P1 = run_1d_problem(rmhd._lib, state1, problem, **run_args)
    P2 = run_1d_problem(rmhd._lib, state2, problem, **run_args)
    P3 = run_1d_problem(rmhd._lib, state3, problem, **run_args)

    rmhd.visual.shocktube(P0, label="Exact", linestyle='--', mfc='None')
    rmhd.visual.shocktube(P1, label="Approx1", linestyle='-.', mfc='None')
    rmhd.visual.shocktube(P2, label="Approx2", linestyle=':', marker='None')
    rmhd.visual.shocktube(P3, label="None", linestyle='-', marker='None')
    show()



def library_dead_unit_test():

    import rmhd
    from numpy import array, zeros

    problem = RMHDShockTube1()
    P = array(problem.get_states()[0])
    U = zeros(8)

    passfail = lambda p: 'Fail' if p else 'Pass'

    print "Testing prim_to_cons_point"
    print "\tOnly option:"       , passfail(rmhd._lib.prim_to_cons_point(P,U))
    print "Testing cons_to_prim_point..."
    rmhd._lib.set_state(rmhd.LibraryState(cons_to_prim_use_estimate=1))
    print "\tWith estimate:"     , passfail(rmhd._lib.cons_to_prim_point(U,P))
    rmhd._lib.set_state(rmhd.LibraryState(cons_to_prim_use_estimate=0))
    print "\tWith good guess:"   , passfail(rmhd._lib.cons_to_prim_point(U,P))
    print "\tWith bad guess:"    , passfail(rmhd._lib.cons_to_prim_point(U,P*0.1))
    print "\tWith aweful guess:" , passfail(rmhd._lib.cons_to_prim_point(U,P*0.0))

    print "Testing prim_to_cons_array..."
    Nx = 100

    U_all = zeros((Nx,8))
    P_all = zeros((Nx,8))
    for i in range(Nx):
        P_all[i,:] = P

    print "Testing prim_to_cons_array"
    print "\tOnly option:"       , passfail(rmhd._lib.prim_to_cons_array(P_all,U_all,Nx))
    print "Testing cons_to_prim_array..."
    rmhd._lib.set_state(rmhd.LibraryState(cons_to_prim_use_estimate=1))
    print "\tWith estimate:"     , passfail(rmhd._lib.cons_to_prim_array(U_all,P_all,Nx))
    rmhd._lib.set_state(rmhd.LibraryState(cons_to_prim_use_estimate=0))
    print "\tWith good guess:"   , passfail(rmhd._lib.cons_to_prim_array(U_all,P_all,Nx))
    print "\tWith bad guess:"    , passfail(rmhd._lib.cons_to_prim_array(U_all,P_all*0.1,Nx))
    print "\tWith aweful guess:" , passfail(rmhd._lib.cons_to_prim_array(U_all,P_all*0.0,Nx))


def test_2d():

    from rmhd import _lib, LibraryState, visual
    from numpy import array, zeros
    from pylab import show, imshow, colorbar

    state = LibraryState(plm_theta=2.0, mode_reconstruct=2, mode_riemann_solver=1)
    problem = SRShockTube1()#RMHDCylindricalA(pre=0.01)
    problem.orientation = 'y'
    P = run_2d_problem(_lib, state, problem, Nx=128, Ny=128, CFL=0.5, tfinal=0.2, verbose=True)
    visual.four_pane_2d(P, x=(-1,1))
    show()


if __name__ == "__main__":

    from rmhd.testbench import *
    #sr_shocktube()
    #compare_riemann_solver()
    #compare_reconstruct()
    #compare_quartic()
    #riemann_wave_patter()
    #library_dead_unit_test()
    test_2d()
