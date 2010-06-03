#!/usr/bin/env python




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

    from rmhd import LibraryState, visual, _lib
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(512,), L=(1.0,))
    problem = RMHDContactWave()

    state0  = LibraryState(mode_riemann_solver=0, mode_reconstruct=2)
    state1  = LibraryState(mode_riemann_solver=1, mode_reconstruct=2)

    run_args = {'CFL':0.1, 'tfinal':0.2}

    P0 = driver.run(_lib, state0, problem, **run_args)
    P1 = driver.run(_lib, state1, problem, **run_args)

    visual.shocktube(P0, label="HLL", linestyle='--', mfc='None')
    visual.shocktube(P1, label="HLLC", linestyle='-', mfc='None')
    visual.show()



def riemann_wave_pattern():

    from rmhd import _lib, visual, LibraryState
    from numpy import zeros, array, linspace


    problem = RMHDRotationalWave()
    Pl, Pr = problem.get_states()

    Pl = array(Pl)
    Pr = array(Pr)

    Nx = 512
    x = linspace(-2.0,2.0,Nx)
    F, U, P, Ul, Ur = zeros(8),zeros(8),zeros(8),zeros(8),zeros(8)
    P_hllc, P_hll = zeros((Nx,8)),zeros((Nx,8))

    _lib.set_state(LibraryState(cons_to_prim_use_estimate=1))

    for i in range(Nx):

        U = zeros(8)

        _lib.hll_flux(Pl,Pr,U,F,x[i])
        if _lib.cons_to_prim_point(U,P_hll[i,:]):
            print "Warning! HLL  generated non-invertible intermediate cons state."

        _lib.hllc_flux(Pl,Pr,U,F,x[i])
        if _lib.cons_to_prim_point(U,P_hllc[i,:]):
            print "Warning! HLLC generated non-invertible intermediate cons state."

    state0 = LibraryState(mode_riemann_solver=0)
    state1 = LibraryState(mode_riemann_solver=1)

    run_args = {'Nx':Nx, 'CFL':0.2, 'tfinal':0.25}
    P_run  = run_1d_problem(_lib, state0, problem, **run_args)
    P_runc = run_1d_problem(_lib, state1, problem, **run_args)

    visual.shocktube(P_hll , x=(0,1), label="HLL" , linestyle='-.', marker='None', lw=6)
    visual.shocktube(P_hllc, x=(0,1), label="HLLC", linestyle='-' , marker='None')
    visual.shocktube(P_run , x=(0,1), label="run hll ", linestyle=':', marker='None', lw=2.5)
    visual.shocktube(P_runc, x=(0,1), label="run hllc", linestyle='-', marker='None', lw=0.5)
    visual.show()



def compare_quartic():

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
    rmhd.visual.show()



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


def divergence_2d_cross_stencil(f,g):

    from numpy import zeros

    Nx = f.shape[0]
    Ny = f.shape[1]

    dx = 1.0/Nx
    dy = 1.0/Ny
    div = zeros((Nx,Ny))

    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            div[i,j] = (f[i+1,j]+f[i+1,j+1]-f[i,j]-f[i,j+1])/(2*dx) + (g[i,j+1]+g[i+1,j+1]-g[i,j]-g[i+1,j])/(2*dy)

    return div



def test_2d():

    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    state = LibraryState(plm_theta=2.0, mode_reconstruct=2, mode_riemann_solver=0,
                         mode_quartic_solver=0)

    driver = ProblemDriver(N=(132,132), L=(2,2))

    problem = RMHDCylindricalA(pre=1.0)

    P = driver.run(_lib, state, problem, RK_order=3, CFL=0.6, tfinal=0.4)

    visual.four_pane_2dA(P, extent=[-1,1,-1,1])
    visual.show()


if __name__ == "__main__":

    from rmhd.testbench import *
    #sr_shocktube()
    #compare_riemann_solver()
    #compare_reconstruct()
    #compare_quartic()
    #riemann_wave_pattern()
    #library_dead_unit_test()
    test_2d()
