#!/usr/bin/env python



def sr_shocktube():

    from rmhd import LibraryState, visual, _lib, riemann
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(1024,), L=(1.0,))
    problem = SRShockTube2()
    state = LibraryState()

    P0 = riemann.exact_sr_vt(    problem, tfinal=0.2)
    P1 = driver.run(_lib, state, problem, tfinal=0.2, CFL=0.5)

    visual.shocktube(P0, label="exact", linestyle='-', marker='None', lw=2)
    visual.shocktube(P1, label="HLL", linestyle='--')
    visual.show()



def compare_reconstruct():

    from pylab import show
    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(1024,), L=(1.0,))
    problem = RMHDShockTube2()

    state0  = LibraryState(mode_reconstruct=0)
    state1  = LibraryState(mode_reconstruct=1)
    state2  = LibraryState(mode_reconstruct=2)

    run_args = {'CFL':0.5, 'tfinal':0.2}

    P0 = driver.run(_lib, state0, problem, **run_args)
    P1 = driver.run(_lib, state1, problem, **run_args)
    P2 = driver.run(_lib, state2, problem, **run_args)

    visual.shocktube(P0, label="Piecewise Constant", linestyle='--', mfc='None')
    visual.shocktube(P1, label="PLM 3-velocity", linestyle='-.', mfc='None')
    visual.shocktube(P2, label="PLM 4-velocity", linestyle='-', marker='None')
    visual.show()



def compare_riemann_solver():

    from rmhd import LibraryState, visual, _lib
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(1024,), L=(1.0,))
    problem = RMHDShockTube1()
    state0  = LibraryState(mode_riemann_solver=0, mode_reconstruct=2)
    state1  = LibraryState(mode_riemann_solver=1, mode_reconstruct=2)

    run_args = {'CFL':0.8, 'tfinal':0.2}

    P0 = driver.run(_lib, state0, problem, **run_args)
    P1 = driver.run(_lib, state1, problem, **run_args)

    visual.shocktube(P0, label="HLL", linestyle='--', mfc='None')
    visual.shocktube(P1, label="HLLC", linestyle='-', mfc='None')
    visual.show()



def riemann_wave_pattern():

    from numpy import zeros, array, linspace
    from rmhd import _lib, visual, LibraryState
    from rmhd.driver import ProblemDriver


    Nx = 512
    x = linspace(-2.0,2.0,Nx)
    F, U, P, Ul, Ur = zeros(8),zeros(8),zeros(8),zeros(8),zeros(8)
    P_hllc, P_hll = zeros((Nx,8)),zeros((Nx,8))

    driver = ProblemDriver(N=(512,), L=(1.0,))
    problem = RMHDRotationalWave()
    Pl, Pr = problem.get_states()

    Pl = array(Pl)
    Pr = array(Pr)

    state = LibraryState(cons_to_prim_use_estimate=1)
    _lib.set_state(state)

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

    run_args = {'CFL':0.4, 'tfinal':0.25}
    P_run  = driver.run(_lib, state0, problem, name="HLL ", **run_args)
    P_runc = driver.run(_lib, state1, problem, name="HLLC", **run_args)

    visual.shocktube(P_hll , x=(0,1), label="HLL" , linestyle='-.', marker='None', lw=6)
    visual.shocktube(P_hllc, x=(0,1), label="HLLC", linestyle='-' , marker='None')
    visual.shocktube(P_run , x=(0,1), label="run hll ", linestyle=':', marker='None', lw=2.5)
    visual.shocktube(P_runc, x=(0,1), label="run hllc", linestyle='-', marker='None', lw=0.5)
    visual.show()



def compare_quartic():

    from rmhd import _lib, visual, LibraryState
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(512,), L=(1.0,))
    problem = RMHDShockTube4()

    state0  = LibraryState(mode_quartic_solver=0)
    state1  = LibraryState(mode_quartic_solver=1)
    state2  = LibraryState(mode_quartic_solver=2)
    state3  = LibraryState(mode_quartic_solver=3)

    run_args = {'CFL':0.5, 'tfinal':0.2}

    print "Testing time of different quartic solver modes..."

    P0 = driver.run(_lib, state0, problem, name="exact   ", **run_args)
    P1 = driver.run(_lib, state1, problem, name="newton-1", **run_args)
    P2 = driver.run(_lib, state2, problem, name="newton-2", **run_args)
    P3 = driver.run(_lib, state3, problem, name="no solve", **run_args)

    visual.shocktube(P0, label="Exact", linestyle='--', mfc='None')
    visual.shocktube(P1, label="Approx1", linestyle='-.', mfc='None')
    visual.shocktube(P2, label="Approx2", linestyle=':', marker='None')
    visual.shocktube(P3, label="None", linestyle='-', marker='None')
    visual.show()



def cylindrical_blast():

    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(132,132), L=(2,2))
    problem = RMHDCylindricalA(pre=100.0)
    state = LibraryState(plm_theta=2.0, mode_reconstruct=2, mode_riemann_solver=1)

    run_args = {'name': "cylindrical blast wave",
                'RK_order': 3, 'CFL': 0.6, 'tfinal': 0.2}

    P = driver.run(_lib, state, problem, **run_args)

    visual.four_pane_2d(P, extent=[-1,1,-1,1])
    visual.show()



def quadrant_test():

    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(64,64), L=(2,2))
    problem = SRQuadrantA()

    states = [LibraryState(mode_reconstruct    = 0,
                           mode_quartic_solver = 0),

              LibraryState(mode_reconstruct    = 0,
                           mode_quartic_solver = 3),

              LibraryState(mode_reconstruct    = 2,
                           mode_quartic_solver = 0),

              LibraryState(mode_reconstruct    = 2,
                           mode_quartic_solver = 3)]

    run_args = {'name': "quadrant test",
                'RK_order': 3, 'CFL': 0.2, 'tfinal': 0.2}

    for n,state in enumerate(states):

        state.plm_theta = 0.0
        P = driver.run(_lib, state, problem, **run_args)
        P.dump('quadrant%d.np' % (n+1))


    visual.four_pane_2d(P, extent=[-1,1,-1,1], do_quiver=False)
    visual.show()




if __name__ == "__main__":

    from rmhd.testbench import *
    from pylab import figure

    #sr_shocktube()
    #riemann_wave_pattern()
    #cylindrical_blast()
    quadrant_test()

    #compare_riemann_solver()
    #compare_reconstruct()
    #compare_quartic()
