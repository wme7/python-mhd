#!/usr/bin/env python


def sr_shocktube():

    from rmhd import LibraryState, visual, _lib, riemann
    from rmhd.driver import ProblemDriver
    from pylab import figure

    driver = ProblemDriver(N=(256,), L=(1.0,))
    problem = SRShockTube3()
    state0 = LibraryState(mode_riemann_solver=1)
    state1 = LibraryState(mode_riemann_solver=1, mode_slope_limiter=0)

    PE = riemann.exact_sr_vt(     problem, tfinal=0.6)
    P0 = driver.run(_lib, state0, problem, tfinal=0.6, CFL=0.5, RK_order=3)
    P1 = driver.run(_lib, state1, problem, tfinal=0.6, CFL=0.5, RK_order=4)

    fig = figure()
    fig.text(.5, .95,
              r"Zhang & MacFadyen (2005) Section 6.1 on 256 zones"+
              r"   $t=0.6s$   CFL=$0.5$   $\theta_{PLM}=2.0$", 
              fontsize=14, horizontalalignment='center')

    visual.shocktube(PE, label="exact", linestyle='-', marker='None', lw=2)
    visual.shocktube(P0, label="HLL +RK3+PLM", linestyle='--')
    visual.shocktube(P1, label="HLLC+CTU+PLM", linestyle=':', lw=3)
    visual.show()



def compare_mlines_ctu():

    from pylab import show
    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(1024,), L=(1.0,))
    problem = RMHDShockTube4()
    state = LibraryState(plm_theta=2.0, mode_slope_limiter=0)

    run_args = {'CFL':0.8, 'tfinal':0.2}

    P0 = driver.run(_lib, state, problem, name='method of lines, RK3', RK_order=2, **run_args)
    P1 = driver.run(_lib, state, problem, name='CTU, first order    ', RK_order=4, **run_args)

    visual.shocktube(P0, label='method of lines, RK3', linestyle='--', mfc='None')
    visual.shocktube(P1, label='CTU, first order    ', linestyle='-.', mfc='None')
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

    run_args = {'CFL':0.5, 'tfinal':0.2, 'RK_order':3}

    P0 = driver.run(_lib, state0, problem, name='piecewise constant', **run_args)
    P1 = driver.run(_lib, state1, problem, name='PLM on 3-velocity ', **run_args)
    P2 = driver.run(_lib, state2, problem, name='PLM on 4-velocity ', **run_args)

    visual.shocktube(P0, label="Piecewise Constant", linestyle='--', mfc='None')
    visual.shocktube(P1, label="PLM 3-velocity", linestyle='-.', mfc='None')
    visual.shocktube(P2, label="PLM 4-velocity", linestyle='-', marker='None')
    visual.show()



def compare_limiter():

    from pylab import show
    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    driver = ProblemDriver(N=(512,), L=(1.0,))
    problem = RMHDShockTube2()

    stateA  = LibraryState(mode_reconstruct=0)
    state0  = LibraryState(mode_slope_limiter=0)
    state1  = LibraryState(mode_slope_limiter=1)
    state2  = LibraryState(mode_slope_limiter=2)

    run_args = {'CFL':0.5, 'tfinal':0.2, 'RK_order':4}

    PA = driver.run(_lib, stateA, problem, **run_args)
    P0 = driver.run(_lib, state0, problem, **run_args)
    P1 = driver.run(_lib, state1, problem, **run_args)
    P2 = driver.run(_lib, state2, problem, **run_args)

    visual.shocktube(PA, label="no reconstruct", linestyle='-', marker='None')
    visual.shocktube(P0, label="minmod", linestyle='--', mfc='None')
    visual.shocktube(P1, label="MC", linestyle='-.', mfc='None')
    visual.shocktube(P2, label="harmonic mean", linestyle=':', lw=3, marker='None')
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

    B = {'B':[0,0,0]}

    driver = ProblemDriver(N=(132,132), L=(2,2))
    problem = RMHDCylindricalA(pre=100.0, I=B, O=B)
    state = LibraryState(plm_theta=2.0, mode_reconstruct=2, mode_riemann_solver=1)

    run_args = {'name': "cylindrical blast wave",
                'RK_order': 3, 'CFL': 0.6, 'tfinal': 0.2}

    P = driver.run(_lib, state, problem, **run_args)

    visual.four_pane_2d(P, extent=[-1,1,-1,1])
    visual.show()



def symmetry_test():

    """
    This test demonstrates an issue in the library where reconstruction on
    4-velocities results in some unexpected asymmetry, in a model with
    initial conditions having mirror symmetry across the diagonal. Note that
    the effect is either fixed or 'masked' when sqrtf is used in place of
    sqrt in 4-velocity reconstruction function.
    """

    from rmhd import _lib, LibraryState
    from rmhd.driver import ProblemDriver
    from pylab import imshow, colorbar, show, subplot, title
    from numpy import flipud

    driver = ProblemDriver(N=(32,32), L=(2,2))
    problem = SRQuadrantA()

    quartic     = [0,3,0,3]
    reconstruct = [1,1,2,2]
    names       = ['3vel, quartic exact', '3vel, quartic none',
                   '4vel, quartic exact', '4vel, quartic none']
    runnum      = range(1,5)
    run_args    = {'RK_order': 3, 'CFL': 0.4, 'tfinal': 0.2}

    def do_sym(P, name):

        tr = lambda x: flipud(x.T)

        Nx, Ny, Nq = P.shape
        q = P[2:Nx-2,2:Ny-2,0]
        imshow(tr((q-q.T)/(q+q.T)))
        colorbar()
        title(name)


    for num,name,q,r in zip(runnum, names, quartic, reconstruct):

        state = LibraryState(mode_quartic_solver=q, mode_reconstruct=r)
        P = driver.run(_lib, state, problem, name=name, **run_args)

        subplot(2,2,num)
        do_sym(P, name)

    show()


def quadrant_problem():

    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver
    from pylab import figure, savefig

    probs = [SRQuadrantA(), SRQuadrantB()]
    names = ['A', 'B']

    driver = ProblemDriver(N=(32,32), L=(2,2))
    state = LibraryState()

    run_args = {'RK_order': 2, 'CFL': 0.6, 'tfinal': 0.8}
    caption = r"Quadrant Problem %s after $t=0.8s$  $N=256^2$  PLM=$2.0$  CFL=$0.6$  RK=$2$  HLL"

    for prob, name in zip(probs, names):

        P = driver.run(_lib, state, prob, **run_args)
        fig = figure(figsize=(10,12))
        visual.contour_2d(P[:,:,0], fig, extent=[-1,1,-1,1], caption=caption % name)
        savefig('quadrant%s.png' % name)



def cross_stencil_div_3d(fx,fy,fz,dx=1,dy=1,dz=1):

    from numpy import zeros
    assert fx.shape == fy.shape and fy.shape == fz.shape
    Nx, Ny, Nz = fx.shape

    div = zeros((Nx,Ny,Nz))
    for i in range(1,Nx-1):
        for j in range(1,Ny-1):
            for k in range(1,Nz-1):

                div[i,j,k] = \
                ((fx[i+1,j,k] + fx[i+1,j+1,k] + fx[i+1,j,k+1] + fx[i+1,j+1,k+1]) - (fx[i,j,k] + fx[i,j+1,k] + fx[i,j,k+1] + fx[i,j+1,k+1])) / (4*dx) + \
                ((fy[i,j+1,k] + fy[i,j+1,k+1] + fy[i+1,j+1,k] + fy[i+1,j+1,k+1]) - (fy[i,j,k] + fy[i,j,k+1] + fy[i+1,j,k] + fy[i+1,j,k+1])) / (4*dy) + \
                ((fz[i,j,k+1] + fz[i+1,j,k+1] + fz[i,j+1,k+1] + fz[i+1,j+1,k+1]) - (fz[i,j,k] + fz[i+1,j,k] + fz[i,j+1,k] + fz[i+1,j+1,k])) / (4*dz)

    return div


def spherical_blast_3d():

    from rmhd import _lib, LibraryState, visual
    from rmhd.driver import ProblemDriver

    B = {'B':[0.8,0,0]}

    problem = RMHDCylindricalA( I=B, O=B )
    driver = ProblemDriver(N=(32,32,32), L=(2,2,2))
    state = LibraryState(mode_riemann_solver=1, mode_reconstruct=1)

    run_args = {'RK_order': 3, 'CFL': 0.05, 'tfinal': 0.3}
    P = driver.run(_lib, state, problem, **run_args)
    visual.four_pane_2d(P[:,:,16], extent=[-1,1,-1,1])
    visual.show()



def analyze_failed_state(pickle_name):

    from pickle import load
    from numpy import zeros_like
    from rmhd import _lib, LibraryState


    e = load(open(pickle_name))

    P = e.SurroundingBlock
    F = zeros_like(P)

    state = LibraryState(mode_riemann_solver=1, mode_reconstruct=1)

    _lib.set_state(state)
    _lib.initialize(P, 5,5,5, 1.0,1.0,1.0, 0)

    _lib.set_dimension(1)
    _lib.Fiph(P,F)
    print "Fx:", F[1,1,1], F[2,2,2]

    _lib.set_dimension(2)
    _lib.Fiph(P,F)
    print "Fy:", F[1,1,1], F[2,2,2]

    _lib.set_dimension(3)
    _lib.Fiph(P,F)
    print "Fz:", F[1,1,1], F[2,2,2]

    _lib.finalize()





if __name__ == "__main__":

    from rmhd.testbench import *
    from pylab import figure, zeros, imshow, show

    from optparse import OptionParser

    parser = OptionParser()
    opt, args = parser.parse_args()

    if len(args) and args[0].endswith('.fail'):
        analyze_failed_state(args[0])

    else:
        #sr_shocktube()
        #riemann_wave_pattern()
        cylindrical_blast()
        #quadrant_problem()
        #spherical_blast_3d()
        #symmetry_test()

        #compare_riemann_solver()
        #compare_reconstruct()
        #compare_limiter()
        #compare_quartic()
        #compare_mlines_ctu()
