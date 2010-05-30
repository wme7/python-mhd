#!/usr/bin/env python




class TestProblemDriver:

    def __init__(self, problem, Nx=128, CFL=0.8):

        from numpy import array

        to_array = lambda S: array([S['Rho'], S['Pre'],
                                    S['v'][0], S['v'][1], S['v'][2],
                                    S['B'][0], S['B'][1], S['B'][2]])

        self.Nx = Nx
        self.CFL = CFL
        self.L_state = to_array(problem.L_state)
        self.R_state = to_array(problem.R_state)
        self.problem = problem


    def run(self, lib, state, tfinal=0.2):

        from numpy import zeros

        self.verbose = False
        Nx = self.Nx
        P = zeros((Nx,8))
        U = zeros((Nx,8))
        L = zeros((Nx,8))

        for i in range(8):
            P[:Nx/2,i] = self.L_state[i]
            P[Nx/2:,i] = self.R_state[i]

        state.adiabatic_gamma = self.problem.adiabatic_gamma

        lib.set_state(state)
        lib.initialize(P,Nx)
        lib.prim_to_cons_array(P,U)

        dx  = 1.0 / Nx
        t   = 0.0
        dt  = self.CFL * dx

        while t < tfinal:
            e1 = lib.dUdt_1d(U,L)
            e2 = lib.dUdt_1d(U + 0.5*dt*L,L)

            if e1 or e2:
                print "Run crashed! Sorry..."
                break

            U += dt*L
            t += dt

            U[ 0:4,   :] = U[   4,:]
            U[Nx-4:Nx,:] = U[Nx-5,:]

            if self.verbose:
                print "t =", t

        lib.cons_to_prim_array(U,P)
        lib.finalize()
        return P




def sr_shocktube():

    from pylab import show
    import rmhd

    problem = SRShockTube2()
    driver  = TestProblemDriver(problem, Nx=100, CFL=0.8)
    state   = rmhd.LibraryState()

    P1 = driver.run(rmhd._lib1, state, tfinal=0.4)
    P2 = rmhd.riemann.exact_sr_vt(problem, tfinal=0.4)

    rmhd.visual.shocktube(P1, label="lib2", linestyle='--')
    rmhd.visual.shocktube(P2, label="exact", linestyle='-', marker='None', lw=2)
    show()


def compare_libs_1_and_2():

    from pylab import show
    import rmhd

    problem = RMHDShockTube1(L={'Pre':10.0})
    driver  = TestProblemDriver(problem, Nx=256, CFL=0.8)
    state1  = rmhd.LibraryState()
    state2  = rmhd.LibraryState(mode_reconstruct=1)

    P1 = driver.run(rmhd._lib1, state1)
    P2 = driver.run(rmhd._lib2, state2)

    rmhd.visual.shocktube(P1, label="run1", linestyle='--')
    rmhd.visual.shocktube(P2, label="run2", linestyle='-', mfc='None')
    show()


if __name__ == "__main__":

    from rmhd.testbench import *
    #sr_shocktube()
    compare_libs_1_and_2()

