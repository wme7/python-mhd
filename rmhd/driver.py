

class ProblemDriver:

    def __init__(self, N=None, L=None):

        self.N = N
        self.L = L
        self.Ng = 2
        self.failures = 0


    def dUdt_2d(self, U):

        from numpy import zeros_like
        L = zeros_like(U)

        Ng = self.Ng
        Nx,Ny = self.N

        for i in range(Ng): # Boundary conditions
            U[   i  ,:] = U[   Ng  ,:]
            U[Nx-i-1,:] = U[Nx-Ng-1,:]

            U[:,   i  ] = U[:,   Ng  ]
            U[:,Ny-i-1] = U[:,Ny-Ng-1]

        self.failures += self.lib.dUdt_2d(U,L)
        return L


    def run(self, lib, state, problem, **kwargs):

        run_func = { 1: self.run_1d, 2: self.run_2d, 3: self.run_3d }
        return run_func[len(self.N)](lib, state, problem, **kwargs)


    def run_1d(self, lib, state, problem, CFL=0.5, tfinal=0.2):
        pass

    def run_3d(self, lib, state, problem, CFL=0.5, tfinal=0.2):
        pass

    def run_2d(self, lib, state, problem, RK_order=2, CFL=0.5, tfinal=0.2):

        from numpy import zeros

        self.lib = lib
        Nx,Ny = self.N
        Lx,Ly = self.L

        P = zeros((Nx,Ny,8))
        U = zeros((Nx,Ny,8))

        state.adiabatic_gamma = problem.adiabatic_gamma
        problem.initial_model(P)

        lib.set_state(state)
        lib.initialize(P,Nx,Ny,1,Lx,Ly,0.0)
        lib.prim_to_cons_array(P,U,Nx*Ny)

        dx = 1.0 * Lx / (Nx-4);
        dy = 1.0 * Ly / (Ny-4);

        t   = 0.0
        dt  = 1e-9

        run_time = 0.0
        n_cycle = 0

        dUdt = self.dUdt_2d

        from time import time
        while t < tfinal:

            n_cycle += 1
            start = time()

            if RK_order is 1:
                U += dt*dUdt(U)

            if RK_order is 2:
                L1 =    dUdt(U);
                U += dt*dUdt(U + 0.5*dt*L1);

            elif RK_order is 3:

                U1 =      U +                  dt * dUdt(U )
                U1 = 3./4*U + 1./4*U1 + 1./4 * dt * dUdt(U1)
                U  = 1./3*U + 2./3*U1 + 2./3 * dt * dUdt(U1)

            t += dt
            step_time = time()-start
            run_time += step_time

            msg_data = (n_cycle, t, dt, 1e6*step_time/U.size, self.failures)
            msg_text = "N: %05d t: %6.4f dt: %6.4e us/zone: %5.4f failures: %d"
            print msg_text % msg_data

            dt = CFL * min([dx,dy])
            self.failures = 0


        lib.cons_to_prim_array(U,P,U.size/8)
        lib.finalize()
        print "Finished! total time = %3.2fs" % run_time
        return P

