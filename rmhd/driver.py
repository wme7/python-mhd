

class ProblemDriver:


    def __init__(self, N=None, L=None):

        self.N = N
        self.L = L
        self.Ng = 2
        self.failures = 0


    def main_loop(self, U):

        dUdt = { 1: self.dUdt_1d, 2: self.dUdt_2d, 3: self.dUdt_3d }[len(self.N)]
        t   = 0.0
        dt  = 1e-9
        n_cycle = 0

        from time import time
        while t < self.tfinal:

            n_cycle += 1
            start = time()

            if self.RK_order is 1:
                U += dt*dUdt(U)

            if self.RK_order is 2:
                L1 =    dUdt(U);
                U += dt*dUdt(U + 0.5*dt*L1);

            elif self.RK_order is 3:

                U1 =      U +                  dt * dUdt(U )
                U1 = 3./4*U + 1./4*U1 + 1./4 * dt * dUdt(U1)
                U  = 1./3*U + 2./3*U1 + 2./3 * dt * dUdt(U1)

            t += dt
            step_time = time()-start

            msg_data = (n_cycle, t, dt, 1e6*step_time/U.size, self.failures)
            msg_text = "N: %05d t: %6.4f dt: %6.4e us/zone: %5.4f failures: %d"
            if not self.quiet: print msg_text % msg_data

            dt = self.CFL * self.min_dx
            self.failures = 0
        return U


    def run(self, lib, state, problem, name=None, **kwargs):

        from time import time

        if name is None: name = problem.__class__

        run_func = { 1: self.run_1d, 2: self.run_2d, 3: self.run_3d }
        start_time = time()
        P = run_func[len(self.N)](lib, state, problem, **kwargs)
        print "Python driver finished '%s'... total time: %f" % (name, time() - start_time)

        return P

    def run_1d(self, lib, state, problem, RK_order=2, CFL=0.5, tfinal=0.2):

        from numpy import zeros

        self.lib = lib
        Nx, = self.N
        Lx, = self.L

        P = zeros((Nx,8))
        U = zeros((Nx,8))

        state.adiabatic_gamma = problem.adiabatic_gamma
        problem.initial_model(P)

        lib.set_state(state)
        lib.initialize(P, Nx, 1, 1, Lx, 0.0, 0.0, 1)
        lib.prim_to_cons_array(P,U,U.size/8)

        dx = 1.0 * Lx / (Nx-2*self.Ng);

        self.CFL = CFL
        self.tfinal = tfinal
        self.RK_order = RK_order
        self.min_dx = dx
        self.quiet = True
        U = self.main_loop(U.copy())

        lib.cons_to_prim_array(U,P,U.size/8)
        lib.finalize()
        return P


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
        lib.initialize(P, Nx, Ny, 1, Lx, Ly, 0.0, 0)
        lib.prim_to_cons_array(P,U,U.size/8)

        dx = 1.0 * Lx / (Nx-2*self.Ng);
        dy = 1.0 * Ly / (Ny-2*self.Ng);

        self.CFL = CFL
        self.tfinal = tfinal
        self.RK_order = RK_order
        self.min_dx = min([dx,dy])
        self.quiet = False
        U = self.main_loop(U)

        lib.cons_to_prim_array(U,P,U.size/8)
        lib.finalize()
        return P


    def run_3d(self, lib, state, problem, RK_order=2, CFL=0.5, tfinal=0.2):

        from numpy import zeros

        self.lib = lib
        Nx,Ny,Nz = self.N
        Lx,Ly,Lz = self.L

        P = zeros((Nx,Ny,Nz,8))
        U = zeros((Nx,Ny,Nz,8))

        state.adiabatic_gamma = problem.adiabatic_gamma
        problem.initial_model(P)

        lib.set_state(state)
        lib.initialize(P, Nx, Ny, Nz, Lx, Ly, Lz, 0)
        lib.prim_to_cons_array(P,U,U.size/8)

        dx = 1.0 * Lx / (Nx-2*self.Ng)
        dy = 1.0 * Ly / (Ny-2*self.Ng)
        dz = 1.0 * Lz / (Nz-2*self.Ng)

        self.CFL = CFL
        self.tfinal = tfinal
        self.RK_order = RK_order
        self.min_dx = min([dx,dy,dz])
        self.quiet = False
        U = self.main_loop(U)

        lib.cons_to_prim_array(U,P,U.size/8)
        lib.finalize()
        return P



    def dUdt_1d(self, U):

        from numpy import zeros_like
        L = zeros_like(U)

        Ng = self.Ng
        Nx, = self.N

        U[ 0:Ng,   ] = U[   Ng  ] # Boundary conditions
        U[Nx-Ng:Nx,] = U[Nx-Ng-1]

        self.failures += self.lib.dUdt_1d(U,L)
        return L


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


    def dUdt_3d(self, U):

        from numpy import zeros_like
        L = zeros_like(U)

        Ng = self.Ng
        Nx,Ny,Nz = self.N

        for i in range(Ng): # Boundary conditions
            U[   i  ,:,:,:] = U[   Ng  ,:,:,:]
            U[Nx-i-1,:,:,:] = U[Nx-Ng-1,:,:,:]

            U[:,   i  ,:,:] = U[:,   Ng  ,:,:]
            U[:,Ny-i-1,:,:] = U[:,Ny-Ng-1,:,:]

            U[:,:,   i  ,:] = U[:,:,   Ng  ,:]
            U[:,:,Nz-i-1,:] = U[:,:,Nz-Ng-1,:]

        self.failures += self.lib.dUdt_3d(U,L)
        return L

