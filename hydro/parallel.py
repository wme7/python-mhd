

class DecomposedDomain():

    def __init__(self, N=(256,), x0=(0.0,), x1=(1.0,)):
        assert len(N) is len(x0) is len(x1)
        from mpi4py.MPI import COMM_WORLD, Compute_dims

        periods = [True for n in N]
        mpi_sizes = Compute_dims(COMM_WORLD.size, len(N))
        cart = COMM_WORLD.Create_cart(mpi_sizes, periods=periods)
        mpi_coord = cart.Get_coords(COMM_WORLD.rank)

        global_shape = [n for n in N]
        global_start = [0 for n in N]

        local_shape = [ ]
        X0, X1 = [ ], [ ]
        dx = [float(l1-l0)/n for n,l0,l1 in zip(N,x0,x1)]

        for i in range(len(N)):

            R = N[i] % mpi_sizes[i]
            normal_size = N[i] / mpi_sizes[i]
            augmnt_size = normal_size + 1
            thisdm_size = augmnt_size if mpi_coord[i]<R else normal_size

            for j in range(mpi_coord[i]):
                global_start[i] += augmnt_size if j<R else normal_size

            local_shape.append(thisdm_size)

            X0.append(x0[i] + dx[i] *  global_start[i])
            X1.append(x0[i] + dx[i] * (global_start[i] + thisdm_size))

        self.N = local_shape
        self.dx = dx
        self.x0 = X0
        self.x1 = X1
        self.cart = cart
        self.rank = COMM_WORLD.rank
        self.mpi_coord = mpi_coord
        self.mpi_sizes = mpi_sizes


    def synchronize(self, A, Ng):
        from mpi4py.MPI import COMM_WORLD, Compute_dims

        if len(A.shape) == 2:
            L,R = self.cart.Shift(0,1)
            A[:+Ng] = self.cart.Sendrecv(A[-2*Ng:-Ng], dest=R, source=L)
            A[-Ng:] = self.cart.Sendrecv(A[+Ng:+2*Ng], dest=L, source=R)

        if len(A.shape) == 3:
            L,R = self.cart.Shift(1,1)
            A[:,:+Ng] = self.cart.Sendrecv(A[:,-2*Ng:-Ng], dest=R, source=L)
            A[:,-Ng:] = self.cart.Sendrecv(A[:,+Ng:+2*Ng], dest=L, source=R)

            L,R = self.cart.Shift(0,1)
            A[:+Ng,:] = self.cart.Sendrecv(A[-2*Ng:-Ng,:], dest=R, source=L)
            A[-Ng:,:] = self.cart.Sendrecv(A[+Ng:+2*Ng,:], dest=L, source=R)

        if len(A.shape) == 4:
            L,R = self.cart.Shift(0,1)
            A[:+Ng,:,:] = self.cart.Sendrecv(A[-2*Ng:-Ng,:,:], dest=R, source=L)
            A[-Ng:,:,:] = self.cart.Sendrecv(A[+Ng:+2*Ng,:,:], dest=L, source=R)

            L,R = self.cart.Shift(1,1)
            A[:,:+Ng,:] = self.cart.Sendrecv(A[:,-2*Ng:-Ng,:], dest=R, source=L)
            A[:,-Ng:,:] = self.cart.Sendrecv(A[:,+Ng:+2*Ng,:], dest=L, source=R)

            L,R = self.cart.Shift(2,1)
            A[:,:,:+Ng] = self.cart.Sendrecv(A[:,:,-2*Ng:-Ng], dest=R, source=L)
            A[:,:,-Ng:] = self.cart.Sendrecv(A[:,:,+Ng:+2*Ng], dest=L, source=R)


    def set_BC(self, A, Ng, BC=None):
        self.synchronize(A, Ng)
        if BC is None: return

        L_BCs = BC.L_wall[len(A.shape)-2]
        R_BCs = BC.R_wall[len(A.shape)-2]

        for i,BC in enumerate(L_BCs):
            if self.mpi_coord[i] == 0:
                BC(A, Ng)

        for i,BC in enumerate(R_BCs):
            if self.mpi_coord[i] == self.mpi_sizes[i]-1:
                BC(A, Ng)


    def dump(self, P, base='dump'):
        fmt = '-'.join(['%03d' for m in self.mpi_coord]) % self.mpi_coord
        P.dump(base+'_'+fmt+'.pk')
