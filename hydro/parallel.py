

class DecomposedDomain():

    def __init__(self, N, x0, x1, Ng):

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
        dx = [1.0*(l1-l0)/(n-2*Ng) for n,l0,l1 in zip(N,x0,x1)]

        for i in range(len(N)):

            R = (N[i] - 2*Ng) % mpi_sizes[i]
            normal_size = (N[i] - 2*Ng) / mpi_sizes[i] + 2*Ng
            augmnt_size = normal_size + 1
            thisdm_size = augmnt_size if mpi_coord[i]<R else normal_size

            for j in range(mpi_coord[i]):
                global_start[i] += augmnt_size-2*Ng if j<R else normal_size-2*Ng

            local_shape.append(thisdm_size)

            X0.append(x0[i] + dx[i] *  global_start[i])
            X1.append(x0[i] + dx[i] * (global_start[i] + thisdm_size-2*Ng))

        self.x0, self.x1 = X0, X1
        self.Ng = Ng
        self.cart = cart
        self.mpi_coord = mpi_coord
        self.mpi_sizes = mpi_sizes


    def synchronize(self, A):

        from mpi4py.MPI import COMM_WORLD, Compute_dims
        Ng = self.Ng

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


    def set_BC(self, A, BC=None):

        self.synchronize(A)
        if BC is None: return

        L_BCs = BC.L_wall[len(A.shape)-2]
        R_BCs = BC.R_wall[len(A.shape)-2]

        for i,BC in enumerate(L_BCs):
            if self.mpi_coord[i] == 0:
                BC(A, self.Ng)

        for i,BC in enumerate(R_BCs):
            if self.mpi_coord[i] == self.mpi_sizes[i]-1:
                BC(A, self.Ng)
