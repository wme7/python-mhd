

class DistributedDomain():

    def __init__(self, N=(256,), x0=(0.0,), x1=(1.0,)):
        assert len(N) is len(x0) is len(x1)
        try:
            from mpi4py.MPI import COMM_WORLD, Compute_dims
        except ImportError:
            print "Error! DistributedDomain requires the mpi4py package."
            exit()

        mpi_sizes = Compute_dims(COMM_WORLD.size, len(N))
        cart = COMM_WORLD.Create_cart(mpi_sizes, periods=[True for n in N])
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
        self.global_start = global_start
        self.global_shape = global_shape
        self.is_distributed = True

    def synchronize(self, A, Ng):
        from mpi4py.MPI import COMM_WORLD

        if len(A.shape) == 2:
            L,R = self.cart.Shift(0,1)
            A[:+Ng] = self.cart.Sendrecv(A[-2*Ng:-Ng], dest=R, source=L)
            A[-Ng:] = self.cart.Sendrecv(A[+Ng:+2*Ng], dest=L, source=R)

        if len(A.shape) == 3:
            L,R = self.cart.Shift(0,1)
            A[:+Ng,:] = self.cart.Sendrecv(A[-2*Ng:-Ng,:], dest=R, source=L)
            A[-Ng:,:] = self.cart.Sendrecv(A[+Ng:+2*Ng,:], dest=L, source=R)

            L,R = self.cart.Shift(1,1)
            A[:,:+Ng] = self.cart.Sendrecv(A[:,-2*Ng:-Ng], dest=R, source=L)
            A[:,-Ng:] = self.cart.Sendrecv(A[:,+Ng:+2*Ng], dest=L, source=R)

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

        L_BCs = BC.L_wall[len(A.shape)-1]
        R_BCs = BC.R_wall[len(A.shape)-1]

        for i,BC in enumerate(L_BCs):
            if self.mpi_coord[i] == 0:
                BC(A, Ng)

        for i,BC in enumerate(R_BCs):
            if self.mpi_coord[i] == self.mpi_sizes[i]-1:
                BC(A, Ng)

    def dump(self, P, base='dump'):
        from pickle import dump
        pickle_stream = open(base+'_%03d.pk' % self.rank, 'w')
        tile = {'data': P, 'global_start': self.global_start}
        dump(tile, pickle_stream)
        pickle_stream.close()
        self.cart.Barrier()

    def read(self, base='dump'):
        from os import listdir, system
        from pickle import load
        from numpy import zeros

        fnames = [x for x in listdir('.') if x.endswith('.pk')]
        tiles = [load(open(fn)) for fn in fnames]
        ndims = len(tiles[0]['data'].shape)-1
        Nq = tiles[0]['data'].shape[-1]

        N = tuple(self.global_shape)
        P = zeros(N+(Nq,))

        for i,t in enumerate(tiles):
            pi = t['data']
            m0 = t['global_start']
            m1 = [g+s for g,s in zip(t['global_start'], pi.shape)]
            if len(N) == 1:
                i0, = m0
                i1, = m1
                P[i0:i1] = pi

            elif len(N) == 2:
                i0,j0 = m0
                i1,j1 = m1
                P[i0:i1,j0:j1] = pi

            elif len(N) == 3:
                i0,j0,k0 = m0
                i1,j1,k1 = m1
                P[i0:i1,j0:j1,k0:k1] = pi
        return P
