

class SimpleCartesianDomain():

    def __init__(self, N=(256,), x0=(0.0,), x1=(1.0,)):
        assert len(N) is len(x0) is len(x1)
        self.N = N
        self.dx = [float(l1-l0)/n for n,l0,l1 in zip(N,x0,x1)]
        self.x0 = x0
        self.x1 = x1
        self.is_distributed = False
        self.rank = 0

    def set_BC(self, A, Ng, BC):
        L_BCs = BC.L_wall[len(A.shape)-1]
        R_BCs = BC.R_wall[len(A.shape)-1]

        for i,BC in enumerate(L_BCs): BC(A, Ng)
        for i,BC in enumerate(R_BCs): BC(A, Ng)
