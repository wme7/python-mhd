

class CartesianDomain:

    def __init__(self, N=(128,), x0=(0.0,), x1=(1.0,)):

        assert len(N) is len(x0) is len(x1)

        self.x0 = x0
        self.x1 = x1
        self.dx = [(r-l)/n for n,l,r in zip(N,x0,x1)]


class OutflowBoundaryConditions:

    def
