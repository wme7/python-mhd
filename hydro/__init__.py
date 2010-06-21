

class HydrodynamicsSolver:

    def __init__(self, scheme='midpoint', solver='hll'):
        from ctypes import CDLL, POINTER, CFUNCTYPE, Structure, c_double, c_int
        from numpy import float64, int32
        from numpy.ctypeslib import ndpointer
        from sys import path
        from os.path import abspath, dirname
        from os import popen

        dbl_arr = ndpointer(dtype=float64, flags=('C_CONTIGUOUS', 'WRITEABLE'))
        dbl_vec = ndpointer(dtype=float64, flags=('C_CONTIGUOUS'))
        int_arr = ndpointer(dtype=int32  , flags=('C_CONTIGUOUS', 'WRITEABLE'))
        int_vec = ndpointer(dtype=int32  , flags=('C_CONTIGUOUS'))

        lib_home = dirname(abspath(popen('find . -name *.so').readline()))
        clib = CDLL(lib_home+'/'+self.libname+'.so')

        self.schemes = ['fwd_euler', 'midpoint', 'RK3', 'ctu_hancock']
        self.ghost_cells = dict(zip(self.schemes,(2,4,6,4)))
        self.advance = { }

        assert scheme in self.schemes
        assert solver in self.solvers

        for sname in self.schemes:
            self.advance[sname] = getattr(clib, 'advance_state_'+sname)
            self.advance[sname].argtypes = [dbl_arr, c_double]

        clib.integrate_init.argtypes = [int_vec, dbl_vec, c_int]
        clib.integrate_free.argtypes = [ ]
        clib.get_failure_mask.argtypes = [int_arr]
        clib.set_riemann_solver.argtypes = [c_int]

        self.clib = clib
        self.scheme = scheme
        self.NumGhostCells = self.ghost_cells[self.scheme]
        self.solver = solver

    def __del__(self):
        self.clib.integrate_free()

    def new_problem(self, domain):
        from numpy import array, float64, int32

        Ng = self.ghost_cells[self.scheme]
        Nq = self.NumComponents
        num_dims = len(domain.N)

        N = [Ng ] + [n+2*Ng for n in domain.N]
        L = [0.0] + [r-l for l,r in zip(domain.x0, domain.x1)]

        while (len(N)<4):
            N.append(1)
            L.append(0.0)

        self.clib.integrate_init(array(N,int32), array(L,float64), Nq, num_dims)
        self.clib.set_riemann_solver(self.solvers[self.solver])
        self.N = N
        self.L = L
        self.num_dims = num_dims

    def advance_state(self, P, dt):
        return self.advance[self.scheme](P, dt)

    def get_failure_mask(self):
        from numpy import zeros, int32
        M = zeros(self.N[1:self.num_dims+1], dtype=int32)
        self.clib.get_failure_mask(M)
        return M


class ScalarEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'scalar'
        self.NumComponents = 1
        self.solvers = {'hll': 0}
        HydrodynamicsSolver.__init__(self, **kwargs)


class EulersEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'euler'
        self.NumComponents = 5
        self.solvers = {'hll': 0}
        HydrodynamicsSolver.__init__(self, **kwargs)


class SRHDEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'srhd'
        self.NumComponents = 5
        self.solvers = {'hll': 0}
        HydrodynamicsSolver.__init__(self, **kwargs)


class RMHDEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'rmhd'
        self.NumComponents = 8
        self.solvers = {'hll': 0, 'hllc': 1}
        HydrodynamicsSolver.__init__(self, **kwargs)
