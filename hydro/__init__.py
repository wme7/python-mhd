

class HydrodynamicsSolver:

    def __init__(self, scheme='midpoint'):
        from ctypes import CDLL, POINTER, Structure, c_double, c_int
        from numpy import float64, int32
        from numpy.ctypeslib import ndpointer
        from sys import path
        from os.path import abspath, dirname
        from os import popen

        dbl_arr = ndpointer(dtype=float64, flags=('C_CONTIGUOUS', 'WRITEABLE'))
        dbl_vec = ndpointer(dtype=float64, flags=('C_CONTIGUOUS'))
        int_vec = ndpointer(dtype=int32  , flags=('C_CONTIGUOUS'))

        lib_home = dirname(abspath(popen('find . -name *.so').readline()))
        clib = CDLL(lib_home+'/'+self.libname+'.so')

        self.schemes = ['fwd_euler', 'midpoint', 'RK3', 'ctu_hancock']
        self.ghost_cells = dict(zip(self.schemes,(2,4,6,4)))
        self.advance = { }

        assert scheme in self.schemes

        for sname in self.schemes:
            self.advance[sname] = clib.__getattr__('advance_state_'+sname)
            self.advance[sname].argtypes = [dbl_arr, c_double]

        clib.integrate_init.argtypes = [ int_vec, dbl_vec, c_int ]

        self.clib = clib
        self.scheme = scheme
        self.NumGhostCells = self.ghost_cells[self.scheme]


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
        self.N = N
        self.L = L

    def advance_state(self, P, dt):
        self.advance[self.scheme](P, dt)


class ScalarEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'scalar'
        self.NumComponents = 1
        HydrodynamicsSolver.__init__(self, **kwargs)


class EulersEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'euler'
        self.NumComponents = 5
        HydrodynamicsSolver.__init__(self, **kwargs)


class SRHDEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'srhd'
        self.NumComponents = 5
        HydrodynamicsSolver.__init__(self, **kwargs)


class RMHDEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):
        self.libname = 'rmhd'
        self.NumComponents = 8
        HydrodynamicsSolver.__init__(self, **kwargs)
