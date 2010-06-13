


class HydrodynamicsSolver:

    def __init__(self):

        raise NotImplementedError


    def _loadlib_(self, **kwargs):

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
        lib = CDLL(lib_home+'/'+self.libname+'.so')

        self.schemes = ['fwd_euler', 'midpoint', 'RK3', 'ctu_hancock']
        self.ghost_cells = dict(zip(self.schemes,(2,4,6,4)))
        self.advance = { }

        for sname in self.schemes:
            self.advance[sname] = lib.__getattr__('advance_state_'+sname)
            self.advance[sname].argtypes = [dbl_arr, c_double]

        lib.integrate_init.argtypes = [ int_vec, dbl_vec, c_int ]

        self._clib = lib
        self.kwargs = kwargs


    def new_problem(self):

        from numpy import zeros, ones, int32

        self.scheme = self.kwargs['scheme']
        assert self.scheme in self.schemes

        N = self.kwargs['N']
        L = self.kwargs['L']

        N_grid, L_grid = ones(4, dtype=int32), zeros(4)
        N_grid[0] = self.ghost_cells[self.scheme]
        num_dims = len(N)

        N_grid[1:1+num_dims] = N
        L_grid[1:1+num_dims] = L

        self.N = N
        self.L = L

        self._clib.integrate_init(N_grid, L_grid, self.NumComponents, num_dims)

    def advance_state(self, P, dt):
        self.advance[self.scheme](P, dt)


    def get_Ng(self):
        return self.ghost_cells[self.scheme]



class EulersEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):

        self.libname = 'euler'
        self.NumComponents = 5
        self._loadlib_(**kwargs)


class SRHDEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):

        self.libname = 'srhd'
        self.NumComponents = 5
        self._loadlib_(**kwargs)


class RMHDEquationsSolver(HydrodynamicsSolver):

    def __init__(self, **kwargs):

        self.libname = 'rmhd'
        self.NumComponents = 8
        self._loadlib_(**kwargs)
