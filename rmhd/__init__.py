

from ctypes import CDLL, POINTER, Structure, c_double, c_int
from numpy import float64
from numpy.ctypeslib import ndpointer
from sys import path
from os.path import abspath, dirname
from os import popen


_lib_home = dirname(abspath(popen('find . -name *.so').readline()))
path.append(abspath(_lib_home))
_type_dict = {'int': c_int, 'double': c_double}
_lib = CDLL(_lib_home+'/librmhd.so')



import testbench
import visual
import driver

try: import riemann
except ImportError: pass


class LibraryState(Structure):

    decl = """
  int cons_to_prim_iter;
  int cons_to_prim_use_estimate;
  int cons_to_prim_verbose;

  double max_lambda;
  double adiabatic_gamma;
  double plm_theta;

  int mode_riemann_solver;
  int mode_reconstruct;
  int mode_slope_limiter;
  int mode_quartic_solver;"""

    _fields_ = [(x.split()[1], _type_dict[x.split()[0]])
                for x in decl.split(';') if x.split()]


    def __init__(self, **kwargs):

        self.cons_to_prim_iter         = 0
        self.cons_to_prim_use_estimate = 0
        self.cons_to_prim_verbose      = 0

        self.max_lambda                = 0.0
        self.adiabatic_gamma           = 1.4
        self.plm_theta                 = 2.0

        self.mode_riemann_solver       = 0 # HLL approximate solver
        self.mode_reconstruct          = 2 # PLM on 4-velocity
        self.mode_slope_limiter        = 0 # Minmod limiter
        self.mode_quartic_solver       = 0 # Exact solver

        for k in kwargs:
            try: self.__getattribute__(k)
            except AttributeError: raise
            self.__setattr__(k,kwargs[k])



const_array = ndpointer(dtype=float64, flags=('C_CONTIGUOUS'))
write_array = ndpointer(dtype=float64, flags=('C_CONTIGUOUS', 'WRITEABLE'))

_lib.initialize        .argtypes = [const_array] + [c_int]*3 + [c_double]*3 + [c_int]
_lib.finalize          .argtypes = [ ]
_lib.dUdt_1d           .argtypes = [const_array, write_array]
_lib.dUdt_2d           .argtypes = [const_array, write_array]
_lib.dUdt_3d           .argtypes = [const_array, write_array]
_lib.Fiph              .argtypes = [const_array, write_array]
_lib.prim_to_cons_array.argtypes = [const_array, write_array, c_int]
_lib.cons_to_prim_array.argtypes = [const_array, write_array, c_int]
_lib.prim_to_cons_point.argtypes = [const_array, write_array]
_lib.cons_to_prim_point.argtypes = [const_array, write_array]


_lib.hll_flux          .argtypes = [const_array]*2 + [write_array]*2 + [c_double]
_lib.hllc_flux         .argtypes = [const_array]*2 + [write_array]*2 + [c_double]
_lib.set_dimension     .argtypes = [c_int]

_lib.set_state.argtypes = [ LibraryState ]
_lib.get_state.restype  =   LibraryState

_lib.get_failed_state.argtypes = [write_array]*2


# Access to the library's internal quartic solver
_lib.new_QuarticEquation.argtypes    = [c_double]*5
_lib.solve_quartic_equation.argtypes = [POINTER(c_double)]*4 + [POINTER(c_int)]*2
