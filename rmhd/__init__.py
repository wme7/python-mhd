

from ctypes import CDLL, Structure, c_double, c_int
from numpy import float64, array, zeros, roll
from numpy.ctypeslib import ndpointer, c_intp
from sys import path
from os.path import abspath, dirname
from os import popen


_lib_home = dirname(abspath(popen('find . -name *.so').readline()))
path.append(abspath(_lib_home))
_type_dict = {'int': c_int, 'double': c_double}
_lib1 = CDLL(_lib_home+'/librmhd-1.so')
_lib2 = CDLL(_lib_home+'/librmhd-2.so')


import testbench
import riemann
import visual


class LibraryState(Structure):

    decl = """
  int cons_to_prim_iter;
  int cons_to_prim_use_estimate;
  int cons_to_prim_verbose;

  double max_lambda;
  double adiabatic_gamma;
  double plm_theta;

  int mode_reconstruct;
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

        self.mode_quartic_solver       = 0
        self.mode_reconstruct          = 1


        for k in kwargs:
            try: self.__getattribute__(k)
            except AttributeError: raise
            self.__setattr__(k,kwargs[k])



const_array = ndpointer(dtype=float64, flags=('C_CONTIGUOUS'))
write_array = ndpointer(dtype=float64, flags=('C_CONTIGUOUS', 'WRITEABLE'))

for lib in [_lib1, _lib2]:

    lib.initialize        .argtypes = [const_array, c_int]
    lib.finalize          .argtypes = [ ]
    lib.dUdt_1d           .argtypes = [const_array, write_array]
    lib.Fiph              .argtypes = [const_array, write_array]
    lib.prim_to_cons_array.argtypes = [const_array, write_array]
    lib.cons_to_prim_array.argtypes = [const_array, write_array]

    lib.initialize        .restype = c_int
    lib.finalize          .restype = c_int
    lib.dUdt_1d           .restype = c_int
    lib.Fiph              .restype = c_int
    lib.prim_to_cons_array.restype = c_int
    lib.cons_to_prim_array.restype = c_int

    lib.set_state.argtypes = [ LibraryState ]
    lib.get_state.restype  =   LibraryState
